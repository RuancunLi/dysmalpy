# How `dysmalpy` Calculates Moment 1 and Moment 2

`dysmalpy` is a 3D forward modeling tool used to fit and simulate the kinematics of galaxies. When providing or fitting 2D kinematic maps—like a velocity map (Moment 1) and a velocity dispersion map (Moment 2)—`dysmalpy` works by first generating an intrinsic high-resolution 3D data cube based on the specified mass, light, geometry, and dispersion models. It then applies observational effects (such as spatial/spectral resolution and instrument beam smearing), and finally extracts the 2D Moment 1 and Moment 2 maps from this degraded 3D model cube.

This document details the exact mathematical formulation and python code used by `dysmalpy` to achieve this.

## 1. Generating the 3D Intrinsic Model Cube

Before any moment maps are calculated, `dysmalpy` builds a 3D data cube of the galaxy's line emission. This process takes place primarily in the `ModelSet.simulate_cube()` function (`dysmalpy/models/model_set.py`).

For each pixel (or spatial "spaxel") defined in the 3D model:
1.  **Coordinates:** It transforms coordinates from the sky frame $(x_{sky}, y_{sky}, z_{sky})$ to the galaxy's internal coordinate frame $(x_{gal}, y_{gal}, z_{gal})$ using the model geometry (inclination, position angle, shifts).
2.  **Velocity ($v$):** The circular rotational velocity ($v_{rot}$) is calculated based on the mass enclosed within a radius $r_{gal} = \sqrt{x_{gal}^2 + y_{gal}^2}$. This intrinsic velocity is then projected along the line of sight (LOS) to generate the observed LOS velocity for that voxel, $v_{LOS}$.
3.  **Dispersion ($\sigma$):** The intrinsic velocity dispersion $\sigma$ is evaluated based on the given dispersion profile evaluated at the spatial position.
4.  **Flux ($f$):** The emission flux $f$ is calculated from the defined 3D light distribution model, multiplied by any extinction and dimming factors.

### Populating the Cube (The Spectral Profile)

Once $v_{LOS}$, $\sigma$, and $f$ are determined for each spatial voxel $(x, y, z)$, `dysmalpy` populates the 3D data cube with a Gaussian spectral profile for each spatial point. This is highly optimized using Cython in `dysmalpy/models/cutils.pyx`.

For a given spectral velocity array $v_{spec}$, the flux density $F$ at a specific channel $s$, and spatial pixel $(x,y)$, is populated as:

$$ F(v_{spec}[s], x, y) = \sum_{z} \frac{f(x, y, z)}{\sqrt{2\pi}\sigma(x, y, z)} \exp \left( -\frac{(v_{spec}[s] - v_{LOS}(x, y, z))^2}{2\sigma(x, y, z)^2} \right) $$

```cython
# From dysmalpy/models/cutils.pyx

# f = flux[z, y, x]
# v = vel[z, y, x]
# sig = sigma[z, y, x]

amp = f / sqrt(2.0 * pi * sig)
for s in range(vspec.shape[0]):
    result[s, y, x] += amp * exp(-0.5 * ((vspec[s] - v) / sig) **2)
```

After the 3D model cube is generated, instrumental effects such as Beam/PSF convolution and Line Spread Function (LSF) convolution are applied using the `Instrument.convolve()` method.

---

## 2. Extracting Moment 1 (Velocity) and Moment 2 (Dispersion)

To get 2D kinematic maps, `dysmalpy` extracts them directly from the instrument-convolved 3D cube. This takes place in the `Observation.create_single_obs_model_data()` function (`dysmalpy/observation.py`), inside the conditional block `elif (ndim_final == 2):`.

There are two primary ways `dysmalpy` performs this extraction, dictated by the `instrument.moment` flag:

1.  **Gaussian Fitting** (`extrac_type = 'gauss'`) -> *The default behavior.*
2.  **Standard Moments** (`extrac_type = 'moment'`) -> *True moments of the distribution.*

### Method A: Gaussian Fitting (`extrac_type = 'gauss'`)

When extracting via Gaussian fits, `dysmalpy` fits a 1D Gaussian to the spectrum of *each individual spatial spaxel* $(x, y)$ in the 3D cube.

For a spectrum $S(v)$ at pixel $(x, y)$, it finds the best-fit parameters $(A, \mu, \sigma)$ for the equation:

$$ S(v) \approx A \cdot \exp \left(-\frac{(v - \mu)^2}{2\sigma^2} \right) $$

-   **Moment 1 (Velocity):** The fitted mean $\mu$.
-   **Moment 2 (Dispersion):** The fitted standard deviation $\sigma$.

**Implementation Details:**
- `dysmalpy` uses the standard true moments (Method B below) as the **initial guesses** for the Gaussian fit.
- It can utilize a fast, multi-threaded C++ fitter (`LeastChiSquares1D`), or fallback to `scipy.optimize.leastsq` via the `gaus_fit_sp_opt_leastsq` utility.

```python
# From dysmalpy/observation.py
# (Initial guesses using true moments)
mom0 = self.model_cube.data.moment0().to(u.km/u.s).value
mom1 = self.model_cube.data.moment1().to(u.km/u.s).value
mom2 = self.model_cube.data.linewidth_sigma().to(u.km/u.s).value

# ...
# Fallback Python fitting method:
for i in range(mom0.shape[0]):
    for j in range(mom0.shape[1]):
        best_fit = gaus_fit_sp_opt_leastsq(
            self.model_cube.data.spectral_axis.to(u.km/u.s).value,
            self.model_cube.data.unmasked_data[:,i,j].value,
            mom0[i,j], mom1[i,j], mom2[i,j]
        )
        flux[i,j] = best_fit[0] * np.sqrt(2 * np.pi) * best_fit[2]
        vel[i,j] = best_fit[1]       # Moment 1
        disp[i,j] = best_fit[2]      # Moment 2
```

### Method B: True Moments (`extrac_type = 'moment'`)

If configured to use true moments, `dysmalpy` leverages the `spectral_cube` Python package to perform direct numerical integration over the spectral axis $v$.

For a spectrum $S_i$ at velocity bin $v_i$ with bin width $\Delta v$:

**Moment 0 (Integrated Flux):**
$$ M_0 = \sum S_i \Delta v $$

**Moment 1 (Intensity-Weighted Mean Velocity):**
$$ v_{mom1} = \frac{\sum v_i S_i \Delta v}{M_0} $$

**Moment 2 (Velocity Dispersion / Line-width Sigma):**
This is the square root of the intensity-weighted variance:
$$ \sigma_{mom2} = \sqrt{\frac{\sum (v_i - v_{mom1})^2 S_i \Delta v}{M_0}} $$

```python
# From dysmalpy/observation.py
elif extrac_type == 'moment':
    # Leveraging astropy/spectral_cube moment calculations
    flux = self.model_cube.data.moment0().to(u.km/u.s).value
    vel = self.model_cube.data.moment1().to(u.km/u.s).value
    disp = self.model_cube.data.linewidth_sigma().to(u.km/u.s).value
```

## 3. Configuring Mass Models (Stellar, Gas, Dark Matter)

In `dysmalpy`, the total mass model is an aggregate of several `MassModel` components combined within a `ModelSet`. These components explicitly define the enclosed mass profile $M_{enc}(r)$, which determines the kinematics. The user configures these models to represent the stellar, gas, and dark matter content.

### A. Baryonic Mass (Stellar & Gas)

Baryonic mass is usually modeled using a **Sersic** or an **Exponential Disk** profile. The `Sersic` profile defines the surface mass density $\Sigma(r)$ as:

$$ \Sigma(r) = \Sigma_e \exp \left\{ -b_n \left[ \left( \frac{r}{r_{\mathrm{eff}}} \right)^{1/n} - 1 \right] \right\} $$

Where:
*   $r_{\mathrm{eff}}$ is the effective (half-mass) radius.
*   $n$ is the Sersic index (where $n=1$ corresponds to an exponential disk, typical for gas or stellar disks, and $n=4$ corresponds to a de Vaucouleurs profile typical for bulges).
*   $b_n$ is a constant determined numerically such that $r_{\mathrm{eff}}$ contains half the total mass.

The 2D projected enclosed mass for a Sersic profile is calculated as an incomplete Gamma function:

$$ M_{enc}(r) = M_{total} \cdot \frac{\gamma\left(2n, b_n (r/r_{\mathrm{eff}})^{1/n}\right)}{\Gamma(2n)} $$

**Configuration Example:**
```python
from dysmalpy.models import Sersic, DiskBulge

# Add a combined stellar disk and bulge
stellar_comp = DiskBulge(total_mass=10.5, r_eff_disk=3.0, n_disk=1.0,
                         r_eff_bulge=1.0, n_bulge=4.0, bt=0.3) # bt = bulge-to-total ratio

# Alternatively, add a gas disk component
gas_comp = Sersic(total_mass=9.5, r_eff=4.0, n=1.0, baryon_type='gas')
```
*Note:* `dysmalpy` uses lookup tables derived from *Noordermeer (2008)* to compute the 3D potential for geometrically "thick" Sersic disks, avoiding the infinite-cylinder approximation.

### B. Dark Matter Halo

Dark matter halos can be modeled using several profiles, most commonly the **NFW (Navarro, Frenk, & White)** profile or the **Burkert** (cored) profile.

**1. NFW Profile:**
The density profile $\rho(r)$ is given by:

$$ \rho(r) = \frac{\rho_0}{(r/r_s)(1 + r/r_s)^2} $$

Where:
*   $r_s = r_{\mathrm{vir}}/c$ is the scale radius ($r_{\mathrm{vir}}$ is the virial radius, $c$ is the concentration).
*   $\rho_0$ is the normalization density tied to the total virial mass $M_{\mathrm{vir}}$.

The spherically enclosed dark matter mass $M_{enc, \mathrm{DM}}(r)$ is:

$$ M_{enc, \mathrm{DM}}(r) = 4\pi \rho_0 r_s^3 \left[ \ln\left(1 + \frac{r}{r_s}\right) - \frac{r}{r+r_s} \right] $$

**2. Burkert Profile:**
A cored profile, often favored in some dark matter dominated galaxies:

$$ \rho(r) = \frac{\rho_0}{(1 + r/r_B)(1 + (r/r_B)^2)} $$

Where $r_B$ is the core radius.

**Configuration Example:**
```python
from dysmalpy.models import NFW

# Add an NFW Dark Matter Halo
halo = NFW(mvirial=12.0, conc=10.0) # log(M_vir) = 12.0 M_sun, c = 10
```

---

## 4. Dynamical Parametrization of Intrinsic Velocity ($v$) and Dispersion ($\sigma$)

With the mass models configured, the 3D model cube requires an intrinsic rotation velocity $v_{rot}(r)$ and intrinsic velocity dispersion $\sigma(r)$ at every point. These properties are derived dynamically from the `ModelSet`.

### A. Intrinsic Dispersion ($\sigma$)

In `dysmalpy`, the intrinsic velocity dispersion of the gas is modeled independently via a `DispersionProfile` component. The default is a constant dispersion across the disk (`DispersionConst`):

$$ \sigma(r) = \sigma_0 $$

This $\sigma_0$ is a free parameter in the model (e.g., `DispersionConst.sigma0`) and represents the isotropic turbulence/dispersion of the gas.

### B. Circular Velocity ($v_{circ}$)

First, the total *circular* velocity $v_{circ}(r)$ is computed from the gravitational potential of all mass components (baryonic disk, bulge, dark matter halo, etc.).

For each component $i$, the circular velocity squared $v_{circ, i}^2(r)$ is derived from the enclosed mass $M_{enc, i}(r)$:
$$ v_{circ, i}^2(r) = \frac{G M_{enc, i}(r)}{r} $$

The total circular velocity squared is the sum of the squares:
$$ v_{circ}^2(r) = \sum_i v_{circ, i}^2(r) $$

*Note:* `dysmalpy` also has an option (`KinematicOptions.adiabatic_contract`) to numerically alter $v_{circ}^2$ to account for adiabatic contraction of the dark matter halo due to the baryonic mass.

### C. Rotation Velocity ($v_{rot}$) with Asymmetric Drift

The actual ordered rotation velocity of the gas, $v_{rot}(r)$, is lower than the theoretical circular velocity $v_{circ}(r)$ due to pressure support (also called *asymmetric drift*). Because the gas has internal turbulence (represented by $\sigma_0$), the pressure gradient partially balances gravity.

In `dysmalpy`, this is handled by `KinematicOptions.apply_pressure_support()`. The relationship is:

$$ v_{rot}^2(r) = v_{circ}^2(r) - v_{asym}^2(r) $$

`dysmalpy` provides three different parametrizations (`pressure_support_type`) for the asymmetric drift $v_{asym}^2(r)$ (following Burkert et al. 2010):

**Type 1 (Default): Exponential Disk Approximation**
Assuming the gas follows an exponential disk profile:
$$ v_{asym}^2(r) = 3.36 \left( \frac{r}{R_e} \right) \sigma_0^2 $$
*(Where $R_e$ is the effective radius).*

**Type 2: Sersic Index Dependence**
Taking into account the specific Sersic index $n$ of the baryonic disk:
$$ v_{asym}^2(r) = 2 \left( \frac{b_n}{n} \right) \left( \frac{r}{R_e} \right)^{1/n} \sigma_0^2 $$
*(Where $b_n \approx 1.9992n - 0.3271$ is the Sersic constant).*

**Type 3: Direct Pressure Gradient**
Calculating the drift explicitly from the local log-density gradient of the gas:
$$ v_{asym}^2(r) = -\left( \frac{d\ln\rho_{gas}(r)}{d\ln r} \right) \sigma_0^2 $$

In `dysmalpy` code, this is executed as:

```python
# From dysmalpy/models/kinematic_options.py
if pressure_support_type == 1:
    vel_asymm_drift_sq = 3.36 * (r / Re) * sigma ** 2
elif pressure_support_type == 2:
    bn = scipy.special.gammaincinv(2. * n, 0.5)
    vel_asymm_drift_sq = 2. * (bn/n) * np.power((r/Re), 1./n) * sigma**2
elif pressure_support_type == 3:
    dlnrhogas_dlnr = model.get_dlnrhogas_dlnr(r)
    vel_asymm_drift_sq = - dlnrhogas_dlnr * sigma**2

# Final rotation velocity applied to the cube
v_rot = np.sqrt(v_circ_sq - vel_asymm_drift_sq)
```

## Summary for `GalfitS` Integration

Because `dysmalpy` explicitly generates the full 3D model cube, computing the moment maps is straightforward.

If integrating into `GalfitS`:
1.  **Build the 3D model:** Integrate over the line of sight ($z$) for each sky pixel $(x, y)$ to produce a velocity spectrum (adding up the flux contributions mapped to their respective line-of-sight velocities and intrinsic dispersions).
2.  **Convolve:** Apply spatial and spectral PSFs.
3.  **Extract:** Collapse the spectral axis. Gaussian fitting is generally preferred in generic forward modeling as it better mirrors the process observers use to extract moment maps from real IFU (Integral Field Unit) data. True moments are faster but are heavily sensitive to noise and asymmetric line profiles.