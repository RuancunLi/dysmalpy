# Skill: Performing a 3D Dynamical Fit on an ALMA Line Cube using `dysmalpy`

This guide provides a step-by-step, "agent-learnable" tutorial on how to use `dysmalpy` to perform a 3D dynamical fit on an ALMA spectral line cube.

`dysmalpy` uses forward modeling, meaning it generates a high-resolution 3D intrinsic cube, applies instrumental effects (beam smearing, spectral resolution), and compares the output directly to the data cube voxel-by-voxel.

---

## Step 1: Data Preparation

To run a 3D fit in `dysmalpy`, you need the following FITS files ready:
1.  **Data Cube (`fdata_cube`):** The 3D ALMA line cube (e.g., [C II] or CO emission). The spectral axis should preferably be in velocity space (`km/s`). If it's in frequency or wavelength, `dysmalpy` can handle it, provided you specify the rest line frequency/wavelength.
2.  **Noise/Error Cube (`fdata_err`):** A 3D FITS cube representing the $1\sigma$ uncertainty per voxel.
3.  **Mask (Optional but Recommended, `fdata_mask`):** A 2D or 3D FITS mask where `1` indicates spaxels/voxels to include in the fit, and `0` indicates pixels to ignore. If omitted, `dysmalpy` can auto-generate a mask.

---

## Step 2: Creating the Parameter File (`alma_3D.params`)

`dysmalpy` uses parameter files (`.params`) to configure the galaxy models, instrumental effects, and fitting routine. Create a file named `alma_3D.params`.

Here is an agent-learnable template customized for an ALMA observation:

```text
# ******************************* OBJECT INFO **********************************
galID,           ALMA_Target_1
z,               2.53             # Redshift of the source

# ****************************** DATA INFO *************************************
datadir,         /path/to/data/   # Optional: path to your FITS files
fdata_cube,      alma_line_cube.fits
fdata_err,       alma_noise_cube.fits
fdata_mask,      alma_mask.fits   # Optional (can use auto_gen_3D_mask, True)

# Spectral setup if cube is not strictly optical velocity
spec_orig_type,  wave             # 'wave', 'freq', or 'velocity'
spec_line_rest,  157.74           # Rest wavelength for [CII]
spec_line_rest_unit, um           # 'um', 'angstrom', 'GHz'

# Masking behavior
auto_gen_3D_mask, True            # Generate a mask automatically if fdata_mask is missing
auto_gen_mask_sig_segmap_thresh, 3.0

# ***************************** OUTPUT *****************************************
outdir,          /path/to/output/ALMA_Target_1_Out/

# ***************************** OBSERVATION SETUP ******************************
pixscale,        0.05             # Pixel scale in arcsec/pixel
use_lsf,         True             # True/False: apply a Line Spread Function
sig_inst_res,    20.0             # Instrumental dispersion (spectral resolution in km/s)

# ALMA Synthesized Beam (PSF) Setup
psf_type,        Gaussian         # ALMA beams are well-approximated by 2D Gaussians
psf_fwhm_major,  0.15             # Beam Major axis FWHM (arcsec)
psf_fwhm_minor,  0.12             # Beam Minor axis FWHM (arcsec)
psf_PA,          45.0             # Beam Position Angle (Degrees East of North)

# **************************** SETUP MODEL *************************************
components_list,       disk+bulge   const_disp_prof   geometry   zheight_gaus  halo
light_components_list, disk

adiabatic_contract,  False
pressure_support,    True         # Apply asymmetric drift correction
noord_flat,          True         # Apply Noordermeer thickness flattening

# ------------- DISK + BULGE -------------
# Initial Values (Log M_sun, kpc)
total_mass,          10.5         # Total mass
bt,                  0.0          # Bulge-to-total ratio (0.0 = pure disk)
r_eff_disk,          2.5          # Disk effective radius (kpc)
n_disk,              1.0          # Sersic index (1 = exponential disk)
invq_disk,           5.0          # R_eff / z_height ratio (thickness)

# Fixed / Free Parameters
total_mass_fixed,    False
r_eff_disk_fixed,    False
bt_fixed,            True
n_disk_fixed,        True

# Prior Bounds
total_mass_bounds,   9.0  12.0
r_eff_disk_bounds,   0.5  10.0

# ------------- DARK MATTER HALO -------------
halo_profile_type,   NFW
mvirial,             11.5         # Halo virial mass log(M_sun)
halo_conc,           5.0          # Concentration
fdm,                 0.5          # DM fraction
mvirial_fixed,       False
halo_conc_fixed,     True
fdm_fixed,           False
mvirial_bounds,      10.0 13.0
fdm_bounds,          0.0  1.0

# Tie the parameters so mvirial and fdm are mathematically linked
fdm_tied,            True
mvirial_tied,        False

# ------------- DISPERSION PROFILE -------------
sigma0,              40.0         # Constant intrinsic dispersion (km/s)
sigma0_fixed,        False
sigma0_bounds,       10.0 150.0

# ------------- GEOMETRY -------------
inc,                 45.0         # Inclination (Degrees)
pa,                  90.0         # Kinematic Position Angle (Degrees)
xshift,              0.0          # Center shift in x (pixels)
yshift,              0.0          # Center shift in y (pixels)
vel_shift,           0.0          # Systemic velocity shift (km/s)

inc_fixed,           False
pa_fixed,            False
xshift_fixed,        False
yshift_fixed,        False
vel_shift_fixed,     False

inc_bounds,          10.0  80.0
pa_bounds,           -180. 180.0
xshift_bounds,       -5.0  5.0
yshift_bounds,       -5.0  5.0
vel_shift_bounds,    -200. 200.

# **************************** FITTING SETTINGS ********************************
fit_method,      mpfit            # 'mpfit' for Levenberg-Marquardt, 'mcmc', or 'nested'
do_plotting,     True             # Generate standard diagnostic plots
maxiter,         200              # Max iterations for mpfit
```

---

## Step 3: Running the Fit in Python

Once the parameter file is defined, you can execute the fit using `dysmalpy`'s built-in fitting wrappers. This python script loads the FITS data, configures the mask, and triggers the `MPFIT` (or MCMC) optimization.

Create a script `run_alma_fit.py`:

```python
import os
import numpy as np
from dysmalpy.fitting_wrappers import dysmalpy_fit_single
from dysmalpy.fitting_wrappers import utils_io
from dysmalpy import plotting

# 1. Path to your param file
param_filename = 'alma_3D.params'

# 2. Read the parameters
params = utils_io.read_fitting_params(fname=param_filename)

# 3. Load the galaxy object (loads FITS cube, err, sets up instrument/models)
gal = utils_io.load_galaxy(params=params, datadir=params.get('datadir', None), skip_mask=True)

# 4. Generate/Load 3D Mask
# This uses parameters like `auto_gen_mask_sig_segmap_thresh` from the .params file
mask, mask_dict = utils_io.generate_3D_mask(obs=gal.observations['OBS'], params=params)

# Optional: Plot the generated mask for visual verification
plotting.plot_3D_data_automask_info(gal.observations['OBS'], mask_dict)
import matplotlib.pyplot as plt
plt.savefig(os.path.join(params['outdir'], 'mask_diagnostic.png'))

# 5. Run the Dynamical Fit
# This step does the heavy lifting: comparing the forward-modeled cube to ALMA data
results = dysmalpy_fit_single.run_fitting(gal, params=params, nproc=4)

# 6. Save Results and Plots
# Save the text report
f_report = os.path.join(params['outdir'], f"{params['galID']}_fit_report.txt")
results.results_report(gal=gal, filename=f_report)

# Generate diagnostic plots (moment maps, PV diagrams, residuals)
results.plot_results(gal)
```

## Step 4: Interpreting Results

After the script finishes, check the output directory (`/path/to/output/ALMA_Target_1_Out/`).

1. **`ALMA_Target_1_fit_report.txt` / `*_bestfit_results.dat`:**
   These files contain the best-fit parameters, including:
   - `total_mass` (log M_sun of the stellar+gas disk)
   - `r_eff_disk` (Disk effective radius)
   - `sigma0` (Intrinsic velocity dispersion)
   - `inc`, `pa` (Geometry)
   - `mvirial`, `fdm` (Dark matter properties)
   - `redchisq` (Reduced Chi-Square goodness-of-fit)

2. **Diagnostic Plots:**
   `dysmalpy` automatically extracts 2D maps from the fitted 3D model and compares them to the 2D maps extracted from the ALMA data. Look for:
   - `*_data_model_residual_flux.png`
   - `*_data_model_residual_velocity.png` (Rotation Curve / Velocity map)
   - `*_data_model_residual_dispersion.png` (Velocity Dispersion map)
   - Position-Velocity (PV) diagrams along the major axis.

These outputs will indicate whether the forward-modeled 3D cube accurately captured the kinematics of your ALMA observation.