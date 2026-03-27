# ALMA 3D Fitting Guide for `dysmalpy` (Wrapper Workflow)

This guide is a practical, up-to-date workflow for fitting ALMA spectral-line cubes with `dysmalpy` using the fitting wrapper (`.params` + `dysmalpy_fit_single`).

It focuses on 3D cube fitting and common failure modes.

---

## 1. Prerequisites

1. Use a Python environment where `dysmalpy` imports from your intended checkout.
2. Confirm the import path before fitting:

```bash
python - <<'PY'
import dysmalpy
print(dysmalpy.__file__)
PY
```

3. Prepare these files:
- Data cube (`fdata_cube`): 3D FITS, shape `(nspec, ny, nx)`.
- Error cube (`fdata_err`): same shape as data cube.
- Optional mask files (`fdata_mask`, `fdata_mask_sky`, `fdata_mask_spec`) if not using auto-mask.

---

## 2. Critical Spectral-Axis Rules

For robust 3D wrapper fitting:

1. Prefer a velocity cube in FITS with:
- `CTYPE3 = VRAD` (or velocity-equivalent)
- `CUNIT3 = m/s` (recommended for wrapper auto-conversion to km/s)

2. Ensure channel direction is increasing in velocity for model stability:
- `CDELT3 > 0` is strongly recommended.

3. In `.params`, do not use `spec_orig_type=freq` for wrapper 3D loading.
- Wrapper conversion logic supports `spec_orig_type = wave` or `velocity`.
- If your original cube is frequency, convert to velocity before fitting.

4. Keep trimming consistent with the line center:
- `spec_vel_trim` should include the actual emission line channels.

---

## 3. Recommended `.params` Template (3D ALMA)

Use this as a base and edit values for your target:

```text
# ******************************* OBJECT INFO **********************************
galID,           ALMA_Target_1
z,               2.530

# ****************************** DATA INFO *************************************
datadir,         /path/to/data/
fdata_cube,      alma_line_cube_vel.fits
fdata_err,       alma_line_errcube.fits
# fdata_mask,    alma_mask_cube.fits

# Auto-mask (recommended if no hand mask)
auto_gen_3D_mask, True
auto_gen_mask_apply_skymask_first, True
auto_gen_mask_snr_thresh_pixel, None
auto_gen_mask_sig_segmap_thresh, 3.0
auto_gen_mask_npix_segmap_min, 5
auto_gen_mask_sky_var_thresh, 2.
auto_gen_mask_snr_int_flux_thresh, 3.

# Optional crop/trim in data cube coordinates
spec_vel_trim,   -500. 500.
# spatial_crop_trim,  40 100 40 100
# xcenter, 70.0
# ycenter, 70.0

data_inst_corr,  False

# ***************************** OUTPUT *****************************************
outdir,          /path/to/output/ALMA_Target_1/
overwrite,       True
do_plotting,     True

# ***************************** OBSERVATION SETUP ******************************
pixscale,        0.05
fov_npix,        60
spec_type,       velocity
spec_start,      -1000.
spec_step,       10.
nspec,           201

smoothing_type,  None
smoothing_npix,  1

slit_width,      None
slit_pa,         0.
integrate_cube,  False

use_lsf,         True
sig_inst_res,    20.0

psf_type,        Gaussian
psf_fwhm_major,  0.15
psf_fwhm_minor,  0.12
psf_PA,          45.0
psf_beta,        -99.

# **************************** SETUP MODEL *************************************
components_list,       disk+bulge   const_disp_prof   geometry   zheight_gaus  halo
light_components_list, disk

adiabatic_contract,  False
pressure_support,    True
noord_flat,          True
oversample,          1
oversize,            1
zcalc_truncate,      True
n_wholepix_z_min,    3

# Disk+bulge
total_mass,          10.5
bt,                  0.0
r_eff_disk,          2.5
n_disk,              1.0
invq_disk,           5.0
n_bulge,             4.0
invq_bulge,          1.0
r_eff_bulge,         1.0

total_mass_fixed,    False
r_eff_disk_fixed,    False
bt_fixed,            True
n_disk_fixed,        True
r_eff_bulge_fixed,   True
n_bulge_fixed,       True

total_mass_bounds,   9.0 12.0
r_eff_disk_bounds,   0.5 10.0

# Halo
halo_profile_type,   NFW
mvirial,             11.5
halo_conc,           5.0
fdm,                 0.5
mvirial_fixed,       False
halo_conc_fixed,     True
fdm_fixed,           False
mvirial_bounds,      10.0 13.0
fdm_bounds,          0.0 1.0
fdm_tied,            True
mvirial_tied,        False

# Dispersion
sigma0,              40.0
sigma0_fixed,        False
sigma0_bounds,       10.0 150.0

# Z-height
sigmaz,              0.9
sigmaz_fixed,        False
sigmaz_bounds,       0.1 5.0
zheight_tied,        True

# Geometry
inc,                 45.0
pa,                  90.0
xshift,              0.0
yshift,              0.0
vel_shift,           0.0

inc_fixed,           False
pa_fixed,            False
xshift_fixed,        False
yshift_fixed,        False
vel_shift_fixed,     False

inc_bounds,          10.0 80.0
pa_bounds,           -180. 180.0
xshift_bounds,       -5.0 5.0
yshift_bounds,       -5.0 5.0
vel_shift_bounds,    -200. 200.

# **************************** FITTING SETTINGS ********************************
fit_method,      mpfit
maxiter,         200
```

Notes:
- `fov_npix/spec_start/spec_step/nspec` are initial instrument settings. For 3D fits, wrapper may reset them to match loaded data after trim/crop.
- Elliptical beam keys (`psf_fwhm_major/minor`, `psf_PA`) are supported.

---

## 4. Run the Fit (Correct Wrapper API)

Use this minimal script:

```python
import os
from dysmalpy.fitting_wrappers import utils_io
from dysmalpy.fitting_wrappers import dysmalpy_fit_single

param_filename = "/path/to/alma_3D.params"
plot_type = "png"
overwrite = True

params = utils_io.read_fitting_params(fname=param_filename)
os.makedirs(params["outdir"], exist_ok=True)

# Official wrapper entry point
dysmalpy_fit_single.dysmalpy_fit_single(
    param_filename=param_filename,
    plot_type=plot_type,
    overwrite=overwrite,
)
```

Do not use old/nonexistent calls like `dysmalpy_fit_single.run_fitting(...)`.

---

## 5. Fast Pre-Fit Sanity Check

Before long fitting runs, test model generation once:

```python
from dysmalpy.fitting_wrappers import utils_io

param_filename = "/path/to/alma_3D.params"
params = utils_io.read_fitting_params(fname=param_filename)
gal, output_options = utils_io.setup_single_galaxy(params=params)

obs = gal.observations['OBS']
print("data shape:", obs.data.data.shape)
print("spec start/step/nspec:", obs.instrument.spec_start, obs.instrument.spec_step, obs.instrument.nspec)

# This should succeed without exception
gal.create_model_data()
print("model shape:", obs.model_data.data.shape)
```

If this fails, fix setup issues before launching MPFIT/MCMC.

---

## 6. Common Errors and Fixes

### `IndexError` in `_create_cube_mask` (tuple index out of range)

Typical cause in 3D ALMA wrapper runs:
- Spectral axis direction is descending (`spec_step < 0`), which can lead to empty model cube after convolution.

Fix:
1. Rewrite velocity cube so spectral axis increases (`CDELT3 > 0`).
2. Regenerate error cube and params.
3. Re-run the pre-fit sanity check.

### `PermissionError` opening `*_mpfit.log`

Cause:
- `outdir` not writable by current process.

Fix:
- Use writable `outdir` and ensure directory exists and permissions are correct.

### `ModuleNotFoundError: dysmalpy`

Cause:
- Wrong Python environment.

Fix:
- Activate the intended environment and confirm import path before fitting.

### Fitting silently reuses old results

Cause:
- Existing result files and `overwrite=False`.

Fix:
- Set `overwrite, True` in params or pass `overwrite=True` in wrapper call.

---

## 7. Output Interpretation

After successful fit, inspect `outdir`:

1. `*_fit_report.txt` and `*_bestfit_results.dat`
- Best-fit values for mass, geometry, dispersion, halo terms, and goodness metrics.

2. Best-fit plots (`png`/`pdf`)
- Data-model comparisons for extracted products and channel/spaxel diagnostics.

3. Pickles (`*_model.pickle`, `*_{fit_method}_results.pickle`)
- Re-load for post-analysis.

---

## 8. ALMA-Specific Practical Tips

1. Start with conservative masking and modest crop around the source to stabilize first fits.
2. Keep `oversample=1`, `oversize=1` initially; increase only after a stable baseline fit.
3. Lock weakly constrained parameters (`bt`, some halo terms, `zheight`) early, then relax.
4. Verify beam (`BMAJ/BMIN/BPA`) and spectral resolution (`sig_inst_res`) from headers/observing setup, not memory.
5. If line center is uncertain, run a few fits with different `vel_shift` priors or center assumptions.

---

## 9. Minimal Checklist

1. Cube/error same shape.
2. Velocity axis valid and increasing.
3. `spec_vel_trim` contains line emission.
4. `outdir` writable.
5. Pre-fit `gal.create_model_data()` succeeds.
6. Then run wrapper fit.
