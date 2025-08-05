# North Atlantic Diagnostics POD Driver
# Last update: 6/30/2025
#   Version & Contact info
#   - Version/revision information: version 1 (7/30/2025)
#   - PIs: Liz Maroon, University of Wisconsin, emaroon@wisc.edu
#          Steve Yeager, NSF National Center for Atmospheric Resarch, yeager@ucar.edu
#   - Developer/point of contact: Liz Maroon, University of Wisconsin, emaroon@wisc.edu
#   - Other contributors: Taydra Low, Brendan Myers, Teagan King
# 
#   Open source copyright agreement
# 
#   The MDTF framework is distributed under the LGPLv3 license (see LICENSE.txt). 
#   Unless you've distirbuted your script elsewhere, you don't need to change this.
# 
#   Functionality
# 
#   1. Preprocessing & Derived Variable Computation:
#        Prepare the model output and compute derived physical variables using POD_utils.py.
#        Key Functions:
#          compute_sigma0(thetao, so): Calculates potential density anomaly (σ₀).
#          compute_mld(sigma0): Computes Mixed Layer Depth based on a density threshold.
#          compute_zavg(ds, var): Calculates thickness-weighted mean of a variable over the upper 200 m.
#          check_depth_units(ds): Ensures depth units are in meters.
#
#   2. Regridding:
#        Interpolate all model and observational fields to a consistent 1×1 lat-lon grid using POD_utils.py.
#        Key Functions:
#          regrid(ds, method='bilinear'): Uses xESMF for bilinear regridding, with NaN handling and global domain definition.
#
#   3. Time Averaging (Climatology)
#        Convert the time series into a monthly climatology over 1989–2018 using user-controlled driver script
#
#   4. Bias Calculation Against Observations
#        Compute model bias relative to observations and calculate error statistics using POD_utils.py.
#        Key Functions:
#          plot_preproc(...): Extracts spatial and temporal slices, sets up metadata.
#          error_stats(...): Computes area-weighted RMSE (map region) and mean bias (focus region).
#
#   5. Plotting
#        Generate spatial and statistical visualizations of model performance using POD_utils.py.
#        Key Functions:
#          SpatialPlot_climo_bias(...): Produces two-panel maps for each variable (bias + bias rank).
#            SpatialBias_panel(...): Shows spatial bias and regional stats.
#            SpatialRank_panel(...): Shows where the target model ranks in bias among OMIP models.
#          ScatterPlot_Error(...): Scatter plot comparing regional bias statistics between two variables.
#            Uses Scatter_panel(...) to visualize model-model comparisons.
#
#   6. Optional: Multi-Cycle OMIP Time Handling
#        Detect and restructure repeating OMIP forcing cycles.
#        Key Functions:
#          forcing_cycles(expid, nt): Identifies number and span of OMIP cycles.
#          reorg_by_cycle(ds, nt, ncyc, yearrange): Restructures the dataset to expose cycles on a new dimension.
#
#   Required programming language and libraries
# 
#   The North Atlantic Ocean Diagnostic recommends python (3.10 or later) because we
#   use xarray. Xarray, matplotlib, os, yaml, intake, numpy, xesmf, xskillscore,
#   scipy, gsw_xarray, numba, cftime, and cartopy are also required.
# 
#   Required model output variables
#     thetao    Potential temperature   degrees Celsius  3D: time × depth × lat × lon
#     so        Salinity                PSU              3D: time × depth × lat × lon
#     lev       Depth level             m or cm          2D
#     lon, lat  Longitude and latitude  degrees          1D or 2D grid coordinates
#     time      Time dimension          datetime         1D
#
#   References  # TODO
# 
#   Here you should cite the journal articles providing the scientific basis for 
#   your diagnostic.

# Import Packages
import xarray as xr
import os
import yaml
import POD_utils

print('Starting North Atlantic Ocean POD')

# User Settings #########################################################

# Plot Lat/Lon Region
plot_region = [360-90, 360-0, 20, 80]

# Focus Lat/Lon Region
focus_region = [360-48, 360-30, 38, 53]

# Analysis Month (Note: set to 13 for 'Annual')
month = 13

# Output
savefig = True

# LOAD IN MODEL VARIABLES NEEDED FOR ALL PARTS  ##########################
case_env_file = os.environ["case_env_file"]
assert os.path.isfile(case_env_file), f"case environment file not found"
with open(case_env_file, 'r') as stream:
    try:
        case_info = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

cat_def_file = case_info['CATALOG_FILE']
case_list = case_info['CASE_LIST']
model_name = list(case_list.keys())[0]
start_year = case_list['CESM2_historical_r1i1p1f1']['startdate'].split('-')[0]
end_year = case_list['CESM2_historical_r1i1p1f1']['enddate'].split('-')[0]

# all cases share variable names and dimension coords in this example, so just get first result for each
volcello_var = [case['volcello_var'] for case in case_list.values()][0]
areacello_var = [case['areacello_var'] for case in case_list.values()][0]
temp_var = [case['thetao_var'] for case in case_list.values()][0]
hfds_var = [case['hfds_var'] for case in case_list.values()][0]
salt_var = [case['so_var'] for case in case_list.values()][0]

for case in case_list.values():
    if 'vsf_var' in case:
        vsf_var = [case['vsf_var'] for case in case_list.values()][0]
        wfo_mod = False
    elif 'wfo_var' in case:
        wfo_var = [case['wfo_var'] for case in case_list.values()][0]
        wfo_mod = True
    else:
        print('vsf_var or wfo_var not found in case')

# Load the files ------------------------------------------------------
# ThetaO
model_temp_dataset = xr.open_dataset(os.environ["THETAO_FILE"])

# Salt
model_salt_dataset = xr.open_dataset(os.environ["SO_FILE"])

# SHF
model_hfds_dataset = xr.open_dataset(os.environ["HFDS_FILE"])

# TArea
model_area_dataset = xr.open_dataset(os.environ["AREACELLO_FILE"])

# Volume
model_vol_dataset = xr.open_dataset(os.environ["VOLCELLO_FILE"])

vol = model_vol_dataset[volcello_var]
area = model_area_dataset[areacello_var]
dz = vol/area
if "lev" not in dz.coords:
    dz = dz.assign_coords(lev=model_vol_dataset["lev"])
dz = dz.assign_coords(lev=dz.lev / 100.0)  # Convert to meters

# ---------------------------------------------------------------------

# set directories
WORK_DIR = os.environ['WORK_DIR']
outmod_dir = os.path.join(WORK_DIR, "model")
outobs_dir = os.path.join(WORK_DIR, "obs")

# PART 1: NORTH ATLANTIC BIAS ASSESSMENT ######################################

print('At Part 1: North Atlantic Bias Assessment')
# Data Ingest from Catalogue
ds_target = model_temp_dataset
ds_target['so'] = model_salt_dataset['so']

ds_target = POD_utils.preprocess_coords(ds_target)

# LOAD IN T/S OBS AND OMIP DATASET ---------------------------------------------
obsdir = os.environ["OBS_DATA"]
# omip_dir = os.environ["OMIP_DATA"]

# Open OMIP data  # TODO: this should probably just load omip above but file not ingested yet!
omip_file = '/glade/work/brendanmy/S_Yeager/Sub2Sub/data_archive/POD_data/omip2.cycle1.1989_2018.mld_sic_t200_s200_sigma200.nc'
# omip_file = omip_dir+'omip2.cycle1.1989_2018.mld_sic_t200_s200_sigma200.nc'
ds_model = xr.open_dataset(omip_file).isel(OMIP=0).load()
#ds_model = xr.open_dataset(omip_file).load()

# Open Obs # TODO: this should probably just use the obsdir above but file not ingested yet!
obs_path = '/glade/campaign/cgd/ccr/yeager/Sub2Sub/POD_data/obs_1x1.nc'
ds_obs = xr.open_dataset(obs_path).load()
# ds_obs = xr.open_dataset(obsdir+'obs_1x1.nc').load()

# Time Subselection:
ds_target = ds_target.sel(time=slice(start_year, end_year))

# PERFORM CALCULATIONS --------------------------------------------------------
# Compute Sigma0 and MLD 
ds_target['sigma0'] = POD_utils.compute_sigma0(ds_target['thetao'], ds_target['so'])
ds_target['mld'] = POD_utils.compute_mld(ds_target['sigma0'])

# Compute Depth-average Fields (hard-wired for 200m-depth average)
zavg_var_list = ['thetao', 'so', 'sigma0']
for var in zavg_var_list:
    if not isinstance(dz, xr.DataArray):
        raise TypeError(f"Expected dz to be a DataArray, got {type(dz)}")
    ds_target = POD_utils.compute_zavg(ds_target, var, dz)

# Drop 3D fields
ds_target = ds_target.drop_vars(['thetao','so','sigma0'])

# Regrid
ds_target = POD_utils.regrid(ds_target, method='bilinear')

# Compute climatology
ds_target = ds_target.groupby('time.month').mean('time', keep_attrs=True)
ds_target = ds_target.assign_coords({'model': model_name})

# CREATE PLOTS ---------------------------------------------------------------------
ds_t200 = POD_utils.SpatialPlot_climo_bias(ds_target, ds_model, ds_obs, 'thetao_zavg', region=plot_region, focus_region=focus_region, month=month, save=savefig, savedir=outmod_dir)
ds_s200 = POD_utils.SpatialPlot_climo_bias(ds_target, ds_model, ds_obs, 'so_zavg', region=plot_region, focus_region=focus_region, month=month, save=savefig, savedir=outmod_dir)
POD_utils.SpatialPlot_climo_bias(ds_target, ds_model, ds_obs, 'sigma0_zavg', region=plot_region, focus_region=focus_region, month=month, save=savefig, savedir=outmod_dir)
POD_utils.SpatialPlot_climo_bias(ds_target, ds_model, ds_obs, 'mld', region=plot_region, focus_region=focus_region, month=month, save=savefig, savedir=outmod_dir)
POD_utils.ScatterPlot_Error(ds_t200, 'thetao_zavg_bias', ds_s200, 'so_zavg_bias', model_name, save=savefig, savedir=outmod_dir)

# Wrap-up by closing datasets that have been opened and informing user of successful completion
model_temp_dataset.close()
model_salt_dataset.close()
model_hfds_dataset.close()
model_area_dataset.close()
ds_target.close()
print('North Atlantic Ocean POD Part 1: North Atlantic Bias Assessment finished successfully!')

# PART 2: AMOC IN SIGMA COORDS #####################################################
# LOAD IN OMIP AMOC(SIGMA)
# ds_omip_amoc = xr.open_dataset(obsdir+'amoc.nc')

# PERFORM CALCULATIONS -----------------------------------------------------------
# CALLING MOC FUNCTIONS FROM PY SCRIPT

# CREATE PLOTS -------------------------------------------------------------------
# AMOC IN SIGMA
# AMOC IN Z (If TIME)
# AMOC at 45 Line Plot

# SAVE FIGS -> HTML
# Wrap-up by closing datasets that have been opened and informing user of successful completion
# print('North Atlantic Ocean POD Part 2: AMOC finished successfully!')

# PART 3: SURFACE-FORCED WATER MASS TRANSFORMATION ###############################
# LOAD IN WMT BENCHMARKS

# PERFORM CALCULATIONS -----------------------------------------------------------

# CREATE PLOTS -------------------------------------------------------------------
# WMT BY REGION
# WMT(45N+) WITH AMOC(SIGMA)

# SAVE FIGS -> HTML
# print('North Atlantic Ocean POD Part 3: Surface-Forced Water Mass Transformation finished successfully!')

# PART 4: SYNTHESIS ##############################################################
# Wrap-up by closing datasets that have been opened and informing user of successful completion

print("North Atlantic Ocean POD finished successfully!")
