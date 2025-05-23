# North Atlantic Diagnostics POD Driver
# Last update: 2/11/2025
#   Version & Contact info
# 
#   Here you should describe who contributed to the diagnostic, and who should be
#   contacted for further information:
# 
#   - Version/revision information: version 1 (5/06/2020)
#   - PI (name, affiliation, email)
#   - Developer/point of contact (name, affiliation, email)
#   - Other contributors
# 
#   Open source copyright agreement
# 
#   The MDTF framework is distributed under the LGPLv3 license (see LICENSE.txt). 
#   Unless you've distirbuted your script elsewhere, you don't need to change this.
# 
#   Functionality
# 
#   In this section you should summarize the stages of the calculations your 
#   diagnostic performs, and how they translate to the individual source code files 
#   provided in your submission. This will, e.g., let maintainers fixing a bug or 
#   people with questions about how your code works know where to look.
# 
#   Required programming language and libraries
# 
#   In this section you should summarize the programming languages and third-party 
#   libraries used by your diagnostic. You also provide this information in the 
#   ``settings.jsonc`` file, but here you can give helpful comments to human 
#   maintainers (eg, "We need at least version 1.5 of this library because we call
#   this function.")
#   Required model output variables
# 
#   In this section you should describe each variable in the input data your 
#   diagnostic uses. You also need to provide this in the ``settings.jsonc`` file, 
#   but here you should go into detail on the assumptions your diagnostic makes 
#   about the structure of the data.
# 
#   References
# 
#   Here you should cite the journal articles providing the scientific basis for 
#   your diagnostic.
# 

######## IMPORT PACKAGES #############
import xarray as xr
import matplotlib.pyplot as plt
import os
import yaml
import intake
import numpy as np

#### LOAD IN MODEL VARIABLES NEEDED FOR ALL PARTS  ####
print("I MADE IT!")

case_env_file = os.environ["case_env_file"]
assert os.path.isfile(case_env_file), f"case environment file not found"
with open(case_env_file, 'r') as stream:
    try:
        case_info = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

cat_def_file = case_info['CATALOG_FILE']
case_list = case_info['CASE_LIST']
# all cases share variable names and dimension coords in this example, so just get first result for each
print(case_list)

sst_var = [case['SST_var'] for case in case_list.values()][0]
temp_var = [case['TEMP_var'] for case in case_list.values()][0]
uvel_var = [case['UVEL_var'] for case in case_list.values()][0]
vvel_var = [case['VVEL_var'] for case in case_list.values()][0]
shf_var = [case['SHF_var'] for case in case_list.values()][0]
wfo_var = [case['SFWF_var'] for case in case_list.values()][0]
tarea_var = [case['TAREA_var'] for case in case_list.values()][0]


time_coord = [case['time_coord'] for case in case_list.values()][0]
#lon_coord = [case['lon_coord'] for case in case_list.values()][0]
#lat_coord = [case['lat_coord'] for case in case_list.values()][0]
z_t_coord = [case['z_t_coord'] for case in case_list.values()][0] 

cat_def_file = case_info['CATALOG_FILE']
print(cat_def_file)

cat = intake.open_esm_datastore(cat_def_file)

sst_subset = cat.search(variable_id=sst_var, frequency="month")

sst_dict = sst_subset.to_dataset_dict(
    xarray_open_kwargs={"decode_times": True, "use_cftime": True}
)


input_path = os.environ["SST_FILE"]
print('SST_FILE is:', input_path)

# command to load the netcdf file
model_sst_dataset = xr.open_dataset(input_path)
print(model_sst_dataset)

time = model_sst_dataset[time_coord]
#lat = model_sst_dataset[lat_coord]
#lon = model_sst_dataset[lon_coord]

sst_tmean = model_sst_dataset[sst_var].isel({z_t_coord:0}).mean(time_coord)

WORK_DIR = os.environ['WORK_DIR']
outmod_dir = os.path.join(WORK_DIR, "model")
outobs_dir = os.path.join(WORK_DIR, "obs")

plt.pcolormesh(sst_tmean)
plt.colorbar()
plt.savefig(outmod_dir+'/example_model_plot.png')

################## PART 1: NORTH ATLANTIC BIAS ASSESSMENT #######################
#LOAD IN T/S OBS AND OMIP DATASET 
obsdir = os.environ["OBS_DATA"]

f=plt.figure()
plt.pcolormesh(sst_tmean)
plt.colorbar()
plt.savefig(outobs_dir+'/example_obs_plot.png')


#ds_en4 = xr.open_dataset(obsdir+'en4.nc')

#CALCULATIONS

#PLOTS 
#POSTAGE STAMPS OF T AND S DIAGRAM
#TAYLOR DIAGRAM?
#OTHER PLOT?

#SAVE FIGS -> HTML

################## PART 2: AMOC IN SIGMA COORDS #############
#LOAD IN OMIP AMOC(SIGMA)

#ds_omip_amoc = xr.open_dataset(obsdir+'amoc.nc')


#CALCULATIONS
#CALLING MOC FUNCTIONS FROM PY SCRIPT

#PLOTS
#AMOC IN SIGMA
#AMOC IN Z (If TIME)
#AMOC at 45 Line Plot

#SAVE FIGS -> HTML

#################  PART 3: SURFACE-FORCED WATER MASS TRANSFORMATION ##########
#LOAD IN WMT BENCHMARKS


#CALCULATIONS

#PLOTS
#WMT BY REGION
#WMT(45N+) WITH AMOC(SIGMA)

#SAVE FIGS -> HTML

################# PART 4: SYNTHESIS #########################################



