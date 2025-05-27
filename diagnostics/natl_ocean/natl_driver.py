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

temp_var = [case['thetao_var'] for case in case_list.values()][0]
#uvel_var = [case['UO_var'] for case in case_list.values()][0]
#vvel_var = [case['VVEL_var'] for case in case_list.values()][0]
hfds_var = [case['hfds_var'] for case in case_list.values()][0]
salt_var = [case['so_var'] for case in case_list.values()][0]

for case in case_list.values():
    if 'vsf_var' in case:
        vsf_var = [case['vsf_var'] for case in case_list.values()][0]
        wfo_mod = False
    elif 'wfo_var' in case:
        wfo_var = [case['wfo_var'] for case in case_list.values()][0]
        wfo_mod = True
    else: 'uhoh'

areacello_var = [case['areacello_var'] for case in case_list.values()][0]

#loading coords; this is currently writtne for multicase mode, but ignoring the other possible cases, 
#because there's only one
#change later to work for single case mode?
time_coord = [case['time_coord'] for case in case_list.values()][0]
lon_coord = [case['lon_coord'] for case in case_list.values()][0]
lat_coord = [case['lat_coord'] for case in case_list.values()][0]
lev_coord = [case['lev_coord'] for case in case_list.values()][0] 


#method 1 for loading in variables: use the catalog
cat_def_file = case_info['CATALOG_FILE']
print(cat_def_file)

cat = intake.open_esm_datastore(cat_def_file)

print('This is the catalog', cat)

##TEMP
temp_subset = cat.search(variable_id=temp_var, frequency="month")
temp_dict = temp_subset.to_dataset_dict(
    xarray_open_kwargs={"decode_times": True, "use_cftime": True}
)


#method 2: load the file directly
#THETAO
input_path = os.environ["THETAO_FILE"]
print('THETAO_FILE is:', input_path)

# command to load the netcdf file
model_temp_dataset = xr.open_dataset(input_path)
print(model_temp_dataset)

#SALT
input_path = os.environ["SO_FILE"]
model_salt_dataset = xr.open_dataset(input_path)

#SHF
input_path = os.environ["HFDS_FILE"]
model_hfds_dataset = xr.open_dataset(input_path)

#TAREA
input_path = os.environ["AREACELLO_FILE"]
model_area_dataset = xr.open_dataset(input_path)


#taking time mean for figures
temp_tmean = model_temp_dataset[temp_var].isel({lev_coord:0}).mean(time_coord)
salt_tmean = model_salt_dataset[salt_var].isel({lev_coord:0}).mean(time_coord)
hfds_tmean = model_hfds_dataset[hfds_var].mean(time_coord)
area = model_area_dataset[areacello_var]



WORK_DIR = os.environ['WORK_DIR']
outmod_dir = os.path.join(WORK_DIR, "model")
outobs_dir = os.path.join(WORK_DIR, "obs")


#TEMP fig
plt.pcolormesh(temp_tmean)
plt.colorbar()
plt.savefig(outmod_dir+'/tmean_toplev_plot.png')

#SALT fig
f=plt.figure()
plt.pcolormesh(salt_tmean)
plt.colorbar()
plt.savefig(outmod_dir+'/smean_toplev_plot.png')

#SHF fig
f=plt.figure()
plt.pcolormesh(hfds_tmean)
plt.colorbar()
plt.savefig(outmod_dir+'/hfds_plot.png')

#TAREA fig
f=plt.figure()
plt.pcolormesh(area)
plt.colorbar()
plt.savefig(outmod_dir+'/area_plot.png')



################## PART 1: NORTH ATLANTIC BIAS ASSESSMENT #######################
#LOAD IN T/S OBS AND OMIP DATASET 
obsdir = os.environ["OBS_DATA"]


print('this is the obsdir', obsdir)

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



