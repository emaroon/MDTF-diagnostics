# Import Analysis Tools
import numpy as np
import xarray as xr
import xesmf as xe
from scipy import stats
import gsw_xarray as gsw
from numba import guvectorize
import cftime

# Import Plotting Tools
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from matplotlib.colors import BoundaryNorm

# Import Colors
import matplotlib.colors as mcolors

# Warnings are hidden with the below code. Comment out if you want warnings
import warnings
warnings.filterwarnings('ignore')

##### PROCESSING  #####
def preprocess_coords(ds):
    rename_coords_dict = {
        'd2': 'bnds',
        'axis_nbounds': 'bnds',
        'olevel': 'lev',
        'olevel_bounds': 'lev_bnds',
        'lev_bounds': 'lev_bnds',
        'time_bounds': 'time_bnds',
        'bounds_lon': 'lon_bnds',
        'bounds_lat': 'lat_bnds',
        'bounds_nav_lon': 'lon_bnds',
        'bounds_nav_lat': 'lat_bnds',
        'nav_lon': 'lon',
        'nav_lat': 'lat',
        'latitude': 'lat',
        'longitude': 'lon',
        'longitude_bnds':'lon_bnds',
        'latitude_bnds':'lat_bnds',
        'i': 'x',
        'j': 'y',
        'nlat':'y',
        'nlon':'x',
        'z_l': 'lev',
        'z_i': 'lev_bnds',
        'xh': 'lon',
        'yh': 'lat',
    }
    
    # rename when possible
    for k, v in rename_coords_dict.items():
        try:
            ds = ds.rename({k: v})
        except:
            pass
            
    ## May not be needed in general (just for CESM?)
    ds = check_depth_units(ds)

    return ds

def compute_sigma0(da_t,da_s):
    """
    Compute potential density anomaly (sigma0) from conservative temperature and salinity.

    Parameters:
    da_t : xarray.DataArray
        Conservative temperature in degrees Celsius ('degC').
    da_s : xarray.DataArray
        Absolute salinity in g/kg. Accepts 'psu' or '0.001' as equivalent.

    Returns:
    xarray.DataArray
        Sigma0 density anomaly (kg/m³ - 1000).

    Raises:
    Exception
        If input units are not recognized.
    """
    if (da_t.attrs['units']!='degC'):
            raise Exception("Check units of thetao!")
    if (da_s.attrs['units']!='g/kg'):
        if (da_s.attrs['units']=='psu') or (da_s.attrs['units']=='0.001'):
            da_s = da_s.assign_attrs({'units':'g/kg'})
        else:
            raise Exception("Check units of so!")
    da_sig0 = gsw.density.sigma0(SA=da_s,CT=da_t)
    return da_sig0


def compute_mld(da_sig0, dsig=0.03, rdep=10.0):
    """
    Compute mixed layer depth (MLD) based on a density threshold criterion.

    Parameters:
    da_sig0 : xarray.DataArray
        Potential density (sigma0) with units 'kg/m^3'.
    dsig : float, optional
        Density increase threshold from reference depth (default: 0.03 kg/m³).
    rdep : float, optional
        Reference depth in meters (default: 10.0 m).

    Returns:
    xarray.DataArray
        Mixed layer depth in meters.

    Raises:
    Exception
        If input units are not 'kg/m^3'.
    """
    if (da_sig0.attrs['units']!='kg/m^3'):
        raise Exception("Check units of sigma0!")
    zdim = 'lev'
    da_sig0_ref = da_sig0.interp(lev=rdep)
    z=da_sig0[zdim]
    k_ref = int(np.min(np.argwhere(z.values >= rdep)))
    da_mld = xr.apply_ufunc(mld_gufunc, da_sig0, z, k_ref, da_sig0_ref, dsig, 
                            input_core_dims=[[zdim], [zdim], [], [], []], dask='parallelized')
    return da_mld


@guvectorize(['void(float64[:], float64[:], intp, float64,  float64,float64[:])'], '(n),(n),(),(),()->()', nopython=True)
def mld_gufunc(phi,  z, k_ref, phi_ref,  phi_step, res):
    #  function from Guillaume Serazin via Anne-Marie Treguier
    #  phi is an xarray of density sigma0
    #  z is the depth
    #  k_ref is the reference level (first index below the reference depth)
    #  phi_ref is the density at the reference level 
    #  phi_step=0.03 (default) is the density threshold to define the mixed layer depth
    #  res is the result = mixed layer depth    
    phi_mlb = phi_ref + phi_step
    if (np.isfinite(phi[k_ref])) :
        jkmax = np.max(np.argwhere(np.isfinite(phi)))
        #  find the indices of points below kref where tabprof and reference 
        #  differ by more than zdel (use absolute value)
        indices=np.argwhere((phi[k_ref:] >= phi_mlb))
        if (indices.size == 0) :
            #  no values differ by more than zdel: mldepth is depth of the last level
            res[0]=z[jkmax]
        else :
            k_mlb = np.min(indices)+k_ref   
            # Make a linear interpolation to find the approximate MLD
            # the logic is to find a point on a 
            delta_z = z[k_mlb - 1] - z[k_mlb]
            alpha = delta_z/(phi[k_mlb - 1] - phi[k_mlb])
            beta = z[k_mlb] - alpha * phi[k_mlb]
            res[0] = alpha * phi_mlb + beta
    else :
        res[0] = np.nan


def get_dlev(lev, lev_bnds, depth_limit=10000):
    """
    Compute vertical layer thickness from depth bounds, with optional depth clipping.

    Parameters:
    lev_bnds : xarray.DataArray
        Depth bounds with shape (..., 2), where the last or first dimension corresponds to layer bounds.
    depth_limit : float, optional
        Maximum depth for clipping bounds (default: 10000 m).

    Returns:
    xarray.DataArray
        Layer thickness (dlev) in meters.
    """
    lev_bnds_lim = lev_bnds.where(lev_bnds < depth_limit, depth_limit)
    if (np.size(lev_bnds_lim.dims)==1):
        dlev = xr.DataArray(lev_bnds_lim[1:].values - lev_bnds_lim[0:-1].values,coords={'lev':lev})
    elif (np.size(lev_bnds_lim.dims)==2):
        dlev = lev_bnds_lim.isel(bnds=1) - lev_bnds_lim.isel(bnds=0)
    elif (np.size(lev_bnds_lim.dims)==3):
        dlev = xr.DataArray(lev_bnds_lim[1:].values - lev_bnds_lim[0:-1].values) #,coords={'lev':lev})
    else:
        raise ValueError('ERROR: could not handle lev_bnds')
    return dlev

def check_depth_units(ds):
    # Convert lev if needed
    if 'lev' in ds.coords:
        lev = ds['lev']
        if lev.max() > 8000:
            lev_converted = lev / 100
            lev_converted.attrs.update(lev.attrs)
            lev_converted.attrs['units'] = 'm'
            ds = ds.assign_coords(lev=lev_converted)

    # Convert lev_bnds if needed
    if 'lev_bnds' in ds.variables:
        lev_bnds = ds['lev_bnds']
        if lev_bnds.max() > 8000:
            lev_bnds_converted = lev_bnds / 100
            lev_bnds_converted.attrs.update(lev_bnds.attrs)
            lev_bnds_converted.attrs['units'] = 'm'
            ds['lev_bnds'] = lev_bnds_converted

    return ds

def compute_zavg(ds, var, dz, depth=200):
    """
    Compute thickness-weighted mean over depth for field.
    """
    # Ensure lev is available
    lev = ds['lev']
    # Find valid levels shallower than depth
    valid_mask = lev <= depth
    valid_levs = lev.where(valid_mask, drop=True)
    if valid_levs.size == 0:
        print(f"Warning: No valid levels found for depth={depth}m in variable '{var}'")
        return ds

    # Slice the variable and weights to these levels
    data = ds[var].sel(lev=valid_levs)
    dz_sel = dz.sel(lev=valid_levs)
    dz_sel = dz_sel.rename({'nlat': 'y', 'nlon': 'x'})

    # Do the weighted mean
    newvar = f"{var}_zavg"

    # Add to dataset and annotate
    ds[newvar] = data.weighted(dz_sel).mean('lev', keep_attrs=True).astype('float32')
    ds[newvar].attrs['zavg'] = f'0-{depth}m'

    return ds


def regrid(ds, var=None, names=None, target_grid=None, dlon=1, dlat=1, method='bilinear'):
    """
    Regrid a dataset or variable to a regular lat-lon grid using xESMF.

    Parameters:
    ds : xarray.Dataset
        Input dataset containing 'lat' and 'lon' coordinates.
    var : str, optional
        Name of the variable to regrid. If None, regrids the entire dataset.
    names : unused
        Placeholder parameter (not currently used).
    target_grid : dict or xarray.Dataset, optional
        Target grid definition. If None, creates a global grid with given resolution.
    dlon : float, optional
        Longitudinal resolution for default grid (default: 1°).
    dlat : float, optional
        Latitudinal resolution for default grid (default: 1°).
    method : str, optional
        Regridding method (e.g., 'bilinear', 'nearest_s2d', 'conservative').

    Returns:
    xarray.DataArray or xarray.Dataset
        Regridded variable or dataset.

    Notes:
    - Filters invalid lat/lon values before regridding.
    - Uses periodic boundary conditions and ignores degenerate grids.
    """
    if target_grid is None:
        target_grid = xe.util.grid_global(dlon, dlat, cf=True, lon1=360)
    
    # ds = utils.regrid(self[k], var=vn, target=target_grid, method=method)
    dsgrid = ds[['lat', 'lon']].copy()

    ## some grids (GFDL,NCAR) have large fill values that should actually be NaN
    lat = dsgrid['lat']
    lon = dsgrid['lon']
    dsgrid['lat'] = lat.where(lat>-90).where(lat<90)
    dsgrid['lon'] = lon.where(lon>-360).where(lon<360)

    regridder = xe.Regridder(dsgrid, target_grid, method=method, periodic=True, ignore_degenerate=True)
    if var is None:
        rgd = regridder(ds, keep_attrs=True, skipna=True, na_thres=0.6)
    else:
        rgd = regridder(ds[var], keep_attrs=True, skipna=True, na_thres=0.6)

    return rgd


def forcing_cycles(expid,nt):
    ''' Function to determine the number and year range of forcing
    cycles included in this OMIP data array.'''
    maxcycles = 10
    if (expid=="omip1" or expid=="omip1-spunup"):
        # For OMIP1, we expect a 1948-2009 (62-year) forcing cycle = 744 mon
        # but some groups submitted a 1948-2007 (60-year) cycle = 720 mon
        # and MIROC submitted a 1958-2009 (52 year) cycle = 624 mon
        if (nt % 744 == 0):
            nyear = 62
            yearrange = (1948,2009)
        elif (nt % 720 == 0):
            nyear = 60
            yearrange = (1948,2007)
        elif (nt % 624 == 0):
            nyear = 52
            yearrange = (1958,2009)
        else:
            raise ValueError('ERROR: could not determine OMIP1 forcing cycle')
    elif (expid=="omip2" or expid=="omip2-spunup"):
        # For OMIP2, we expect a 1958-2018 (61-year) forcing cycle = 732 mon
        if (nt % 732 == 0):
            nyear = 61
            yearrange = (1958,2018)
        elif (nt % 624 == 0):
            nyear = 52
            yearrange = (1958,2009)
        else:
            raise ValueError('ERROR: could not determine OMIP2 forcing cycle')
    else:
        raise ValueError('ERROR: experiment_id not recognized')
    ntmon = nyear*12*np.arange(1,maxcycles+1,1)
    findcyc = np.where(nt==ntmon)
    ncyc = findcyc[0][0] + 1
    print('  found '+str(ncyc)+' forcing cyles spanning '+str(yearrange))
    return ncyc,yearrange


##### PLOTTING #####
def blue2red_cmap(n):
    """ combine two existing color maps to create a diverging color map with white in the middle
    n = the number of contour intervals
    """
    if (int(n/2) == n/2):
        # even number of contours
        nwhite=1
        nneg=n/2
        npos=n/2
    else:
        nwhite=2
        nneg = (n-1)/2
        npos = (n-1)/2
    colors1 = plt.cm.Blues_r(np.linspace(0,1, int(nneg)))
    colors2 = plt.cm.YlOrRd(np.linspace(0,1, int(npos)))
    colorsw = np.ones((nwhite,4))
    colors = np.vstack((colors1, colorsw, colors2))
    cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    return cmap

def get_varname(var):
    vardict = {'thetao_zavg':'T200','so_zavg':'S200','sigma0_zavg':r'${\sigma_0}$200','mld':'MLD'}
    return vardict[var]

def SpatialBias_panel(da, stats, focus_region, focus_model, ax):
    var = da.name
    var_name = get_varname(var)
    da_model = da.sel(model=focus_model)  #Create DataArray of the target model alone
    rmse_model = stats[var+'_rmse'].sel(model=focus_model).values
    bias_model = stats[var+'_bias'].sel(model=focus_model).values

    # Fontsizes
    font_title = 18
    font_label = 16
    font_tick = 16
    font_stat = 14

    # Settings
    proj = ccrs.PlateCarree()
    facecolor="grey"
        
    # set up contour levels and color map
    if var=='so_zavg':
        ci = 0.25; cmin = -2.5; cmax = 2.5
    elif var=='thetao_zavg':
        ci = 0.5; cmin = -8; cmax = 8
    elif var=='sigma0_zavg':
        ci = 0.1; cmin = -2; cmax = 2
    else:
        ci = 50; cmin = -500; cmax = 500
    nlevs = (cmax-cmin)/ci + 1
    levs = np.arange(cmin, cmax+ci, ci)
    cmap = cmap = blue2red_cmap(nlevs)
    cmap.set_over('magenta')
    cmap.set_under('cyan')
    norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)

    # Setup
    ax.set_aspect('auto')
    ax.set_facecolor(facecolor)
    ax.add_feature(cfeature.COASTLINE)

    # Plot
    cntr1 = ax.pcolormesh(da_model['lon'], da_model['lat'], da_model, shading='nearest', cmap=cmap, norm=norm, rasterized=True, transform=proj)
    cbar = plt.colorbar(cntr1, ax=ax, orientation='vertical', spacing='proportional', pad=.01)
    cbar.set_label(label=r'{}'.format(da.units), size=font_label, rotation=270, labelpad=5)
    cbar.ax.tick_params(labelsize=font_tick)

    # Show focus region
    r = focus_region
    ax.add_patch(mpatches.Rectangle(xy=[r[0], r[2]], width=abs(r[0]-r[1]), height=abs(r[3]-r[2]),
                                    facecolor=None, edgecolor='lime', fill=False, linewidth=2,
                                    transform=ccrs.PlateCarree()))

    gl = ax.gridlines(draw_labels=True, dms=True,  linewidth=0.5, color='k', alpha=0.8)
    gl.xlabel_style = {'fontsize': font_tick}
    gl.ylabel_style = {'fontsize': font_tick}
    gl.top_labels = False
    gl.right_labels = False
    ax.set_title('{}:\nClimatological {} Bias ({})'.format(focus_model,var_name,da.time_avg), fontsize=font_title)
    ax.set_title(r'rmse={:.1f}'.format(rmse_model) + '\n' + r'bias={:.1f}'.format(bias_model), fontsize=font_stat, loc='right')

def SpatialRank_panel(da, focus_region, focus_model, ax):
    var = da.name
    var_name = get_varname(var)
    da_model = da.sel(model=focus_model)  # Create DataArray of the target model alone

    # Fontsizes
    font_title = 18
    font_label = 16
    font_tick = 16

    # Settings
    proj = ccrs.PlateCarree()
    facecolor="grey"
    levs = np.arange(0, 110, 10)  # Colorbar Levels
    nlevs = levs.size
    cmap = cmap = blue2red_cmap(nlevs)
    norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)

    # Setup
    ax.set_aspect('auto')
    ax.set_facecolor(facecolor)
    ax.add_feature(cfeature.COASTLINE)

    # Process Rank Data
    rank_da = xr.apply_ufunc(stats.percentileofscore, # I think default is kind='rank' & nan_policy='propogate'
                             abs(da.where(da.model!=focus_model, drop=True)),
                             abs(da_model),
                             'rank',
                             'omit',
                             input_core_dims=[['model'], [], [], []],
                             vectorize=True,
                            )

    # Plot
    cntr1 = ax.pcolormesh(rank_da['lon'], rank_da['lat'], rank_da, shading='nearest', cmap=cmap, norm=norm, rasterized=True, transform=proj)
    cbar = plt.colorbar(cntr1, ax=ax, orientation='vertical', spacing='proportional', pad=0.01)
    cbar.set_label(label='%', size=font_label, rotation=270, labelpad=5)
    cbar.ax.tick_params(labelsize=font_tick)

    # Show focus region
    r = focus_region
    ax.add_patch(mpatches.Rectangle(xy=[r[0], r[2]], width=abs(r[0]-r[1]), height=abs(r[3]-r[2]),
                                    facecolor=None, edgecolor='lime', fill=False, linewidth=2,
                                    transform=ccrs.PlateCarree()))

    # Add lat/lon grid
    gl = ax.gridlines(draw_labels=True, dms=True,  linewidth=0.5, color='k', alpha=0.8)
    gl.xlabel_style = {'fontsize': font_tick}
    gl.ylabel_style = {'fontsize': font_tick}
    gl.top_labels = False
    gl.right_labels = False
    ax.set_title('{}:\nClimatological {} Bias Rank ({})'.format(focus_model,var_name,da.time_avg), fontsize=font_title)

def Scatter_panel(da1, da2, focus_model, ax):
    # Fontsizes
    font_title = 14
    font_label = 12
    font_tick = 12
    font_stat = 12

    # Strip out focus model
    da1_foc = da1.sel(model=focus_model)
    da2_foc = da2.sel(model=focus_model)
    da1 = da1.drop_sel(model=focus_model)
    da2 = da2.drop_sel(model=focus_model)

    mymarkers=[".","o","v","^","<",">","s","p","P","h","H","+","x","X","D","d","1","2"]
    mycolors = ['r','b','g','m']

    for i, key in enumerate(da1.model):  # Add each model to the scatter plot
        thismodel = da1.sel(model=key).model.item()  # Label for legend
        thisplot = ax.scatter(da1.sel(model=key), da2.sel(model=key), marker=mymarkers[i], color=mycolors[i % 4],
                                  s=150, linewidths=1, label=thismodel)
    ax.scatter(da1_foc,da2_foc,marker='*',color='yellow',edgecolor='k',s=250,label=da1_foc.model.item())

    # Plot the line of best fit
    fitcorr = stats.pearsonr(da1.values, da2.values)
    yfit = np.poly1d(np.polyfit(da1.values, da2.values, 1))
    ax.plot(da1.values, yfit(da1.values), color='k', label='r={:.2f} (p={:.1e})'.format(fitcorr[0].item(),fitcorr[1].item()))
    handles, labels = ax.get_legend_handles_labels()
    
    #Labels
    ax.set_title('Climatological Bias compared to OMIP{}'.format(da1.OMIP.values), fontsize=font_title)
    ax.set_xlabel('{} ({}), [{}]'.format(da1.name,da1.units,da1.region), fontsize=font_label)
    ax.set_ylabel('{} ({}), [{}]'.format(da2.name,da2.units,da2.region), fontsize=font_label)
    ax.tick_params(axis='both', labelsize=font_tick)
    ax.grid(True)
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

    return handles,labels

def region_str(region):
    if (region[0]>180):
        x0 = str(360-region[0])+r'$^{\circ}$W'
    else:
        x0 = str(region[0])+r'$^{\circ}$E'
    if (region[1]>180):
        x1 = str(360-region[1])+r'$^{\circ}$W'
    else:
        x1 = str(region[1])+r'$^{\circ}$E'
    if (region[2]>0):
        y0 = str(region[2])+r'$^{\circ}$N'
    else:
        y0 = str(abs(region[2]))+r'$^{\circ}$S'
    if (region[3]>0):
        y1 = str(region[3])+r'$^{\circ}$N'
    else:
        y1 = str(abs(region[3]))+r'$^{\circ}$S'
    regstr = '{}-{}, {}-{}'.format(x0,x1,y0,y1)
    return regstr

def error_stats(da,region,fregion):
    varname = da.name
    ds_out = da.rename(varname+'_bias').to_dataset()
    regweights = np.cos(np.deg2rad(da.lat)) #.drop_attrs()
    da_freg = da.sel(lat=slice(fregion[2], fregion[3]), lon=slice(fregion[0], fregion[1]))
    fregweights = np.cos(np.deg2rad(da_freg.lat)) #.drop_attrs()
    ds_out[varname+'_bias'] = da_freg.weighted(fregweights).mean(['lat', 'lon'], skipna=True, keep_attrs=True)
    ds_out[varname+'_bias'] = ds_out[varname+'_bias'].assign_attrs({'description':'Focus region mean bias',
                                                                   'region':region_str(fregion)})
    ds_out[varname+'_rmse'] = ((da**2).weighted(regweights).mean(['lat', 'lon'], skipna=True))**0.5
    ds_out[varname+'_rmse'] = ds_out[varname+'_rmse'].assign_attrs({'description':'Plot region RMSE','units':da.units,'time_avg':da.time_avg,
                                                                   'region':region_str(region)})
    return ds_out

def plot_preproc(da,region,month):
    monstr = {1:'JAN',2:'FEB',3:'MAR',4:'APR',5:'MAY',6:'JUN',7:'JUL',8:'AUG',9:'SEP',10:'OCT',11:'NOV',12:'DEC'}
    da = da.sel(lat=slice(region[2], region[3]), lon=slice(region[0], region[1]))
    
    # Add 'month' coordinate from 'time'
    if 'time' in da.coords:
        da = da.assign_coords(month=da['time'].dt.month)
        if (month==13):
            # annual mean
            da = da.groupby('month').mean('time').mean('month').assign_attrs({'time_avg': 'Annual'})
            # da = da.mean('month').assign_attrs({'time_avg':'Annual'})
        else:
            # monthly mean
            da = da.groupby('month').mean('time').sel(month=month).assign_attrs({'time_avg': monstr[month]})
            # da = da.sel(month=month).assign_attrs({'time_avg':monstr[month]})
    elif 'month' in da.dims:
        # Already has monthly mean as a dimension
        if month == 13:
            da = da.mean('month').assign_attrs({'time_avg': 'Annual'})
        else:
            da = da.sel(month=month).assign_attrs({'time_avg': monstr[month]})

    else:
        raise ValueError("DataArray must have either 'time' coordinate or 'month' dimension for time averaging")

    return da

def get_units(var):
    unit_dict = {'thetao_zavg':r'$^{\circ}$C','thetao':r'$^{\circ}$C','so_zavg':'psu','so':'psu','sigma0':r'kg m$^{-3}$','sigma0_zavg':r'kg m$^{-3}$','mld':'m'}
    return unit_dict[var]

def SpatialPlot_climo_bias(ds_target, ds_model, ds_obs, var, region=[360-90, 360-0, 20, 80], focus_region=[360-48, 360-30, 38, 53], month=13, save=False, savedir='./'):
    '''
    SpatialPlot_climo_bias generates spatial plots of climatological bias and bias rank relative to a CMIP6 OMIP2 simulation library.

    Parameters:
        ds_target (xarray.Dataset) : Dataset containing the target model data to be analyzed. Must have dimensions: 'model', 'lat', 'lon', 'month'
        ds_model (xarray.Dataset) : Dataset containing the OMIP2 model data library. Must have dimensions: 'model', 'lat', 'lon', 'month'
        ds_obs (xarray.Dataset) : Dataset containing the benchmark observational data. Must have dimensions: 'lat', 'lon', 'month'
        var (str) : The variable to be analyzed. (Select from: ['thetao_zavg', 'so_zavg', 'sigma0_zavg', 'mld'])
        region (list) : The spatial domain for the plot in format [x0,x1,y0,y1] where x0,x1 span [0..360] and y0,y1 span [-90..90]
        focus_region (list) : The spatial subdomain computing regional bias, in format [x0,x1,y0,y1] where x0,x1 span [0..360] and y0,y1 span [-90..90]
        month (int) : The month of climatology to analyze. Set to 13 for annual-average.
        save (bool) : Save figures if True.

    Returns:
        ds_stats (xarray.Dataset) : Dataset containing area-weighted RMSE (over full plot domain) and regional bias (over focus region domain).
    '''
    da_target = plot_preproc(ds_target[var],region,month)
    da_model = plot_preproc(ds_model[var],region,month)
    da_obs = plot_preproc(ds_obs[var],region,month)
    
    # Add target to models
    focus_model = ds_target.model.item()
    da_model = xr.concat([da_model, da_target], dim='model')

    # Compute Bias
    da_model = (da_model - da_obs).assign_attrs({'units':get_units(var),'time_avg':da_obs.time_avg})
    ds_stats = error_stats(da_model,region,focus_region)

    # Plot
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(20, 10), layout='constrained')
    gs = GridSpec(1, 2, figure=fig)

    axs11 = fig.add_subplot(gs[0, 0], projection=proj)
    SpatialBias_panel(da_model, ds_stats, focus_region=focus_region, focus_model=focus_model, ax=axs11)
        
    axs12 = fig.add_subplot(gs[0, 1], projection=proj)
    SpatialRank_panel(da_model, focus_region=focus_region, focus_model=focus_model, ax=axs12)

    # Save Plots
    if save:
        plotname = savedir+'/SpatialPlot_climo_bias.{}.{}.{}.png'.format(focus_model.replace(" ", "_"),var,da_model.time_avg)
        fig.savefig(plotname)

    return ds_stats

def ScatterPlot_Error(ds_x, var_x, ds_y, var_y, focus_model, save=False, savedir='./'):
    '''
    ScatterPlot_Error generates a scatter plot relating error statistics generated through calls to SpatialPlot_climo_bias().

    Parameters:
        ds_x (xarray.Dataset) : Dataset containing the error statistics to be analyzed. Must have dimensions: 'model'
        ds_y (xarray.Dataset) : Dataset containing the error statistics to be analyzed. Must have dimensions: 'model'
        focus_model (str) : Name of focus model 
    '''

    fig = plt.figure(figsize=(10, 5), layout='constrained')
    gs = GridSpec(1, 2, figure=fig)
    
    # Add Spatial Bias Plots
    axs11 = fig.add_subplot(gs[0, 0], label='var1 bias')
    han,lab = Scatter_panel(ds_x[var_x], ds_y[var_y], focus_model, axs11)
    axs12 = fig.add_subplot(gs[0, 1], label='legend')
    axs12.legend(han, lab, ncols=1, loc='center', fontsize=12, markerscale=0.75)  # Add legend to the dedicated axes
    axs12.axis('off')

    # Save Plots
    if save:
        plotname = savedir+'/ScatterPlot_Error.{}.{}_{}.{}_{}.png'.format(focus_model.replace(" ", "_"),var_x,ds_x[var_x].time_avg,var_y,ds_y[var_y].time_avg)
        fig.savefig(plotname)

    return 
