# Michael Wasserstein
# Plot_Synoptic_Analysis.py
# 10/10/2023
# Script takes in WRF 4 panel outputs and plots synoptic analysis

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts/Plot_Synoptic_Analysis.py -r 2 -p 2
# -r represents the run number you want to plot
# -p represents the wrf path

import wrf
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim, ALL_TIMES)
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#import proplot as pplt
import glob
import xarray
import metpy.calc as mpcalc
import matplotlib.gridspec as gridspec
from metpy.units import units
import os, sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from map_script import *
import datetime
from scipy.ndimage import gaussian_filter
import pyart

######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")


args = parser.parse_args()

# Get user inputs
run = str(args.run)
path = int(args.path)
print('Plotting data for run', run)

run_number = '{}'.format(run).zfill(2)
# ---------------- END USER SPECIFIED VARIABLES ----------------



# IVT contour fill colors                                                                                                                                                                                                                                 
def ivtcfillcolors():                                                                                                                                                                                                                                     
    ivt_fill_colors = ['#3488c4',                                                                                                                                                                                                                         
                           '#68bca8',                                                                                                                                                                                                                                
                           '#b6ec8c',                                                                                                                                                                                                                                
                           '#ebf15e',
                           '#dfae31',
                           '#d06932',
                           '#be2b32',                                                                                                                                                                                                                                
                           '#950002',
                           '#750546',
                           '#45062b']                                                                                                                                                                                                                                 
    return ivt_fill_colors   

# Precip contour fill colors
def precipcfillcolors():
    precip_fill_colors = ['#3488c4',
                   '#68bca8',
                   '#b6ec8c',
                   '#ebf15e',
                   '#dfae31',
                   '#d06932',
                   '#be2b32',
                   '#950002']
    return precip_fill_colors

g = 9.80665 # m /s 

########### End of user inputs ############
# Base paths
base = '/uufs/chpc.utah.edu/common/home/'
home = base + 'u1371671/'
# paths for data
if path ==1:
    base_path = base + 'steenburgh-group12/michael/wrf/'
else:
    base_path = base + 'steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number)
WPS_path = base_path + 'WPS/'

# paths for saving fig
parent_dir = home + 'WRF'
Fig_dir = parent_dir + '/Figures_{}/wrf_{}/Synoptic_Scale_4_Panel/'.format(path,run_number)
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

# load in all the wrf output data files
data_files_d01 = glob.glob(WRF_path + '*wrfout_d01*') # for the outermost domain
data_files_d01.sort()

# extract the first and last times for the simulation as strings from the file names
start_time_str = data_files_d01[0][-19:]
end_time_str = data_files_d01[-1][-19:]

date_format = '%Y-%m-%d_%H:%M:%S' # format of the string

# Get the start and end times of the simulation as datetimes
start_dt = datetime.datetime.strptime(start_time_str, date_format)
end_dt = datetime.datetime.strptime(end_time_str, date_format)

# User inputs - for what times do you want to plot the results (potentially start and end time of model run)
time_of_interest_list = pd.date_range(start = start_dt, end = end_dt, freq = '1H')

for ind, f in enumerate(data_files_d01[:]):

    # Stuff for the times
    time_of_interest = time_of_interest_list[ind]
    time_of_interest_str = datetime.datetime.strftime(time_of_interest, '%Y-%m-%d_%H:%M:%S')
    time_of_interest_save = datetime.datetime.strftime(time_of_interest, '%Y%m%d%H%M')

    wrfin_d01 = Dataset(f) # Extract the actual file

    wrflist = [Dataset(file) for file in data_files_d01]


    # Initialization times
    init_time = wrfin_d01.SIMULATION_START_DATE
    init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
    init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')
    time_to_plot = time_of_interest
    valid_time_str = datetime.datetime.strftime(time_to_plot, '%b %-d, %Y %H:%M UTC')

    # Extract Necessary variables
    tc = getvar(wrfin_d01, "tc", )
    td = getvar(wrfin_d01, "td", )
    p = getvar(wrfin_d01,"pressure",)
    ua = wrf.getvar(wrfin_d01,"ua", )
    va = wrf.getvar(wrfin_d01,"va")
    wa = wrf.getvar(wrfin_d01,"wa")
    rh = wrf.getvar(wrfin_d01,"rh",)
    rain = wrf.getvar(wrflist, "RAINNC", timeidx=ALL_TIMES, method="cat") + wrf.getvar(wrflist, "RAINC", timeidx=ALL_TIMES, method="cat")
    mslp = wrf.getvar(wrfin_d01, 'slp')
    avo = getvar(wrfin_d01, "avo")
    z = getvar(wrfin_d01, "z", units="dm")
    OLR = getvar(wrfin_d01, 'OLR')
    T2 = getvar(wrfin_d01, 'T2') - 273.15 # convert from K to C

    # Use metpy to calculate specificy humidity
    q = mpcalc.specific_humidity_from_dewpoint(p * units.hPa, td * units.degC)


    # Compute IVT and IVT Vector
    q_times_u = q * ua
    q_times_v = q * va

    q_times_v['bottom_top'] = q_times_v['bottom_top'] * 1000 / g # times 1000 for unit conversion to kg/kg
    q_times_u['bottom_top'] = q_times_u['bottom_top'] * 1000 / g

    qu_int =  q_times_u.integrate('bottom_top')
    qv_int =  q_times_v.integrate('bottom_top')

    IVT = np.sqrt(qu_int ** 2 + qv_int ** 2)


    # Interpolate data for pressure levels
    # 700 mbar
    u_700 = gaussian_filter(interplevel(ua, p, 700), sigma = 2)
    v_700 = gaussian_filter(interplevel(va, p, 700), sigma = 2)
    tc_700 = gaussian_filter(interplevel(tc, p, 700), sigma = 2)
    rh_700 = interplevel(rh, p, 700)
    w_700 = gaussian_filter(interplevel(wa, p, 700), sigma = 6) # Units of m / s

    # 500 mbar
    z_500 = gaussian_filter(interplevel(z, p, 500), sigma = 2) * 10
    avo_500 = gaussian_filter(interplevel(avo, p, 500), sigma = 2)    

    # Gaussian filters for surface stuff
    slp = gaussian_filter(mslp, sigma = 2)
    tc_2m = gaussian_filter(T2, sigma = 4)

    # Stuff for rain ( doing 1 hr rain )
    if ind > 0:
        rain_plot = rain[ind] - rain[ind-1]
    else:
        rain_plot = rain[ind]

    # Get the map projection information
    cart_proj = get_cartopy(rh_700)
    lati, long = latlon_coords(rh_700)


    # some map settings
    zorder_bounds = 50
    plt.rcParams.update({'font.size': 5})
    datacrs = ccrs.PlateCarree()
    shrink = 0.2


    ########################################################################
    ########################## Plot Results ################################
    ########################################################################
    #-----------------------------------------------------
    #---------- Plot 500-mb diagnostic upper left ----------
    #-----------------------------------------------------
    # Set plot axes and figure extent
    fig,ax = plt.subplots(2, 2, figsize=(12, 10))
    ax = plt.subplot(221, projection=cart_proj)

    # Add geopolitical boundaries
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = zorder_bounds)
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = zorder_bounds)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = zorder_bounds)

    # Plot absolute vorticity (multiplying by 10^5)
    clevs_500_avor = list(range(-40, 42, 2))
    colors1 = plt.cm.YlOrRd(np.linspace(0, 1, 20))
    colors2 = plt.cm.BuPu(np.linspace(1, 0, 21))
    colors = np.vstack((colors2, (1, 1, 1, 1), colors1))

    # Plot 500 mbar absolute vorticity
    cf = ax.contourf(long, lati, avo_500, clevs_500_avor, colors=colors, extend='both', 
                     transform=datacrs, transform_first=True)
    cb = plt.colorbar(cf, pad=0.01, aspect=40, extend='both', shrink=shrink,)
    cb.ax.tick_params(length=2, width=.25, pad = 1)
    cb.set_label('Abs Vort (x10$^{-5}$$s^{-1}$)', labelpad=-14, y=1.4)

    #Plot 750-mb upward vertical velocity (red shades) 
    clevs_vv = [1.5,3.,4.5]
    cw = ax.contour(long, lati, w_700*100, clevs_vv, colors=['#fb6a4a','#de2d26','#a50f15'], linewidths=0.5,
                    transform=datacrs, transform_first=True)
    plt.clabel(cw, fmt='%1.1f')

    # Plot 750-mb downward vertical velocity (blue shades)
    clevs_vv = [-4.5,-3,-1.5]
    cw = ax.contour(long, lati, w_700*100, clevs_vv, colors=['#08519c','#3182bd','#6baed6'], linewidths=0.5,
                    transform=datacrs, transform_first=True)
    plt.clabel(cw, fmt='%1.1f')

    # Plot 500-mb geopotential height
    clevs_500_hght = np.arange(0, 8000, 30)
    cs = ax.contour(long, lati, z_500, clevs_500_hght, colors='black', linewidths=0.5,
                    transform=datacrs, transform_first=True)
    plt.clabel(cs, fmt='%d')

    # Plot titles
    plt.title('WRF{} Run {} d01 Init '.format(path, run)+init_time_str+
        '\n500-mb Height (m), 500-mb Abs Vorticity (x10$^{-5}$$s^{-1}$), and 700-mb Vertical Velocity (cm $s^{-1}$)', loc='left', pad=2)

    #-----------------------------------------------------
    #---------- Plot SLP diagnostic upper right ----------
    #-----------------------------------------------------
    
    # Set plot axes with proper projection
    ax = plt.subplot(222, projection=cart_proj)

    # Set background to grey if hr 0 as no OLR or Precip Available
    ax.patch.set_facecolor('#6e6d6d')

    # Add geopolitical boundaries
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = zorder_bounds)
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = zorder_bounds)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = zorder_bounds)

    # Plot olr and precip if not forecast hour 0
    clevs_olr = [*range(160, 380, 4)]
    clevs_precip = np.array([.01,.02,.05,.10,.25,.5,1.,2.]) *4
    colors_precip = precipcfillcolors()

    # Plot OLR
    cf = ax.contourf(long, lati, OLR, cmap = 'gray_r', extend='both', 
                 transform=datacrs, transform_first=True)
    # Plot 3 hr precip
    cf = ax.contourf(long, lati, rain_plot, clevs_precip, colors=colors_precip, extend='max',
                 transform=datacrs, transform_first=True)
    cb = plt.colorbar(cf, orientation='vertical', pad=0.01, aspect=40, extend='max', shrink=shrink)
    cb.ax.tick_params(length=2, width=.25, pad=0.01)
    cb.set_label('1-h Precip (mm)', labelpad=-12, y=1.3)


    # Plot 2-m temperature
    clevs_2_t = np.arange(-48, 48, 4)
    cs = ax.contour(long, lati, tc_2m, clevs_2_t, colors='red', linewidths=0.5,
                    transform=datacrs, transform_first=True)
    plt.clabel(cs, cs.levels[::2], fmt='%d')     # cs.levels[::2] will label every other contour

    # # Plot SLP 
    clevs_slp = np.arange(920, 1060, 4)
    cs = ax.contour(long, lati, slp, clevs_slp, colors='black', linewidths=0.5,
            transform=datacrs, transform_first=True)

    plt.clabel(cs, cs.levels[::2], fmt='%d')     # cs.levels[::2] will label every other contour

    # Plot titles
    plt.title('Valid '+valid_time_str+'\nSea Level Pressure (mb), 2-m Temperature (˚C), OLR, and 1-h Precip (mm)', loc='right', pad=2)

    #-----------------------------------------------------
    # ---------------- 700 mbar --------------------------
    #-----------------------------------------------------
    # Create plot axes with proper projection
    ax = plt.subplot(223, projection=cart_proj)

    # Add geopolitical boundaries
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = zorder_bounds)
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = zorder_bounds)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = zorder_bounds)

    # Plot rh
    clevs_rh = list(range(0, 101, 10))
    colors = plt.cm.BrBG(np.linspace(.2, .8, 11))
    cf = ax.contourf(long, lati, rh_700, clevs_rh, colors=colors, extend='max',
                     transform=datacrs, transform_first=True)
    cb = plt.colorbar(cf, pad=0.01, aspect=40, extend='max', shrink=shrink)
    cb.ax.tick_params(length=2, width=.25, pad=0.01)
    cb.set_label('RH (%)', labelpad=-12, y=1.3)

    # Plot 700-hPa temperature
    clevs_700_t = np.arange(-48, 48, 2)
    cs = ax.contour(long, lati, tc_700, clevs_700_t, colors='red', linewidths=.5,
                    transform=datacrs, transform_first=True)
    plt.clabel(cs, cs.levels[::1], fmt='%d')     # cs.levels[::1] will label every contour

    skip = 1
    # Plot 700-hPa Wind Barbs with regridding to put them on a regular grid (this will also thin the barbs)                                                                                                                                                                     
    plt.barbs(long.values[::skip], lati.values[::skip], u_700[::skip], v_700[::skip], transform=datacrs, 
          length=4, linewidth=0.2, zorder=3, color='black', alpha=1.0,regrid_shape = 20)           

    # Plot titles
    plt.title('700-mbar Temperature (˚C), RH (%), and Wind (Full Barb = 10 kt)', loc='left', pad=2)

    #-----------------------------------------------------
    #---------- Plot IVT diagnostic lower right ----------
    #-----------------------------------------------------
    # Create plot axes with proper projection
    ax = plt.subplot(224, projection=cart_proj)

    # Add geopolitical boundaries
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = zorder_bounds)
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = zorder_bounds)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = zorder_bounds)

    # Plot IVT
    clevs_ivt = [0,50,100,250,500,750,1000,1250,1500,1750]
    colors_ivt = ivtcfillcolors()
    cf = ax.contourf(long, lati, IVT, clevs_ivt, colors=colors_ivt, extend='max',
                     transform=datacrs, transform_first=True)
    cb = plt.colorbar(cf, orientation='vertical',pad=0.01, aspect=40, extend='max', shrink=shrink)
    cb.ax.tick_params(length=2, width=.25, pad=0.01)
    cb.set_label('IVT (kg $m^{-1}$$s^{-1}$)', labelpad=-14, y=1.3)

    # Plot IVT contours (these correspond to Cat 1, 2, 3, 4, 5, and if ever needed 6 and 7)
    clevs_ivt = [250, 500, 750, 1000, 1250, 1500, 1750]
    cs = ax.contour(long, lati, IVT, clevs_ivt, colors='black', linewidths=0.25,
                    transform=datacrs, transform_first=True)
    plt.clabel(cs, cs.levels[::2], fmt='%d') # cs.levels[::2] will label every other contour

    # Plot IVT vectors with regridding to put them on a regular grid (this will also thin the vectors)
    cq = ax.quiver(long.values, lati.values, qu_int.values, qv_int.values, scale=12000, transform=ccrs.PlateCarree(), zorder=2, regrid_shape=25)
    qk = plt.quiverkey(cq, .05, 1.02, 500, '500 kg $m^{-1}$$s^{-1}$', coordinates='axes', labelpos = 'E', labelsep=.01, zorder=5)

    # PLot titles
    plt.title('IVT (kg $m^{-1}$ $s^{-1}$) Magnitude and Vectors', loc='right', pad=2)

    # Tighten up layout
    fig.tight_layout()
    plt.subplots_adjust(hspace = -0.63, wspace = -0.1)

    # Save figure
    #plt.savefig(Fig_dir + '4_panel_{}.png'.format(time_of_interest_save), dpi = 300, bbox_inches = 'tight')
    plt.show()
    plt.close()
        