# Michael Wasserstein
# Plot_Reflectivity_2_Runs.py
# 11/15/2024
# Script takes in WRF outputs and plots the reflectivity and winds at a given level for 2 runs

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts/Plot_Reflectivity_2_Runs.py -r1 22 -r2 19 -p 12
# -r1 respresents the run number to plot (the sensitivity)
# -r2 is the run number for the control
# -p is the path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import glob
import wrf
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim, ALL_TIMES)
from netCDF4 import Dataset
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime
import os, sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from map_script import *
import matplotlib as mpl
import pyart
from matplotlib.colors import ListedColormap


######## User input arguments #############
import argparse
parser = argparse.ArgumentParser(description="Specify WRF runs and their paths.")

parser.add_argument("-r1", "--run1", help="WRF run 1 of interest. The sensitivity")
parser.add_argument("-r2", "--run2", help="WRF run 2 of interest. Typically the control")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")

args = parser.parse_args()

# Get user inputs
run1 = str(args.run1)
run2 = str(args.run2)
path = int(args.path)
print('Plotting data for run', run1, 'and', run2)

# Format runs with leading zeros
run_number1 = '{}'.format(run1).zfill(2)
run_number2 = '{}'.format(run2).zfill(2)

lev = 700 # LEvel at which to do interpolation

# Plot settings
vmin = -20
vmax = 80
# Get colormap
cmap = plt.get_cmap('pyart_ChaseSpectral')  # Radar colormap - colorblind friendly
cmap.set_under('none')
pad = 0.05
barb_c = '#3d3811' # color for wind barbs

# skip for barbs
sknum = 15
skip=(slice(None,None,sknum),slice(None,None,sknum))

if path == 12:
    # List of times to analyze
    time_of_interest_list = pd.date_range(datetime.datetime(2019,3,22,19), datetime.datetime(2019,3,23,0), freq = '15min')

# paths for data
if path ==1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
else:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
    
# Path containing the WRF data
WRF_path1 = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number1)
WRF_path2 = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number2)

# paths for saving fig
Fig_dir1 = '/uufs/chpc.utah.edu/common/home/u1371671/WRF/Figures_{}/wrf_{}/Ref_Comparison_wrf_{}_wrf_{}/'.format(path, run_number1, run_number1, run_number2)
Fig_dir2 = '/uufs/chpc.utah.edu/common/home/u1371671/WRF/Figures_{}/wrf_{}/Ref_Comparison_wrf_{}_wrf_{}/'.format(path, run_number2,run_number1, run_number2)

# Make figure directories
if os.path.exists(Fig_dir1) == False:
    os.mkdir(Fig_dir1)
if os.path.exists(Fig_dir2) == False:
    os.mkdir(Fig_dir2)

# Loop through times of interest
for time_of_interest in time_of_interest_list:
    time_of_interest = pd.to_datetime(time_of_interest)
    time_of_interest_file = datetime.datetime.strftime(time_of_interest, '%Y-%m-%d_%H:%M:%S')
    valid_time_title = datetime.datetime.strftime(time_of_interest, '%b %-d, %Y %H:%M UTC')
    valid_time_save = datetime.datetime.strftime(time_of_interest, '%Y%m%d%H%M')
    
    # Load in information for run 1
    filename1 = WRF_path1 + f'wrfout_d03_{time_of_interest_file}'
    wrfin_run1 = Dataset(filename1)
    
    # Load in information for run 2
    filename2 = WRF_path2 + f'wrfout_d03_{time_of_interest_file}'
    wrfin_run2 = Dataset(filename2)
    
    ##### run 1
    # Get the useful data
    ter_d03_run1 = getvar(wrfin_run1, "ter", timeidx = 0)
    LANDMASK_d03_run1 = getvar(wrfin_run1, "LANDMASK", timeidx = 0)
    U_run1 = getvar(wrfin_run1, "ua", timeidx=0)
    V_run1 = getvar(wrfin_run1, "va", timeidx=0)
    p_run1 = getvar(wrfin_run1, "p", timeidx=0, units = 'hPa')
    dbz_run1 = getvar(wrfin_run1, 'dbz', timeidx=0)

    # Interpolate the data to a set level
    dbz_700_run1 = wrf.interplevel(field3d=dbz_run1, vert=p_run1, desiredlev=lev)
    U_700_run1 = wrf.interplevel(field3d=U_run1, vert=p_run1, desiredlev=lev)
    V_700_run1 = wrf.interplevel(field3d=V_run1, vert=p_run1, desiredlev=lev)
    
    ###### run 2
    # Get the useful data
    ter_d03_run2 = getvar(wrfin_run2, "ter", timeidx = 0)
    LANDMASK_d03_run2 = getvar(wrfin_run2, "LANDMASK", timeidx = 0)
    U_run2 = getvar(wrfin_run2, "ua", timeidx=ALL_TIMES)
    V_run2 = getvar(wrfin_run2, "va", timeidx=ALL_TIMES)
    p_run2 = getvar(wrfin_run2, "p", timeidx=ALL_TIMES, units = 'hPa')
    dbz_run2 = getvar(wrfin_run2, 'dbz', timeidx=ALL_TIMES)

    # Interpolate the data to a set level
    dbz_700_run2 = wrf.interplevel(field3d=dbz_run2, vert=p_run2, desiredlev=lev)
    U_700_run2 = wrf.interplevel(field3d=U_run2, vert=p_run2, desiredlev=lev)
    V_700_run2 = wrf.interplevel(field3d=V_run2, vert=p_run2, desiredlev=lev)

    # Extract geobounds (should be the same for both runs you're comparign)
    geobounds = wrf.geo_bounds(wrfin=wrfin_run2)
    bottom_latitude_d03 = geobounds.bottom_left.lat
    left_longitude_d03 = geobounds.bottom_left.lon
    top_latitude_d03 = geobounds.top_right.lat
    right_longitude_d03 = geobounds.top_right.lon
    
    
    # Cartopy projection
    cart_proj_run2 = wrf.get_cartopy(var = None, wrfin=wrfin_run2, timeidx = -1)

    # Get information about times
    init_time = wrfin_run2.SIMULATION_START_DATE
    init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
    init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')
    
    ####################################################################################
    ####################################### Plotting ###################################
    ####################################################################################
    fig, ((ax1, ax2)) = plt.subplots(1,2,facecolor = 'white',edgecolor = 'k', figsize = (16, 8), subplot_kw = {'projection' : cart_proj_run2})


    ################################ ax1 - run1 d03 ############################################
    ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
    ax1.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
    ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)
    
    # Set facecolor as gainsboro
    ax1.set_facecolor('#9e9d9d')

    # d03 set limits
    ax1.set_extent([left_longitude_d03, right_longitude_d03, bottom_latitude_d03, top_latitude_d03])

    # Plot wind barbs
    ax1.barbs(U_700_run1.XLONG.values[skip], U_700_run1.XLAT.values[skip], U_700_run1.values[skip], V_700_run1.values[skip],
        length=4, zorder=500, facecolor = 'black', edgecolor = '#f7f3f2',transform=ccrs.PlateCarree(),)

    # Plot reflectivity
    cf = ax1.pcolormesh(dbz_700_run1.XLONG, dbz_700_run1.XLAT, dbz_700_run1.values, cmap = cmap, transform=ccrs.PlateCarree(), vmin = vmin, vmax = vmax, zorder = 15)

    # Add terrain contours
    ax1.contour(ter_d03_run1.XLONG.values, ter_d03_run1.XLAT.values, ter_d03_run1.values, colors = 'k',
               levels = np.arange(1000,4000,500), transform=ccrs.PlateCarree(),zorder = 21, linewdiths = 0.5)

    # Add lake to  map
    ax1.contour(LANDMASK_d03_run1.XLONG.values, LANDMASK_d03_run1.XLAT.values, LANDMASK_d03_run1.values, colors = 'black', levels = np.arange(0.0,0.6,0.5),
                 transform = ccrs.PlateCarree(), zorder = 23, linewidths = 2)

    # Title
    ax1.set_title(f'WRF Run {run1} d03\n{lev}-hPa Reflectivity and Wind', loc = 'left')
    ax1.set_title(f'Valid {valid_time_title}', loc = 'right')


    # ################################ ax2 - run2 d03 ############################################
    ax2.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
    ax2.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
    ax2.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)
    
    # Set facecolor as gainsboro
    ax2.set_facecolor('#9e9d9d')

    # d03 set limits
    ax2.set_extent([left_longitude_d03, right_longitude_d03, bottom_latitude_d03, top_latitude_d03])

    # Plot wind barbs
    ax2.barbs(U_700_run2.XLONG.values[skip], U_700_run2.XLAT.values[skip], U_700_run2.values[skip], V_700_run2.values[skip],
        length=4, zorder=500, facecolor = 'black', edgecolor = '#f7f3f2',transform=ccrs.PlateCarree(),)

    # Plot reflectivity
    cf = ax2.pcolormesh(dbz_700_run2.XLONG, dbz_700_run2.XLAT, dbz_700_run2.values, cmap = cmap, transform=ccrs.PlateCarree(), vmin = vmin, vmax = vmax, zorder = 15)

    # Add terrain contours
    ax2.contour(ter_d03_run2.XLONG.values, ter_d03_run2.XLAT.values, ter_d03_run2.values, colors = 'k',
               levels = np.arange(1000,4000,500), transform=ccrs.PlateCarree(),zorder = 21, linewdiths = 0.5)

    # Add lake to  map
    ax2.contour(LANDMASK_d03_run2.XLONG.values, LANDMASK_d03_run2.XLAT.values, LANDMASK_d03_run2.values, colors = 'black', levels = np.arange(0.0,0.6,0.5),
                 transform = ccrs.PlateCarree(), zorder = 23, linewidths = 2)

    # Title
    ax2.set_title(f'WRF Run {run2} d03\n{lev}-hPa Reflectivity and Wind', loc = 'left')
    ax2.set_title(f'Valid {valid_time_title}', loc = 'right')

    
    ################ Colorbar #################
    cax = plt.axes([0.125,0.08, 0.77, 0.04])
    cb = plt.colorbar(cf, cax = cax, pad=0.01, aspect=40, extend='max', shrink=0.5, orientation = 'horizontal')
    cb.ax.tick_params(length=8, width=.25, pad=0.01)
    cb.set_label(f'{lev}-hPa Reflectivity (dBZ)', labelpad=8, y=0.5, fontsize = 14)
    cb.ax.set_xticks(np.arange(-20, 80+.01, 10).astype(int), labels = np.arange(-20, 80+.011, 10).astype(str), fontsize=14)
    
    # Adjust spacing between subplots
    plt.subplots_adjust(wspace = 0.1)

    # Save figure and show
    plt.savefig(Fig_dir1 + f'{lev}_wind_ref_Diff_WRF{run2}_WRF{run1}_d03_{valid_time_save}.png', dpi = 300, bbox_inches = 'tight')
    plt.savefig(Fig_dir2 + f'{lev}_wind_ref_Diff_WRF{run2}_WRF{run1}_d03_{valid_time_save}.png', dpi = 300, bbox_inches = 'tight')
    #plt.show()
    plt.close()