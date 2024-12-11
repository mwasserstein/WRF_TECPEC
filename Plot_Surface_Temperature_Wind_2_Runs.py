# Michael Wasserstein
# Plot_Wind_Temperature_2_runs.py
# 11/15/2024
# Script takes in WRF outputs and plots the surface wind and temperature for two runs

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts/Plot_Wind_Temperature_2_runs.py -r 2 -p 2
# -r represents the run number you want to plot

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


# paths for data
if path ==1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
else:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
    
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number1)

# paths for saving fig
Fig_dir1 = '/uufs/chpc.utah.edu/common/home/u1371671/WRF/Figures_{}/wrf_{}/Precip_Comparison_wrf_{}_wrf_{}'.format(path, run_number1, run_number1, run_number2)

# load in all the wrf output data files
data_files_d03_run1 = glob.glob(WRF_path + '*wrfout_d03*') # for the 3rd
data_files_d03_run1.sort()

# Load in all wrf files as datasets
wrflist_d03_run1 = [Dataset(data_files_d03_run1[i]) for i in range(len(data_files_d03_run1))]

# Get the useful data
ter_d03_run1 = getvar(wrflist_d03_run1, "ter", timeidx = 0)
LANDMASK_d03_run1 = getvar(wrflist_d03_run1, "LANDMASK", timeidx = 0)
U_V_10_d03_run1 = getvar(wrflist_d03_run1, "uvmet10", timeidx=ALL_TIMES)
T2_run1 = getvar(wrflist_d03_run1, 'T2', timeidx=ALL_TIMES) - 273.15 # convert from Kelvin to celcius

# Extract u and v
U10_run1 = U_V_10_d03_run1.sel(u_v = 'u')
V10_run1 = U_V_10_d03_run1.sel(u_v = 'v')

# Cartopy projection
cart_proj_run1 = wrf.get_cartopy(var = None, wrfin=wrflist_d03_run1, timeidx = -1)

# Extract geobounds (should be the same for both runs you're comparign)
geobounds = wrf.geo_bounds(wrfin=wrflist_d03_run1[0])
bottom_latitude_d03 = geobounds.bottom_left.lat
left_longitude_d03 = geobounds.bottom_left.lon
top_latitude_d03 = geobounds.top_right.lat
right_longitude_d03 = geobounds.top_right.lon


# paths for data
if path ==1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
else:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
    
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number2)

# paths for saving fig
Fig_dir2 = '/uufs/chpc.utah.edu/common/home/u1371671/WRF/Figures_{}/wrf_{}/Precip_Comparison_wrf_{}_wrf_{}'.format(path, run_number2, run_number1, run_number2)

# load in all the wrf output data files
data_files_d03_run2 = glob.glob(WRF_path + '*wrfout_d03*') # for the 3rd
data_files_d03_run2.sort()

# Load in all wrf files as datasets
wrflist_d03_run2 = [Dataset(data_files_d03_run2[i]) for i in range(len(data_files_d03_run2))]

# Get the useful data
ter_d03_run2 = getvar(wrflist_d03_run2, "ter", timeidx = 0)
LANDMASK_d03_run2 = getvar(wrflist_d03_run2, "LANDMASK", timeidx = 0)
U_V_10_d03_run2 = getvar(wrflist_d03_run2, "uvmet10", timeidx=ALL_TIMES)
T2_run2 = getvar(wrflist_d03_run2, 'T2', timeidx=ALL_TIMES) - 273.15 # convert from Kelvin to celcius

# Extract u and v
U10_run2 = U_V_10_d03_run2.sel(u_v = 'u')
V10_run2 = U_V_10_d03_run2.sel(u_v = 'v')

# Cartopy projection
cart_proj_run2 = wrf.get_cartopy(var = None, wrfin=wrflist_d03_run2, timeidx = -1)

# Get information about times
init_time = wrflist_d03_run2[0].SIMULATION_START_DATE
init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')



# Plot settings
levels = np.arange(-6,15.1,1) # Specify levels - this could take some work
cmap = matplotlib.cm.get_cmap('pyart_HomeyerRainbow')
pad = 0.05

# skip for barbs
sknum = 15
skip=(slice(None,None,sknum),slice(None,None,sknum))

for timeidx in range(0,len(U10_run1.Time.values),4): # Loop through all time indices in the run
    fig, ((ax1, ax2)) = plt.subplots(1,2,facecolor = 'white',edgecolor = 'k', figsize = (16, 8), subplot_kw = {'projection' : cart_proj_run2})


    ################################ ax1 - run1 d03 ############################################
    ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
    ax1.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
    ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)

    # d03 set limits
    ax1.set_extent([left_longitude_d03, right_longitude_d03, bottom_latitude_d03, top_latitude_d03])

    # Plot wind barbs
    ax1.barbs(U10_run1.XLONG.values[skip], U10_run1.XLAT.values[skip], U10_run1[timeidx].values[skip], V10_run1[timeidx].values[skip],
        length=4, zorder=21, color = 'black',transform=ccrs.PlateCarree(),)

    # Plot Temperature
    plot = ax1.contourf(T2_run1.XLONG.values, T2_run1.XLAT.values, T2_run1[timeidx].values, cmap = cmap, 
                         levels = levels, extend = 'both',zorder = 20,transform=ccrs.PlateCarree(),)

    # Add terrain contours
    ax1.contour(ter_d03_run1.XLONG.values, ter_d03_run1.XLAT.values, ter_d03_run1.values, colors = 'k',
               levels = np.arange(1000,4000,500), transform=ccrs.PlateCarree(),zorder = 21, linewdiths = 0.5)

    # Add lake to  map
    ax1.contour(LANDMASK_d03_run1.XLONG.values, LANDMASK_d03_run1.XLAT.values, LANDMASK_d03_run1.values, colors = 'black', levels = np.arange(0.0,0.6,0.5),
                 transform = ccrs.PlateCarree(), zorder = 23, linewidths = 2)
    
    # TIme information
    valid_time = pd.to_datetime(T2_run1.Time.values[timeidx])
    valid_time_title = datetime.datetime.strftime(valid_time, '%b %-d, %Y %H:%M UTC')
    valid_time_save = datetime.datetime.strftime(valid_time, '%Y%m%d%H%M')

    # Title
    ax1.set_title(f'WRF{path} Run {run1} d03', loc = 'left')
    ax1.set_title(f'Valid {valid_time_title}', loc = 'right')

    # Add colorbar
    plt.colorbar(plot, ax = ax1, orientation = 'horizontal', label = '2-m Temperature (°C)', pad = pad)


    # ################################ ax2 - run2 d03 ############################################
    ax2.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
    ax2.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
    ax2.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)

    # d03 set limits
    ax2.set_extent([left_longitude_d03, right_longitude_d03, bottom_latitude_d03, top_latitude_d03])

    # Plot wind barbs
    ax2.barbs(U10_run2.XLONG.values[skip], U10_run2.XLAT.values[skip], U10_run2[timeidx].values[skip], V10_run2[timeidx].values[skip],
        length=4, zorder=21, color = 'black',transform=ccrs.PlateCarree(),)

    # Plot Temperature
    plot = ax2.contourf(T2_run2.XLONG.values, T2_run2.XLAT.values, T2_run2[timeidx].values, cmap = cmap, 
                         levels = levels, extend = 'both',zorder = 20,transform=ccrs.PlateCarree(),)

    # Add terrain contours
    ax2.contour(ter_d03_run2.XLONG.values, ter_d03_run2.XLAT.values, ter_d03_run2.values, colors = 'k',
               levels = np.arange(1000,4000,500), transform=ccrs.PlateCarree(),zorder = 21, linewdiths = 0.5)

    # Add lake to  map
    ax2.contour(LANDMASK_d03_run2.XLONG.values, LANDMASK_d03_run2.XLAT.values, LANDMASK_d03_run2.values, colors = 'black', levels = np.arange(0.0,0.6,0.5),
                 transform = ccrs.PlateCarree(), zorder = 23, linewidths = 2)

    # TIme information
    valid_time = pd.to_datetime(T2_run2.Time.values[timeidx])
    valid_time_title = datetime.datetime.strftime(valid_time, '%b %-d, %Y %H:%M UTC')

    # Title
    ax2.set_title(f'WRF{path} Run {run1} d03', loc = 'left')
    ax2.set_title(f'Valid {valid_time_title}', loc = 'right')

    # Add colorbar
    plt.colorbar(plot, ax = ax2, orientation = 'horizontal', label = '2-m Temperature (°C)', pad = pad)

    # Save, show, and close
    plt.savefig(Fig_dir1 + f'Sfc_Temp_WRF{run2}_WRF{run1}_d03_{valid_time_save}.png', dpi = 300, bbox_inches = 'tight')
    plt.savefig(Fig_dir2 + f'Sfc_Temp_WRF{run2}_WRF{run1}_d03_{valid_time_save}.png', dpi = 300, bbox_inches = 'tight')
    #plt.show()
    plt.close()