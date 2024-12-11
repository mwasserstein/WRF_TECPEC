# Michael Wasserstein
# Plot_Precipitation_Difference_2_runs_time_step_d03.py
# 11/22/2024
# Script takes in WRF outputs for two different runs and plots hourly precipitation for each as well as the precip difference
# FOr domain 3

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts/Plot_Precipitation_Difference_2_runs_time_step_d03.py -r1 20 -r2 19 -p 12
# -r1 represents the run number you want to plot (typically the sensitivity i.e. no terrain)
# -r2 represents the run number you want to plot (typically the control)
# -p represents the path of the data, both runs should have same path

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
import cartopy.feature as cfeature


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

print('Plotting data for run', run1, 'and', run2)

# paths for data
if path ==1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
else:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
    

#################################### Run 1 Information #########################
# wrf data path
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number1)

# paths for saving fig
Fig_dir1 = '/uufs/chpc.utah.edu/common/home/u1371671/WRF/Figures_{}/wrf_{}/Precipitation_Difference_d03_{}_{}_time_step/'.format(path, run_number1, run_number1, run_number2)

# Make fig directory if it doesnt already exist
if os.path.exists(Fig_dir1) == False:
    os.mkdir(Fig_dir1)

###################### Domain 3 #######################
# load in all the wrf output data files
data_files_d03_run1 = glob.glob(WRF_path + '*wrfout_d03*') # for the 3rd
data_files_d03_run1.sort()
data_files_d03_run1 = data_files_d03_run1[::4] # Select every 4th file because you are using 15 minute time resolution

# Load in all wrf files as datasets
wrflist_d03_run1 = [Dataset(data_files_d03_run1[i]) for i in range(len(data_files_d03_run1))]

# Get the useful data
ter_d03_run1 = getvar(wrflist_d03_run1, "ter", timeidx = 0)
LANDMASK_d03_run1 = getvar(wrflist_d03_run1, "LANDMASK", timeidx = 0)
RAINNC_d03_run1 = getvar(wrflist_d03_run1, "RAINNC", timeidx = ALL_TIMES) # Don't need RAINC cause you use no convective parameterization
cart_proj_run1 = wrf.get_cartopy(var = None, wrfin=wrflist_d03_run1, timeidx = -1)

# Extract geobounds (should be the same for both runs you're comparing so don't need to do this twice)
geobounds = wrf.geo_bounds(wrfin=wrflist_d03_run1[0])
bottom_latitude_d03 = geobounds.bottom_left.lat
left_longitude_d03 = geobounds.bottom_left.lon
top_latitude_d03 = geobounds.top_right.lat
right_longitude_d03 = geobounds.top_right.lon


#################################### Run 2 Information #########################
# wrf data path
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number2)

# paths for saving fig
Fig_dir2 = '/uufs/chpc.utah.edu/common/home/u1371671/WRF/Figures_{}/wrf_{}/Precipitation_Difference_d03_{}_{}_time_step/'.format(path, run_number2, run_number1, run_number2)

# Make fig directory if it doesnt already exist
if os.path.exists(Fig_dir2) == False:
    os.mkdir(Fig_dir2)

###################### Domain 3 #######################
# load in all the wrf output data files
data_files_d03_run2 = glob.glob(WRF_path + '*wrfout_d03*') # for the 3rd
data_files_d03_run2.sort()
data_files_d03_run2 = data_files_d03_run2[::4] # Select every 4th file because you are using 15 minute time resolution

# Load in all wrf files as datasets
wrflist_d03_run2 = [Dataset(data_files_d03_run2[i]) for i in range(len(data_files_d03_run2))]

# Get the useful data
ter_d03_run2 = getvar(wrflist_d03_run2, "ter", timeidx = 0)
LANDMASK_d03_run2 = getvar(wrflist_d03_run2, "LANDMASK", timeidx = 0)
RAINNC_d03_run2 = getvar(wrflist_d03_run2, "RAINNC", timeidx = ALL_TIMES)
cart_proj_run2 = wrf.get_cartopy(var = None, wrfin=wrflist_d03_run2, timeidx = -1)

# Plot settings
levels = np.arange(0.01,5.1,0.25) # Specify levels - this could take some work
cmap = matplotlib.cm.get_cmap('pyart_HomeyerRainbow')
pad = 0.05

# Loop through all time steps (don't include very last cause you cant do an i+1 for that one)
for i in range(RAINNC_d03_run2.shape[0]-1):  
    
    # Times for title
    start_time = datetime.datetime.strftime(pd.to_datetime(RAINNC_d03_run2.Time[i].values), '%b %-d, %Y %H:%M UTC')
    end_time = datetime.datetime.strftime(pd.to_datetime(RAINNC_d03_run2.Time[i+1].values), '%b %-d, %Y %H:%M UTC')
    
    # Times for saving
    start_time_save = datetime.datetime.strftime(pd.to_datetime(RAINNC_d03_run2.Time[i].values), '%Y%m%d%H%M')
    end_time_save = datetime.datetime.strftime(pd.to_datetime(RAINNC_d03_run2.Time[i+1].values), '%Y%m%d%H%M')
    
    # Determine amount of precip in the hour
    RAINC_d03_run2_timestep_precip = RAINNC_d03_run2[i+1].values - RAINNC_d03_run2[i].values
    RAINC_d03_run1_timestep_precip = RAINNC_d03_run1[i+1].values - RAINNC_d03_run1[i].values
    
    # Compute the precipitation difference for the time period of interest
    timestep_precipitation_difference = RAINC_d03_run2_timestep_precip - RAINC_d03_run1_timestep_precip
    
    ################################################ Plotting ###############################
    fig, ((ax1, ax2, ax3)) = plt.subplots(1,3,facecolor = 'white',edgecolor = 'k', figsize = (16, 8), subplot_kw = {'projection' : cart_proj_run2})


    ################################################################################################
    ################################ ax1 - run1 d03-d03 ############################################
    ################################################################################################
    # Add geopolitical boundaries
    ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
    ax1.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
    ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)

    ############################# Domain 2 ##########################
    # d03 set limits
    ax1.set_extent([left_longitude_d03, right_longitude_d03, bottom_latitude_d03, top_latitude_d03])

    # Plot PRecipitation
    plot = ax1.contourf(RAINNC_d03_run1.XLONG.values, RAINNC_d03_run1.XLAT.values, RAINC_d03_run1_timestep_precip, cmap = cmap, 
                         levels = levels, extend = 'max',zorder = 20,transform=ccrs.PlateCarree(), alpha = 1)

    # Add terrain contours
    ax1.contour(ter_d03_run1.XLONG.values, ter_d03_run1.XLAT.values, ter_d03_run1.values, colors = 'k',
               levels = np.arange(1000,4001,500), transform=ccrs.PlateCarree(),zorder = 21, linewdiths = 0.1, alpha = 1)

    # Add lake to  map
    ax1.contour(LANDMASK_d03_run1.XLONG.values, LANDMASK_d03_run1.XLAT.values, LANDMASK_d03_run1.values, colors = 'black', levels = np.arange(0.0,0.6,0.5),
                 transform = ccrs.PlateCarree(), zorder = 23, linewidths = 2, alpha = 1)
    
    # TItle
    ax1.set_title(f'WRF Run {run1} d03')

    # Add colorbar
    plt.colorbar(plot, ax = ax1, orientation = 'horizontal', label = 'Accumulated Precipitation (mm)', pad = pad)

    ##############################################################################################
    # ################################ ax2 - run2 d03-d03 ########################################
    ##############################################################################################
    # Add geopolitical boundaries
    ax2.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
    ax2.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
    ax2.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)

    ############################# Domain 2 ##########################
    # d03 set limits
    ax2.set_extent([left_longitude_d03, right_longitude_d03, bottom_latitude_d03, top_latitude_d03])

    # Plot PRecipitation
    plot = ax2.contourf(RAINNC_d03_run2.XLONG.values, RAINNC_d03_run2.XLAT.values, RAINC_d03_run2_timestep_precip, cmap = cmap, 
                         levels = levels, extend = 'max',zorder = 20,transform=ccrs.PlateCarree(),alpha = 1)

    # Add terrain contours
    ax2.contour(ter_d03_run2.XLONG.values, ter_d03_run2.XLAT.values, ter_d03_run2.values, colors = 'k',
               levels = np.arange(1000,4001,500), transform=ccrs.PlateCarree(),zorder = 21, linewdiths = 0.1,alpha = 1)

    # Add lake to  map
    ax2.contour(LANDMASK_d03_run2.XLONG.values, LANDMASK_d03_run2.XLAT.values, LANDMASK_d03_run2.values, colors = 'black', levels = np.arange(0.0,0.6,0.5),
                 transform = ccrs.PlateCarree(), zorder = 23, linewidths = 2,alpha = 1)

    # TItle
    ax2.set_title(f'WRF Run {run2} d03')

    # Add colorbar
    plt.colorbar(plot, ax = ax2, orientation = 'horizontal', label = 'Accumulated Precipitation (mm)', pad = pad)

    ########################################################################################################
    # ################################ ax3 - difference d03-d03 ############################################
    ########################################################################################################
    # Add geopolitical boundaries
    ax3.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
    ax3.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
    ax3.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)

    ############################# Domain 2 ##########################
    # # d03 set limits
    ax3.set_extent([left_longitude_d03, right_longitude_d03, bottom_latitude_d03, top_latitude_d03])


    # # Plot precipitation difference
    plot = ax3.contourf(RAINNC_d03_run2.XLONG.values, RAINNC_d03_run2.XLAT.values, timestep_precipitation_difference, cmap = 'bwr', 
                        transform=ccrs.PlateCarree(), levels = np.arange(-6,6.2,0.5), extend = 'both')
    
    # Add lake to  map
    ax3.contour(LANDMASK_d03_run2.XLONG.values, LANDMASK_d03_run2.XLAT.values, LANDMASK_d03_run2.values, colors = 'black', levels = np.arange(0.0,0.6,0.5),
                 transform = ccrs.PlateCarree(), zorder = 23, linewidths = 2,alpha = 1)

    # # Title
    ax3.set_title(f'Run {run2} \u2212 Run {run1} precipitation')

    # Colorbar
    cb = plt.colorbar(plot, ax = ax3, orientation = 'horizontal', label = 'Difference', pad = pad)

    # Add a suptitle for the whole plot
    plt.suptitle(f'Hourly Precipitation {start_time}\u2013{end_time}', y = 0.80)

    # Save and show and close
    plt.savefig(Fig_dir1 + f'Precip_Diff_WRF{run2}_WRF{run1}_d03_{start_time_save}_{end_time_save}.png', dpi = 300, bbox_inches = 'tight')
    plt.savefig(Fig_dir2 + f'Precip_Diff_WRF{run2}_WRF{run1}_d03_{start_time_save}_{end_time_save}.png', dpi = 300, bbox_inches = 'tight') 
    #plt.show()
    plt.close()