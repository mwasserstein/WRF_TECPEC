# Michael Wasserstein
# Compare_Precip_2_runs.py
# 11/5/2024

# Script to plot precipitation difference in accumulated precip between 2 wrf runs

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Compare_Precip_2_runs.py -r1 16 -r2 14 -p 2
# -r1 represents the run number you subtract from
# -r2 represents the run number you subtract
# -p represents the path of interest

import os, sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from map_script import *
import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors, ListedColormap, LinearSegmentedColormap, BoundaryNorm
import matplotlib.cm as cm
import matplotlib.dates as mdates
import pyart
import metpy.calc as mpcalc
from metpy.units import units

import numpy as np
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair, get_cartopy, ALL_TIMES)
import wrf
import glob
import pandas as pd
import datetime
from datetime import datetime


import urllib.request as req
import os.path
import json
import urllib

# ######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r1", "--run1", help="First WRF run of interest, the one you subtract from")
parser.add_argument("-r2", "--run2", help="Second WRF run of interest, the one you subtract")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")

args = parser.parse_args()

###########################
##### stuff for WRF  ######
###########################
# Get user inputs
path = int(args.path)
run_1 = str(args.run1) # this is the run that you subtract from when plotting
run_2 = str(args.run2) # THis is the subtracted run

# RUn numbers with leading zeros
run_number_1 = f'{run_1}'.zfill(2)
run_number_2 = f'{run_2}'.zfill(2)


# Get the timing for the run (user input)
if path == 1:
    time_start = '202212120000'
    time_end = '202212150000'
elif (path == 2) or (path == 6) or (path == 12):
    time_start = '201903221200'  # verify this is right
    time_end = '201903230600'
    if run_1 in ['14', '16']:
        time_start = '201903221200'  # verify this is right
        time_end = '201903230200'
    if run_1 == '15':
        time_start = '201903220000'  # verify this is right
        time_end = '201903230200'
elif (path == 9) or (path == 8):
    time_start = '202212130000'  # verify this is right
    time_end = '202212140600'

# TImes for the title
start_time_title = datetime.strftime(pd.to_datetime(time_start), '%b %-d, %Y %H:%M UTC')
end_time_title = datetime.strftime(pd.to_datetime(time_start), '%b %-d, %Y %H:%M UTC')

# paths for saving fig
parent_dir = '/uufs/chpc.utah.edu/common/home/u1371671/WRF/'
Fig_dir_1 = parent_dir + 'Figures_{}/wrf_{}/Accumulated_Precipitation/'.format(path,run_number_1)
Fig_dir_2 = parent_dir + 'Figures_{}/wrf_{}/Accumulated_Precipitation/'.format(path,run_number_2)

# Make fig dir if it doesn't exist
if os.path.exists(Fig_dir_1) == False:
    os.mkdir(Fig_dir_1)
    
# Make fig dir if it doesn't exist
if os.path.exists(Fig_dir_2) == False:
    os.mkdir(Fig_dir_2)

##############################################################
################ WRF STUFF for run #1 ########################
##############################################################
domain = 3

# paths for data
if path ==1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
else:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number_1)

# load in all the wrf output data files
data_files = glob.glob(WRF_path + '*wrfout_d0{}*'.format(domain)) # for the innermost domain
data_files.sort()

# Get the init times
init_time = Dataset(data_files[0]).SIMULATION_START_DATE
init_time = datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
init_time_str = datetime.strftime(init_time,  '%b %-d, %Y %H:%M UTC')

# Stuff for boundaries ( should be same for both runs)
geobounds = wrf.geo_bounds(wrfin=Dataset(data_files[0]))
bottom_latitude = geobounds.bottom_left.lat
left_longitude = geobounds.bottom_left.lon
top_latitude = geobounds.top_right.lat
right_longitude = geobounds.top_right.lon

# Load in all the wrf files
wrflist_d03 = [Dataset(file) for file in data_files]

# Extract variables of interest
rain_d03 = wrf.getvar(wrflist_d03, "RAINNC", timeidx=-1) # Get last time cause youre interested in the total
rainc_d03 = wrf.getvar(wrflist_d03, "RAINC", timeidx=-1) # Get last time cause youre interested in the total
terr_d03 = wrf.getvar(wrflist_d03, "ter", timeidx=-1)    # only needed for map
XLAND_d03 = wrf.getvar(wrflist_d03, 'XLAND', timeidx=-1) # only needed for map

# Compute total precipitation
rain_d03_total_1 = rain_d03 + rainc_d03

##############################################################
################ WRF STUFF for run #2 ########################
##############################################################
domain = 3

run_number = '{}'.format(run_2).zfill(2)

# paths for data
if path ==1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
else:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number_2)

# load in all the wrf output data files
data_files = glob.glob(WRF_path + '*wrfout_d0{}*'.format(domain)) # for the innermost domain
data_files.sort()

# Load in all the wrf files
wrflist_d03 = [Dataset(file) for file in data_files]

# Extract variables of interest
rain_d03 = wrf.getvar(wrflist_d03, "RAINNC", timeidx=-1) # Get last time cause youre interested in the total
rainc_d03 = wrf.getvar(wrflist_d03, "RAINC", timeidx=-1) # Get last time cause youre interested in the total

# Compute total precipitation
rain_d03_total_2 = rain_d03 + rainc_d03

# Path for saving figure
save_path_1 = Fig_dir_1 + f'Compare_Precip_{run_1}_{run_2}_d0{domain}_{time_start}_{time_end}.png'
save_path_2 = Fig_dir_2 + f'Compare_Precip_{run_1}_{run_2}_d0{domain}_{time_start}_{time_end}.png'

# Levels for plotting
levels = np.arange(-4,4.01,0.25)

#################### Plotting ########################
# Make figure and axis
fig, ax = plt.subplots(1,1,figsize = (14,12),facecolor = 'white', edgecolor = 'k', subplot_kw = {'projection' : ccrs.PlateCarree()})

# Set projection and extent of plot
ax.set_extent([left_longitude, right_longitude, bottom_latitude, top_latitude])

# # Add terrain to locator map
ax.contour(terr_d03.XLONG.values, terr_d03.XLAT.values, terr_d03.values, colors = 'k', levels = np.arange(0,3500,500),
             transform = ccrs.PlateCarree(), zorder = 2, linewidths = 0.75)

# # Add lake to locator map
ax.contour(XLAND_d03.XLONG.values, XLAND_d03.XLAT.values, XLAND_d03.values, colors = 'black', levels = np.arange(1.5,2.6,1),
             transform = ccrs.PlateCarree(), zorder = 2, linewidths = 2)

# Plot precipitation difference
im = ax.contourf(rain_d03_total_2.XLONG.values, rain_d03_total_2.XLAT.values, rain_d03_total_1.values - rain_d03_total_2.values,
                 cmap = 'pyart_balance', levels = levels,
             transform = ccrs.PlateCarree(), zorder = 1, extend = 'both')

# Add colorbar
cax = plt.axes([0.14,0.25, 0.75, 0.025])
plt.colorbar(im, cax = cax, orientation = 'horizontal', label = f'WRF Run {run_1} \u2212 WRF Run {run_2} Precipitation Difference (mm)', ax = ax, fraction = 0.05,
             ticks = levels[::4], drawedges=True)

# Add a title
ax.set_title("WRF{} Init {}".format(path, init_time_str), fontsize = 10, loc = 'left',)
ax.set_title("WRF Simulated Preciptiation\n{}\u2014{}".format(start_time_title,
                                                            end_time_title), fontsize = 10, loc = 'right',)

# Save figure in path of both runs, show and close
plt.savefig(save_path_1, dpi = 300, bbox_inches = 'tight')
plt.savefig(save_path_2, dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()