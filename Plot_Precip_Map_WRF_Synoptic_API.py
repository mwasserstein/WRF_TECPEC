# Michael Wasserstein
# Plot_Precip_Map_WRF_Synoptic_API.py
# 11/5/2024

# Script to plot precipitation accumulation from domain 3 as a map
# As well as observations from synoptic API
# https://docs.synopticdata.com/services/

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_Precip_Map_WRF_Synoptic_API.py -r 2 -p 2
# -r represents the run number you want to plot
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

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")

args = parser.parse_args()

###########################
##### stuff for WRF  ######
###########################
# Get user inputs
run = str(args.run)
path = int(args.path)
print('Plotting data for run', run,)

# Leading zeros
run_number = '{}'.format(run).zfill(2)

# Get the timing for the run (user input)
if path == 1:
    time_start = '202212120000'
    time_end = '202212150000'
elif (path == 2) or (path == 6) or (path == 12):
    time_start = '201903221200'  # verify this is right
    time_end = '201903230600'
    if run in ['14', '16', '18']:
        time_start = '201903221200'  # verify this is right
        time_end = '201903230200'
    if run == '15':
        time_start = '201903220000'  # verify this is right
        time_end = '201903230200'
elif (path == 9) or (path == 8):
    time_start = '202212130000'  # verify this is right
    time_end = '202212140600'

plot_d04 = False # do you want to plot data for d04
time_offset = False # do you want time offsetting?
time_offset_minutes = 45 # Positive value indicatees that all the wrf times will become 45 minutes later
time_offset_hours = 0

# TImes for the title
start_time_title = datetime.strftime(pd.to_datetime(time_start), '%b %-d, %Y %H:%M UTC')
end_time_title = datetime.strftime(pd.to_datetime(time_end), '%b %-d, %Y %H:%M UTC')


############################
## stuff for synoptic API ##
############################
API_ROOT = "https://api.synopticdata.com/v2/"
API_TOKEN = "0d0f87d395244930af38c1460e0c1a0f"
within = 30 # how many minutes within the time of interest

# paths for data
if path ==1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
else:
    
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number)

# paths for saving fig
parent_dir = '/uufs/chpc.utah.edu/common/home/u1371671/WRF/'
Fig_dir = parent_dir + 'Figures_{}/wrf_{}/Accumulated_Precipitation/'.format(path,run_number)

# Make fig dir if it doesn't exist
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

##################################################
#################WRF STUFF for domain 3########################
##################################################
domain = 3

# load in all the wrf output data files
data_files = glob.glob(WRF_path + '*wrfout_d0{}*'.format(domain)) # for the innermost domain
data_files.sort()

# Get the init times
init_time = Dataset(data_files[0]).SIMULATION_START_DATE
init_time = datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
init_time_str = datetime.strftime(init_time,  '%Y-%m-%d %H:%M:%SZ')

# Stuff for boundaries
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
rain_d03_total = rain_d03 + rainc_d03

##################################################
#################API STUFF########################
##################################################

# String to use for api request - all stations in requested domain
bbox_string = f'{left_longitude},{bottom_latitude},{right_longitude},{top_latitude}'

# Get the total precipitation accumulation during the period of interest
api_request_url = os.path.join(API_ROOT, "stations/precip")            # https://docs.synopticdata.com/services/precipitation
api_request_url += "?bbox={}&start={}&end={}&pmode=totals&token={}".format(bbox_string, time_start, time_end, API_TOKEN)

# Read the api request
response = req.urlopen(api_request_url)
api_text_data = response.read() # read the api request

use_data = json.loads(api_text_data) # Now you can work with use_data because it is a dictionary of the data the API returned.


# Empty list to store data
stn_precip_accums, stn_lats, stn_lons = [],[],[],

# Loop through all the stations in the bounding box
for i in range(len(use_data['STATION'])): 
    
    # Get station data
    stn_data = use_data['STATION'][i]
    
    # Extract precipitation
    stn_precip_accum  = stn_data['OBSERVATIONS']['precipitation'][0]['total']
    
    # Extract latitude and longitude
    stn_lat, stn_lon = float(stn_data['LATITUDE']), float(stn_data['LONGITUDE'])
    
    # Add the data to your lists
    stn_precip_accums.append(stn_precip_accum)
    stn_lats.append(stn_lat)
    stn_lons.append(stn_lon)
    
# Convert to np.arrays
stn_precip_accums = np.array(stn_precip_accums)
stn_lats = np.array(stn_lats)
stn_lons = np.array(stn_lons)

# Get indicies that would sort the array
arg_sort = np.argsort(stn_precip_accums)

# Sort precipitation, lat, and lon # This will be beneficial for plotting cause highest amounts will plot on top
stn_precip_accums = stn_precip_accums[arg_sort]
stn_lats = stn_lats[arg_sort]
stn_lons = stn_lons[arg_sort]

################################################################################
##################### Scatter Plot #############################################
################################################################################
# Path for saving figure
save_path = Fig_dir + f'WRF_Synoptic_precip_accum_scatt_d0{domain}_{time_start}_{time_end}.png'

# Plot settings
levels = np.arange(0,12.1,0.25) # Specify levels - this could take some work
cmap = matplotlib.cm.get_cmap('pyart_HomeyerRainbow')
vals = np.linspace(0,1,len(levels))
cols = [cmap(val) for val in vals]

#################### Plotting ########################
fig, ax = plt.subplots(1,1,figsize = (14,12),facecolor = 'white', edgecolor = 'k', subplot_kw = {'projection' : ccrs.PlateCarree()})

# Set projection and extent of plot
ax.set_extent([left_longitude, right_longitude, bottom_latitude, top_latitude])

# # Add terrain to locator map
ax.contour(terr_d03.XLONG.values, terr_d03.XLAT.values, terr_d03.values, colors = 'k', levels = np.arange(0,3500,500),
             transform = ccrs.PlateCarree(), zorder = 2, linewidths = 0.75)

# # Add lake to locator map
ax.contour(XLAND_d03.XLONG.values, XLAND_d03.XLAT.values, XLAND_d03.values, colors = 'black', levels = np.arange(1.5,2.6,1),
             transform = ccrs.PlateCarree(), zorder = 2, linewidths = 2)

# Plot total precipitation accumulation
im = ax.contourf(rain_d03_total.XLONG.values, rain_d03_total.XLAT.values, rain_d03_total.values, cmap = cmap, levels = levels,
             transform = ccrs.PlateCarree(), zorder = 1, extend = 'max')

# Get info about plot so that colorbar and levels for scatter is correct
cmap = ListedColormap(im.get_cmap()(im.norm(im.cvalues)))
norm = BoundaryNorm(levels, len(levels) - 1)

# Plot station precipitation amounts
scatter = ax.scatter(stn_lons, stn_lats, c = stn_precip_accums, s = 75, edgecolor = 'black', transform=ccrs.PlateCarree(), zorder = 100,
                 cmap = cmap, norm = norm, linewidth= 2)

# Add colorbar
cax = plt.axes([0.14,0.1, 0.75, 0.025])
plt.colorbar(im, cax = cax, orientation = 'horizontal', label = 'Precipitation accumulation (mm)', ax = ax, fraction = 0.05,
             ticks = levels[::8], drawedges=True)

# Add a title
ax.set_title("WRF{} Run {}\nInit {}".format(path, run, init_time_str), fontsize = 10, loc = 'left',)
ax.set_title("WRF simulated and Synoptic API observed precipitation accumulation\n{}\u2014{}".format(start_time_title,
                                                                                                    end_time_title), fontsize = 10, loc = 'right',)

# Save figure, show and close
plt.savefig(save_path, dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()

################################################################################
##################### Text  Plot ###############################################
################################################################################
# Path for saving figure
save_path = Fig_dir + f'WRF_Synoptic_precip_accum_text_d0{domain}_{time_start}_{time_end}.png'
print(save_path)

# Plot settings
levels = np.arange(0,12.1,0.25) # Specify levels - this could take some work
cmap = matplotlib.cm.get_cmap('pyart_HomeyerRainbow')
vals = np.linspace(0,1,len(levels))
cols = [cmap(val) for val in vals]

#################### Plotting ########################
fig, ax = plt.subplots(1,1,figsize = (14,12),facecolor = 'white', edgecolor = 'k', subplot_kw = {'projection' : ccrs.PlateCarree()})

# Set projection and extent of plot
ax.set_extent([left_longitude, right_longitude, bottom_latitude, top_latitude])

# # Add terrain to locator map
ax.contour(terr_d03.XLONG.values, terr_d03.XLAT.values, terr_d03.values, colors = 'k', levels = np.arange(0,3500,500),
             transform = ccrs.PlateCarree(), zorder = 2, linewidths = 0.75)

# # Add lake to locator map
ax.contour(XLAND_d03.XLONG.values, XLAND_d03.XLAT.values, XLAND_d03.values, colors = 'black', levels = np.arange(1.5,2.6,1),
             transform = ccrs.PlateCarree(), zorder = 2, linewidths = 2)

# Plot total precipitation accumulation
im = ax.contourf(rain_d03_total.XLONG.values, rain_d03_total.XLAT.values, rain_d03_total.values, cmap = cmap, levels = levels,
             transform = ccrs.PlateCarree(), zorder = 1, extend = 'max')

# Get info about plot so that colorbar and levels for scatter is correct
cmap = ListedColormap(im.get_cmap()(im.norm(im.cvalues)))
norm = BoundaryNorm(levels, len(levels) - 1)

# Plot station precipitation amounts as a text
for i in range(len(stn_precip_accums)):
    text = ax.text(stn_lons[i], stn_lats[i], s = round(stn_precip_accums[i], 1), size = 8,
                      transform=ccrs.PlateCarree(), zorder = 100,
                     horizontalalignment='center',
                     verticalalignment='center', 
                      bbox=dict(boxstyle="round",
                   ec='k',
                   fc=cmap(stn_precip_accums[i] / 12), # This code gets the color in the colormap corresponding with the precipitation accumulation
                   ))                                  # It works because the max levels is 12. ###### Use 'gainsboro' if you just want grey

# Add colorbar
cax = plt.axes([0.14,0.1, 0.75, 0.025])
plt.colorbar(im, cax = cax, orientation = 'horizontal', label = 'Precipitation accumulation (mm)', ax = ax, fraction = 0.05,
             ticks = levels[::8], drawedges=True)

# Add a title
ax.set_title("WRF{} Run {}\nInit {}".format(path, run, init_time_str), fontsize = 10, loc = 'left',)
ax.set_title("WRF simulated and Synoptic API observed precipitation accumulation\n{}\u2014{}".format(start_time_title,
                                                                                                    end_time_title), fontsize = 10, loc = 'right',)

# Save figure, show and close
plt.savefig(save_path, dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()