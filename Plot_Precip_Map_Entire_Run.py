#####################################
# Plot_Accumulated_Precip_Map.py
# Michael Wasserstein 03/28/2024
# Updated 10/18/2024
# Script to precipitaiton map for wrf runs
# Usage:    /uufs/chpc.utah.edu/common/home/u1371671/software/pkg/miniconda3/envs/Radar_env_2/bin/python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_Precip_Map_Entire_run.py -r 12 -p 12
# -r represents the run number you want to plot
# -p represents the path of interest
#####################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os, sys
import wrf
import matplotlib
import matplotlib.patheffects as path_effects
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.img_tiles import GoogleTiles
from scipy.ndimage import gaussian_filter

from metpy.plots import USCOUNTIES

import datetime
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim, ALL_TIMES)
from netCDF4 import Dataset
import pyart

sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')

from Useful_Python_Functions import arrow_points, truncate_colormap
#from map_script_2 import *

from PIL import Image
import glob

def timedelta_to_hours(td):
    total_seconds = td.total_seconds()
    hours = total_seconds / 3600
    return hours

############# Get the User inputs ##########
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")

args = parser.parse_args()

# Get user inputs
run = str(args.run)
path = int(args.path)
domain     = 3
run_number = '{}'.format(run).zfill(2)
plot_d04 = False # Do you want to plot data for domain 4 also?

# Star tand end times
if (path == 6) or (path == 12):
    start_time = datetime.datetime(2019,3,22,12,0)
    end_time   = datetime.datetime(2019,3,23,6,0)
    if run in ['14', '16']:
        start_time = datetime.datetime(2019,3,22,12,0)
        end_time   = datetime.datetime(2019,3,23,2,0)
    if run in ['15',]:
        start_time = datetime.datetime(2019,3,22,0,0)
        end_time   = datetime.datetime(2019,3,23,2,0)
if (path == 1) or (path == 9) or (path == 8):
    start_time = datetime.datetime(2022,12,12,0,0)
    end_time   =  datetime.datetime(2022,12,15,0,0)
    
    if run == '13':
        start_time = datetime.datetime(2022,12,13,6,0)
        end_time   =  datetime.datetime(2022,12,14,6,0)

###################################################################
################## Extract info from user inputs ##################
###################################################################
# Time of the accumulation
accumulation_time_td = end_time - start_time
accumulation_time_hours = timedelta_to_hours(accumulation_time_td)

# Formatting for strings in figures
# Reference: https://www.ametsoc.org/index.cfm/ams/publications/author-information/formatting-and-manuscript-components/mathematical-formulas-units-and-time-and-date/#:~:text=Dates%20must%20be%20written%20in,1500%20UTC%203%20May%202015.
#Dates must be written in “Day Month Year” format. Abbreviations for months are permitted only in figures, tables, and their captions.
#For example, AMS style for May 3, 2015 at 3:00 p.m. GMT is 1500 UTC 3 May 2015.
date_str_format = '%H%M UTC %b %d %Y'

# For plotting and writing times
start_time_str = start_time.strftime(date_str_format)
end_time_str = end_time.strftime(date_str_format)

# For extracting WRF files
start_time_WRF = start_time.strftime('%Y-%m-%d_%H:%M:%S')
end_time_WRF = end_time.strftime('%Y-%m-%d_%H:%M:%S')

# For saving files
start_time_save = start_time.strftime('%Y%m%d%H%M')
end_time_save = end_time.strftime('%Y%m%d%H%M')

print('Plotting {} hour precipitation accumulation for WRF{} run {} domain {}'.format(accumulation_time_hours, path,run,domain))
print(f'For {start_time_str} to {end_time_str}')

###################################################################
################## Extact WRF Information/path ####################
###################################################################
# paths for data
base = '/uufs/chpc.utah.edu/common/home/'
home = base + 'u1371671/'
if path ==1:
    base_path = base + 'steenburgh-group12/michael/wrf/'
else:
    base_path = base + 'steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number)
WPS_path = base_path + 'WPS/'

# paths for saving fig
Fig_dir = home + 'WRF/Figures_{}/wrf_{}/Accumulated_Precipitation/'.format(path,run_number)

# Directories for figures if they dont exist
if os.path.exists(os.path.dirname(Fig_dir[:-1])) == False:
    os.mkdir(os.path.dirname(Fig_dir[:-1]))
if os.path.exists(os.path.dirname(Fig_dir)) == False:
    os.mkdir(Fig_dir)

# Path for saving figure
save_path = Fig_dir + 'WRF_Accum_Precip_d{:02}_{:02}h_{}_{}.png'.format(domain, np.round(accumulation_time_hours, 1), start_time_save, end_time_save)

if plot_d04:
    save_path = Fig_dir + 'WRF_Accum_Precip_d03_d04_{:02}h_{}_{}.png'.format(np.round(accumulation_time_hours, 1), start_time_save, end_time_save)

######################## Domain 4 #######################
domain = 3
# The start and the end file for the period of interest
start_file = WRF_path + 'wrfout_d{:02}_{}'.format(domain, start_time_WRF)
end_file = WRF_path + 'wrfout_d{:02}_{}'.format(domain, end_time_WRF)

# Load in the 2 data files
data_files = [start_file, end_file]
wrflist = [Dataset(file) for file in data_files] # get netcdf dataset for each file in a list

# info about the inititalization time
init_time = wrflist[0].SIMULATION_START_DATE
init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
init_time_str = datetime.datetime.strftime(init_time,  date_str_format)
init_time_save = datetime.datetime.strftime(init_time, '%Y%m%d%H%M')

wind_files = glob.glob(WRF_path + '*d0{}*'.format(domain))[24:] # 24 correstponds with hr 6
wrflist_wind = [Dataset(file) for file in wind_files]

# Non-cumulus accumulated rain (total grid-scale precipitation)
rainnc_d03 = wrf.getvar(wrflist, "RAINNC", timeidx=ALL_TIMES, method="cat")
accumulated_rainnc_d03 = rainnc_d03

# Non-cumulus accumulated rain (should be zero for inner domain, since you have no cumulus parameterization)
rainc_d03 = wrf.getvar(wrflist, "RAINC", timeidx=ALL_TIMES, method="cat")
accumulated_rainc_d03 = rainc_d03

# Find total accumulated precip
accumulated_precip_d03 = (accumulated_rainnc_d03 + accumulated_rainc_d03) # extract the zeroth time index
accumulated_precip_d03 = accumulated_precip_d03[1] - accumulated_precip_d03[0] # accumulated_precip_d03[0] should be zero, but we do this as a check


# Extract latitudes and longitudes for the data
latitudes_d03 = rainnc_d03.XLAT.values
longitudes_d03 = rainnc_d03.XLONG.values

# Extract geobounds
geobounds = wrf.geo_bounds(wrfin=Dataset(data_files[0]))
bottom_latitude = geobounds.bottom_left.lat
left_longitude = geobounds.bottom_left.lon
top_latitude = geobounds.top_right.lat
right_longitude = geobounds.top_right.lon

######################## Domain 4 #######################
if plot_d04:
    domain     = 4
    start_file = WRF_path + 'wrfout_d{:02}_{}'.format(domain, start_time_WRF)
    end_file = WRF_path + 'wrfout_d{:02}_{}'.format(domain, end_time_WRF)

    data_files = [start_file, end_file]
    wrflist = [Dataset(file) for file in data_files] # get netcdf dataset for each file in a list

    wind_files = glob.glob(WRF_path + '*d0{}*'.format(domain))[24:] # 24 correstponds with hr 6
    wrflist_wind = [Dataset(file) for file in wind_files]

    # info about the inititalization time
    init_time = wrflist[0].SIMULATION_START_DATE
    init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
    init_time_str = datetime.datetime.strftime(init_time,  date_str_format)
    init_time_save = datetime.datetime.strftime(init_time, '%Y%m%d%H%M')

    # Non-cumulus accumulated rain (total grid-scale precipitation)
    rainnc_d04 = wrf.getvar(wrflist, "RAINNC", timeidx=ALL_TIMES, method="cat")
    accumulated_rainnc_d04 = rainnc_d04

    # Non-cumulus accumulated rain (should be zero for inner domain, since you have no cumulus parameterization)
    rainc_d04 = wrf.getvar(wrflist, "RAINC", timeidx=ALL_TIMES, method="cat")
    accumulated_rainc_d04 = rainc_d04

    # Find total accumulated precip
    accumulated_precip_d04 = (accumulated_rainnc_d04 + accumulated_rainc_d04) # extract the zeroth time index
    accumulated_precip_d04 = accumulated_precip_d04[1] - accumulated_precip_d04[0]

    # Extract latitudes and longitudes for the data
    latitudes_d04 = rainnc_d04.XLAT.values
    longitudes_d04 = rainnc_d04.XLONG.values

    # Extract geobounds
    geobounds = wrf.geo_bounds(wrfin=Dataset(data_files[0]))
    bottom_latitude_d04 = geobounds.bottom_left.lat
    left_longitude_d04 = geobounds.bottom_left.lon
    top_latitude_d04 = geobounds.top_right.lat
    right_longitude_d04 = geobounds.top_right.lon

    # Mask out the precipitation for d03 for where you'll be showing d04
    mask = (
        (latitudes_d03 > bottom_latitude_d04) & (latitudes_d03 < top_latitude_d04) &
        (longitudes_d03 > left_longitude_d04) & (longitudes_d03 < right_longitude_d04)
    )
    accumulated_precip_d03 = np.where(mask == False, accumulated_precip_d03, np.nan)

# Define variables for plotting
datacrs = ccrs.PlateCarree()
fsize = 12
image_res = 12
adj = 0
lw = 1

# Levels and colormap
levels = [1,2,3,4,5,6,7,9,11,13,15,18,] # Specify levels - this could take some work
cmap = matplotlib.cm.get_cmap('pyart_HomeyerRainbow')
vals = np.linspace(0,1,len(levels))
cols = [cmap(val) for val in vals]

# User specified info. Terrain file must be named region+'terrain.nc'
region = 'nutah'
map_background = 'hillshade'   # 'topo' 'hillshade' or 'hillshade-dark'

# Define information for opentopmap or ESRI hillshade depending on map_background requested
# Terms of use for hillshade and hillshade-dark at https://www.arcgis.com/home/item.html?id=9c5370d0b54f4de1b48a3792d7377ff2
class ShadedReliefESRI(GoogleTiles):
    # shaded relief
    def _image_url(self, tile):
        x, y, z = tile
        if map_background == 'topo':
            url = ('https://a.tile.opentopomap.org/{z}/{x}/{y}.png').format(z=z, y=y, x=x)
        elif map_background == 'hillshade':
            url = ('https://server.arcgisonline.com/ArcGIS/rest/services/' \
               'World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg').format(
               z=z, y=y, x=x)
        elif map_background == 'hillshade-dark':
            url = ('https://services.arcgisonline.com/arcgis/rest/services/' \
                'Elevation/World_Hillshade_Dark/MapServer/tile/{z}/{y}/{x}.jpg').format(
                z=z, y=y, x=x)
        else: 
            sys.exit()
        return url

################################################################################
####################      Plotting Figure       ################################
################################################################################
# Set size 
fig = plt.figure(figsize=(10,8))

# Set projection and extent 
ax1 = plt.axes(projection=datacrs)
ax1.set_extent([left_longitude+adj, right_longitude-adj, bottom_latitude+adj, top_latitude-adj])

# Add Background image
ax1.add_image(ShadedReliefESRI(), image_res)#, zorder=1)

# Add geopolitical boundaries
ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
ax1.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)
ax1.add_feature(USCOUNTIES.with_scale('500k'), linewidth=0.5, edgecolor='black', zorder=100)

# Contour fill
cf = ax1.contourf(longitudes_d03, latitudes_d03, accumulated_precip_d03, transform=datacrs, zorder = 69,levels = levels, colors = cols, #cmap = 'pyart_HomeyerRainbow',
                 alpha = 0.6, extend = 'max')

# Line Contours
ax1.contour(longitudes_d03, latitudes_d03, accumulated_precip_d03,
             transform=datacrs, extend='max', transform_first=True, alpha=1, levels = levels[::], # every other level
             colors = 'black',zorder=50, linewidths = lw)

# If the user wants d04, plot it
if plot_d04:
    # Contour fill
    cf = ax1.contourf(longitudes_d04, latitudes_d04, accumulated_precip_d04, transform=datacrs, zorder = 69,levels = levels, colors = cols, #cmap = 'pyart_HomeyerRainbow',
                     alpha = 0.6, extend = 'max')

    # Line Contours
    ax1.contour(longitudes_d04, latitudes_d04, accumulated_precip_d04,
                 transform=datacrs, extend='max', transform_first=True, alpha=1, levels = levels[::], # every other level
                 colors = 'black',zorder=50, linewidths = lw)
    
    # d04 box
    ax1.add_patch(matplotlib.patches.Rectangle((left_longitude_d04, bottom_latitude_d04), right_longitude_d04-left_longitude_d04, top_latitude_d04-bottom_latitude_d04,
             fill=None, lw=2, edgecolor='black', zorder=2000,transform = ccrs.PlateCarree(), alpha = 0.5))
    

# Plot titles
if plot_d04:
    ax1.set_title( 'Domain 3-4'+ '\n' +"WRF initialized ".format(path, run)+init_time_str ,loc='left', fontsize=fsize)
else:
    ax1.set_title( 'Domain {}'.format(3)+ '\n' +"WRF initialized ".format(path, run)+init_time_str ,loc='left', fontsize=fsize)
ax1.set_title(f"{int(accumulation_time_hours)}-h precipitation\n{start_time_str}\n\u2014{end_time_str}", loc='right', fontsize=fsize)

# Colorbar stuff
cb = plt.colorbar(cf, ax=ax1, shrink=.6, pad=0.02, drawedges=True, orientation='horizontal', ticks = levels)
cb.ax.tick_params(labelsize=fsize)
cb.set_label("Accumulated Precipitation (mm)", fontsize=fsize)

plt.savefig(save_path, dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()