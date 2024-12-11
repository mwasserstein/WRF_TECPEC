# Michael Wasserstein
# Plot_Radar_Map.py
# 12/6/2024
# Script takes in wrf data and plots radar reflectivity for domain of interest

####### Usage #########
# Conda environment - Radar_env_2
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_Radar_reflectivity_map.py -r 2 -p 2
# -r represents the run number you want to plot
# -p represents the path of the wrf run

# Imports
import datetime
import tempfile
import pyart

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
import matplotlib.patheffects as PathEffects
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas
from cartopy.io.img_tiles import GoogleTiles
from shapely import geometry
from shapely.geometry import Point, Polygon
import xarray as xr

from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair)
import wrf
import urllib.request
import json
import math
import os, sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')

from map_script_2 import *
import pyproj

# Function to use
def find_time_idx(time, files, domain):
    '''
    Function will find the index in a list of files corresponding to the time you want
    '''
    
    # Loop through the file paths and extract the datetime part
    for idx, file_path in enumerate(files):
        # Extract the datetime part from the file path
        date_str = file_path.split('_d0{}_'.format(domain))[1]
        file_datetime = datetime.datetime.strptime(date_str, '%Y-%m-%d_%H:%M:%S')

        # Compare with the target datetime
        if file_datetime == time:
            break
            
    return idx

######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")
parser.add_argument("-d", "--domain", help="domain to plot")

args = parser.parse_args()

# # Get user inputs
run = str(args.run)
path = int(args.path)
domain = int(args.domain)
print('Plotting data for run', run)

time_offset = False
#different way to format run number
run_number = '{}'.format(run).zfill(2)

# paths for data
if path ==1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
elif path in [2,6,12]:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number)
WPS_path = base_path + 'WPS/'

# paths for saving fig
parent_dir = '/uufs/chpc.utah.edu/common/home/u1371671/WRF'
if time_offset:
    Fig_dir = parent_dir + '/Figures_{}/wrf_{}/Radar_Map_d0{}_time_offset_{}h_{}min/'.format(path,run_number, domain, time_offset_hours, time_offset_minutes)
else:
    Fig_dir = parent_dir + '/Figures_{}/wrf_{}/Radar_Map_d0{}/'.format(path,run_number, domain)
    
# Make the figure directory if it doesn't exist
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

# Paths for the radar analysis
if path == 1:
    KMTX_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group10/michael/Radar/KMTX_lvl2_download2/2022/12/'
elif path in [2,12]:
    KMTX_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group10/michael/Radar/KMTX_lvl2_download2/2019/03/'
    # Start and end time for the period you want to analyze
    start_time_analysis = datetime.datetime(2019,3,22,18,0)
    end_time_analysis = datetime.datetime(2019,3,23,2,0) 
else:
    print('There is in error in your path!')

# Open up all KMTZ files for the date of interest    
KMTX_files = glob.glob(KMTX_path + '*')

# Convert the list of file paths to a DataFrame for easier manipulation
df = pd.DataFrame(KMTX_files, columns=['file_path'])


########### End of user inputs ############

# load in all the wrf output data files
data_files = glob.glob(WRF_path + '*wrfout_d0{}*'.format(domain)) # for the outermost domain
data_files.sort()

# Find the start and end indicies for the data files of the period of interest for your study
start_ind = find_time_idx(time=start_time_analysis, files=data_files, domain=domain)
end_ind = find_time_idx(time=end_time_analysis, files=data_files, domain=domain)


# Extract only the data files for the WRF run that you are interested in
data_files = data_files[start_ind:end_ind+1]

# Get colormap
cmap = plt.get_cmap('pyart_ChaseSpectral')  # Radar colormap - colorblind friendly
cmap.set_under('none')
vmin = -20
vmax = 80
ticks = np.arange(vmin, vmax+.01, 10).astype(int)

# Loop through all d03 files
for ind, file in enumerate(data_files):
    # open file 
    wrf_file = Dataset(file)

    
    # Stuff for boundaries
    geobounds = wrf.geo_bounds(wrfin=wrf_file)
    bottom_latitude = geobounds.bottom_left.lat
    left_longitude = geobounds.bottom_left.lon
    top_latitude = geobounds.top_right.lat
    right_longitude = geobounds.top_right.lon

    # Load necessary wrf data
    max_dbz = getvar(wrf_file, "mdbz", timeidx=-1)
    terr = wrf.getvar(wrf_file, 'ter')
    lake = wrf.getvar(wrf_file, 'LANDMASK')
    
    # Information about the valid wrf time
    valid_time = pd.to_datetime(max_dbz.Time.values)
    valid_time_for_KMTX = valid_time.strftime('KMTX%Y%m%d_%H%-M')
    valid_time_save = datetime.datetime.strftime(valid_time, '%Y%m%d%H%M')
    valid_time_str = datetime.datetime.strftime(valid_time, '%Y-%m-%d %H:%M:%SZ')

    # Get the model initialization time information
    init_time = wrf_file.SIMULATION_START_DATE
    init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
    init_time_str = datetime.datetime.strftime(init_time,  '%Y-%m-%d %H:%M:%SZ')



    topo_levels = np.arange(1200,3000,500) # specify levels for plotting topograpy as a contour

    datacrs = ccrs.PlateCarree()

    ########################################## Plotting ##########################################
    # Create figure and axis
    fig, (ax1) = plt.subplots(1, 1, figsize = (14, 8), subplot_kw= {'projection' : ccrs.PlateCarree()}, facecolor = 'white', edgecolor = 'k')

    
    ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
    ax1.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
    ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)

    ################### Ax2 WRF #################
    # Plot composite reflectivity
    # domain 3
    cf = ax1.pcolormesh(max_dbz.XLONG, max_dbz.XLAT, max_dbz.values, cmap = cmap, transform=datacrs, vmin = vmin, vmax = vmax, zorder = 15)
    # Plot the terrain in 500 m increments
    terrain = ax1.contour(terr.XLONG, terr.XLAT, terr,transform=ccrs.PlateCarree(),zorder = 100, levels = topo_levels, cmap = 'binary', linewidths = 1.3, )

    # Add lake
    lake = ax1.contour(lake.XLONG, lake.XLAT, lake,transform=ccrs.PlateCarree(),zorder = 100, levels = [-0.5,0.5], colors = 'black',)

    # Set titles for WRF
    ax1.set_title('WRF{} Run {}\nInit: {}'.format(path, run, init_time_str), loc = 'left')
    ax1.set_title('d0{} Composite Reflectivity\nValid: {}'.format(domain, valid_time_str), loc = 'right')

    # Set limits
    ax1.set_xlim(left_longitude, right_longitude)
    ax1.set_ylim(bottom_latitude, top_latitude)
    
    ################ Colorbar #################
    cax = plt.axes([0.17,0.04, 0.7, 0.04])
    cb = plt.colorbar(cf, cax = cax, pad=0.01, aspect=40, extend='max', shrink=0.5, orientation = 'horizontal')
    cb.ax.tick_params(length=8, width=.25, pad=0.01)
    cb.set_label('Composite Reflectivity (dBZ)', labelpad=8, y=0.5, fontsize = 14)
    cb.ax.set_xticks(ticks, labels = ticks, fontsize=14)

    # Save figure
    plt.savefig(Fig_dir + 'Radar_composite_d0{}_{}.png'.format(domain, valid_time_save), dpi = 200, bbox_inches = 'tight')
    #plt.show()
    plt.close()
    
