# Michael Wasserstein
# Plot_Radar_Validation_d03_d04.py
# 10/9/2024
# Script takes in radar data from KMTX as well as WRF data from your WRF simulations for d03 and/or d04.
# It then plots them side by side
# Scrpt gives users the option to do a time offset for the model run.

####### Usage #########
# Conda environment - Radar_env_2
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_Radar_Validation_d03_d04.py -r 2 -p 2
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

args = parser.parse_args()

# # Get user inputs
run = str(args.run)
path = int(args.path)
print('Plotting data for run', run)

#different way to format run number
run_number = '{}'.format(run).zfill(2)

# Plot domain 4?
plot_d04 = False

# How far off from the wrf time do you want the radar time to be?
time_offset = False
time_offset_minutes = -45 # minus 45 indicates Radar will be 45 minutes before WRF
time_offset_hours = 0

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
    Fig_dir = parent_dir + '/Figures_{}/wrf_{}/Radar_Validation_d03_d04_time_offset_{}h_{}min/'.format(path,run_number, time_offset_hours, time_offset_minutes)
else:
    Fig_dir = parent_dir + '/Figures_{}/wrf_{}/Radar_Validation_d03_d04/'.format(path,run_number,)
    if plot_d04 == False:
        Fig_dir = parent_dir + '/Figures_{}/wrf_{}/Radar_Validation_d03/'.format(path,run_number,)

# Make the figure directory if it doesn't exist
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

# Paths for the radar analysis
if path == 1:
    KMTX_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group10/michael/Radar/KMTX_lvl2_download2/2022/12/'
elif path in [2,12]:
    KMTX_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group10/michael/Radar/KMTX_lvl2_download2/2019/03/'
    # Start and end time for the period you want to analyze
    start_time_analysis = datetime.datetime(2019,3,22,19,30)
    end_time_analysis = datetime.datetime(2019,3,23,0,15) 
else:
    print('There is in error in your path!')

# Open up all KMTZ files for the date of interest    
KMTX_files = glob.glob(KMTX_path + '*')

# Convert the list of file paths to a DataFrame for easier manipulation
df = pd.DataFrame(KMTX_files, columns=['file_path'])


########### End of user inputs ############

# load in all the wrf output data files
domain = 3
data_files_d03 = glob.glob(WRF_path + '*wrfout_d0{}*'.format(domain)) # for the outermost domain
data_files_d03.sort()

# Find the start and end indicies for the data files of the period of interest for your study
start_ind_d03 = find_time_idx(time=start_time_analysis, files=data_files_d03, domain=3)
end_ind_d03 = find_time_idx(time=end_time_analysis, files=data_files_d03, domain=3)

if plot_d04 == True: # do the same things for domain 4
    domain = 4
    data_files_d04 = glob.glob(WRF_path + '*wrfout_d0{}*'.format(domain)) # for the outermost domain
    data_files_d04.sort()

    start_ind_d04 = find_time_idx(time=start_time_analysis, files=data_files_d04, domain=4)
    end_ind_d04 = find_time_idx(time=end_time_analysis, files=data_files_d04, domain=4)

# Extract only the data files for the WRF run that you are interested in
data_files_d03 = data_files_d03[start_ind_d03:end_ind_d03+1]

if plot_d04:
    data_files_d04 = data_files_d04[start_ind_d04:end_ind_d04+1]


# Loop through all d03 files
for ind, file in enumerate(data_files_d03):
    # open d03 and d04 files
    wrf_file_d03 = Dataset(file)
    if plot_d04:
        wrf_file_d04 = Dataset(data_files_d04[ind]) 
    
    # Stuff for boundaries
    geobounds = wrf.geo_bounds(wrfin=wrf_file_d03)
    bottom_latitude_d03 = geobounds.bottom_left.lat
    left_longitude_d03 = geobounds.bottom_left.lon
    top_latitude_d03 = geobounds.top_right.lat
    right_longitude_d03 = geobounds.top_right.lon
    
    if plot_d04:
        geobounds = wrf.geo_bounds(wrfin=wrf_file_d04)
        bottom_latitude_d04 = geobounds.bottom_left.lat
        left_longitude_d04 = geobounds.bottom_left.lon
        top_latitude_d04 = geobounds.top_right.lat
        right_longitude_d04 = geobounds.top_right.lon

    # Load necessary wrf data
    max_dbz_d03 = getvar(wrf_file_d03, "mdbz", timeidx=-1)
    if plot_d04:
        max_dbz_d04 = getvar(wrf_file_d04, "mdbz", timeidx=-1)
    terr = wrf.getvar(wrf_file_d03, 'ter')
    lake = wrf.getvar(wrf_file_d03, 'LANDMASK')
    
    # Information about the valid wrf time
    valid_time = pd.to_datetime(max_dbz_d03.Time.values)
    valid_time_for_KMTX = valid_time.strftime('KMTX%Y%m%d_%H%-M')
    valid_time_save = datetime.datetime.strftime(valid_time, '%Y%m%d%H%M')
    valid_time_str = datetime.datetime.strftime(valid_time, '%Y-%m-%d %H:%M:%SZ')

    # Get the model initialization time information
    init_time = wrf_file_d03.SIMULATION_START_DATE
    init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
    init_time_str = datetime.datetime.strftime(init_time,  '%Y-%m-%d %H:%M:%SZ')

    # Extract timestamp for KMTX files from file path and add that to the dataaray
    df['timestamp'] = pd.to_datetime(df['file_path'].str.extract(r'(\d{8}_\d{6})')[0], format='%Y%m%d_%H%M%S')

    # Calculate time difference between the time stamp and the wrf file time
    if time_offset:
        valid_time_offset = valid_time + datetime.timedelta(hours = time_offset_hours, minutes = time_offset_minutes) # do time offset
        df['time_difference'] = abs(df['timestamp'] - valid_time_offset)
    else:
        df['time_difference'] = abs(df['timestamp'] - valid_time)

    # Find the file path with the minimum time difference
    f = df.loc[df['time_difference'].idxmin(), 'file_path']
    print(valid_time, f)

    #### DO the NECESSARY PYART STUFF ########
    radar = pyart.io.read_nexrad_archive(f)  # Read in file using pyart

    # Extract date and time for the radar, and create a string of them
    date = radar.time['units'][14:].split('T')[0]
    time = radar.time['units'][14:].split('T')[1]

    # String of KMTX time
    valid_time_KMTX = date + ' ' + time

    radar_location = radar.metadata['instrument_name']  # Get the location of the radar - should be KMTX

    type_of_data = radar.metadata['original_container'] # Get type of data - should be NEXRAD Level II

    radar_lat, radar_lon = radar.latitude['data'][0], radar.longitude['data'][0] # Get lat and lon of the radar

    rad_altitude = radar.altitude['data'][0] # Get the altitude of the radar

    composite = pyart.retrieve.composite_reflectivity(radar) # Pyart function to determine the composite reflectivity - creates a new radar object

    # Extract lat, lon, and composite reflectivity data
    lat_composite = composite.get_gate_lat_lon_alt(sweep = 0)[0]  # only one sweep for this one
    lon_composite = composite.get_gate_lat_lon_alt(sweep = 0)[1]  # only one sweep for this one
    composite_ref = composite.get_field(field_name = 'composite_reflectivity', sweep = 0)  # Only one sweep

    # Get colormap
    cmap = plt.get_cmap('pyart_LangRainbow12')  # Radar colormap - colorblind friendly
    cmap.set_under('none')

    topo_levels = np.arange(1200,3000,500) # specify levels for plotting topograpy as a contour

    datacrs = ccrs.PlateCarree()

    ########################################## Plotting ##########################################
    # Create figure and axis
    fig, (ax1, ax2) = plt.subplots(1,2,figsize = (14, 8), subplot_kw= {'projection' : ccrs.PlateCarree()}, facecolor = 'white', edgecolor = 'k')

    ################### Ax1 KMTX #################
    # Plot Composite Reflectivity data
    plot = ax1.pcolormesh(lon_composite, lat_composite, composite_ref, cmap = cmap, vmin = 0, vmax = 32, zorder = 15)

    # Plot Lakes
    ax1.add_geometries(lakes_gdf.geometry, ccrs.PlateCarree(), zorder = 150, facecolor='none',
                edgecolor = 'black', linewidth = 2.5)

    # Plot Topography
    contour = ax1.contour(lons,lats, topo, zorder = 100, levels = topo_levels, cmap = 'binary', linewidths = 1.3, )#colors = 'black',)

    # Titles
    ax1.set_title('Valid: ' + valid_time_KMTX, loc = 'right')
    ax1.set_title(radar_location +' Composite Reflectivity', loc = 'left')

    # Set limits corresponding with WRF domain
    ax1.set_xlim(left_longitude_d03, right_longitude_d03)
    ax1.set_ylim(bottom_latitude_d03, top_latitude_d03)
    
    if plot_d04: 
        # d04 box
        ax1.add_patch(matplotlib.patches.Rectangle((left_longitude_d04, bottom_latitude_d04), right_longitude_d04-left_longitude_d04, top_latitude_d04-bottom_latitude_d04,
                 fill=None, lw=2, edgecolor='grey', zorder=2000,transform = ccrs.PlateCarree(), alpha = 0.5))

    ################### Ax2 WRF #################
    # Plot composite reflectivity
    # domain 3
    cf = ax2.pcolormesh(max_dbz_d03.XLONG, max_dbz_d03.XLAT, max_dbz_d03.values, cmap = cmap, transform=datacrs, vmin = 0, vmax = 32, zorder = 15)
    if plot_d04:
        # Domain 4
        cf_d04 = ax2.pcolormesh(max_dbz_d04.XLONG, max_dbz_d04.XLAT, max_dbz_d04.values, cmap = cmap, transform=datacrs, vmin = 0, vmax = 32, zorder = 15)

    # Plot the terrain in 500 m increments
    terrain = ax2.contour(terr.XLONG, terr.XLAT, terr,transform=ccrs.PlateCarree(),zorder = 100, levels = topo_levels, cmap = 'binary', linewidths = 1.3, )

    # Add lake
    lake = ax2.contour(lake.XLONG, lake.XLAT, lake,transform=ccrs.PlateCarree(),zorder = 100, levels = [-0.5,0.5], colors = 'black',)

    # Set titles for WRF
    ax2.set_title('WRF{} Run {}\nInit: {}'.format(path, run, init_time_str), loc = 'left')
    if plot_d04:
        ax2.set_title('d03-d04\nValid: {}'.format(valid_time_str), loc = 'right')
    else:
        ax2.set_title('d03\nValid: {}'.format(valid_time_str), loc = 'right')

    # Set limits
    ax2.set_xlim(left_longitude_d03, right_longitude_d03)
    ax2.set_ylim(bottom_latitude_d03, top_latitude_d03)
    
    # d04 box
    if plot_d04:
        ax2.add_patch(matplotlib.patches.Rectangle((left_longitude_d04, bottom_latitude_d04), right_longitude_d04-left_longitude_d04, top_latitude_d04-bottom_latitude_d04,
                 fill=None, lw=2, edgecolor='grey', zorder=2000,transform = ccrs.PlateCarree(), alpha = 0.5))
    
    
    ################ Colorbar #################
    cax = plt.axes([0.125,0.17, 0.77, 0.04])
    cb = plt.colorbar(cf, cax = cax, pad=0.01, aspect=40, extend='max', shrink=0.5, orientation = 'horizontal')
    cb.ax.tick_params(length=8, width=.25, pad=0.01)
    cb.set_label('Composite Reflectivity (dBZ)', labelpad=8, y=0.5, fontsize = 14)
    cb.ax.set_xticks(np.arange(0, 32+.01, 4).astype(int), labels = np.arange(0, 32+.01, 4).astype(str), fontsize=14)

    # Save figure
    plt.savefig(Fig_dir + 'Radar_validation_{}.png'.format(valid_time_save), dpi = 200, bbox_inches = 'tight')
    #plt.show()
    plt.close()
    
