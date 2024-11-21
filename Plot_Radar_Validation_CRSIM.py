# Michael Wasserstein
# Plot_Radar_Validation_CRSIM.py
# 11/21/2024
# Script takes in radar data from KMTX as well as WRF data post-processed through CR-SIM
# It then plots them side by side for a given scanning angle of interest
# Scrpt gives users the option to do a time offset for the model run.

####### Usage #########
# Conda environment - Radar_env_2
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_Radar_Validation_CRSIM.py -r 2 -p 2
# -r represents the run number you want to plot
# -p represents the path of the wrf run


import numpy as np
import xarray as xr
import pandas as pd
import datetime
import matplotlib.pyplot as plt
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair)
import wrf
import glob
from netCDF4 import Dataset
import cartopy
from cartopy import crs
import pyart
import os, sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')

from map_script_2 import *

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

scan = 0.5 # radar degree scan you want to use

# How far off from the wrf time do you want the radar time to be?
time_offset = True
time_offset_minutes = -45 # minus 45 indicates Radar will be 45 minutes before WRF
time_offset_hours = 0


print('Plotting data for run', run,)
base = '/uufs/chpc.utah.edu/common/home/'
home = base + 'u1371671/'
# paths for data
if path ==1:
    base_path = base + 'steenburgh-group12/michael/wrf/'
else:
    base_path = base + 'steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number)

# paths for saving fig
Fig_dir = home + 'WRF/Figures_{}/wrf_{}/Radar_Validation_CRSIM_d03/'.format(path,run_number)
if time_offset:
    Fig_dir = home + 'WRF/Figures_{}/wrf_{}/Radar_Validation_CRSIM_d03_time_offset_{}h_{}min/'.format(path,run_number, time_offset_hours, time_offset_minutes)

# If it doesn't exist, make figure
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

# CRSIM Data path
CRSIM_out_path = base + 'steenburgh-group12/michael/CRSIM_4/wrf{}/wrf_{}/KMTX/'.format(path, run_number)

# User Inputs
gatefilter = None
field = 'reflectivity'

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

# load in all the wrf output data files
domain = 3
data_files_d03 = glob.glob(WRF_path + '*wrfout_d0{}*'.format(domain)) # for the outermost domain
data_files_d03.sort()

# Find the start and end indicies for the data files of the period of interest for your study
start_ind_d03 = find_time_idx(time=start_time_analysis, files=data_files_d03, domain=3)
end_ind_d03 = find_time_idx(time=end_time_analysis, files=data_files_d03, domain=3)

# Extract only the data files for the WRF run that you are interested in
data_files_d03 = data_files_d03[start_ind_d03:end_ind_d03+1]

# Loop through all d03 files
for ind, file in enumerate(data_files_d03):
    # open d03 and d04 files
    wrf_file_d03 = Dataset(file)
    
    # Stuff for boundaries
    geobounds = wrf.geo_bounds(wrfin=wrf_file_d03)
    bottom_latitude = geobounds.bottom_left.lat
    left_longitude = geobounds.bottom_left.lon
    top_latitude = geobounds.top_right.lat
    right_longitude = geobounds.top_right.lon
    
    # Extract wrf data about the land and water
    terr_d03 = wrf.getvar(wrf_file_d03, 'ter')
    lake = wrf.getvar(wrf_file_d03, 'LANDMASK')
    
    # Information about the valid wrf time
    valid_time = pd.to_datetime(terr_d03.Time.values)
    valid_time_for_KMTX = valid_time.strftime('KMTX%Y%m%d_%H%-M')
    valid_time_save = datetime.datetime.strftime(valid_time, '%Y%m%d%H%M')
    valid_time_str = datetime.datetime.strftime(valid_time, '%Y-%m-%d %H:%M:%SZ')
    
    # Get the times
    init_time = wrf_file_d03.SIMULATION_START_DATE
    init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
    init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')
    
    
    #################### CRSIM Stuff ####################
    # Get the CRSIM file name
    CRSIM_file = CRSIM_out_path + f'crsimout_d03_{datetime.datetime.strftime(valid_time, "%Y-%m-%d_%H:%M:%S")}'

    # Open dataset
    ds_CRSIM = xr.open_mfdataset(CRSIM_file)

    # Extract the relevant data for analysiss and plotting
    elev_deg = ds_CRSIM.elev.values
    xlat = ds_CRSIM.xlat
    xlon = ds_CRSIM.xlong
    Zhh = ds_CRSIM.Zhh.values
    
    # THe way the CRSIM output is structured is not on a set fixed elevation scan - so we need to find the 0.5 degree for each location
    # I do that here:
    # Step 1: Compute the absolute difference with scan in degrees
    abs_diff = np.abs(elev_deg - scan)

    # Step 2: Find the index corresponding to the minimum difference between elevation in degrees and scan along the z-axis
    min_idx = np.argmin(abs_diff, axis=0)

    # Step 3: Use the indices of the elevations for 0.5 degree to gather the reflectivity values
    # Closest to the 0.5 degree scan
    # Take along axis takes the zhh data corresponding to the indices for the 0.5 degree scan
    Zhh_scan = np.take_along_axis(Zhh, min_idx[None, :, :], axis=0).squeeze(0) # Squeeze removes z-axis from data

    ########################## KMTX Stuff ###############
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

    for sweep, angle in enumerate(radar.fixed_angle['data']):        # Loop through sweeps
        print(sweep, angle)
        if np.round(angle, 1) == scan:                                # See if the given sweep is the 0.5 degree scan
            print(f'Sweep for {scan} degree radar scan: ' + str(sweep))
            break                                                    # Stop the loop at the 0.5 deg scan
    print('hello')
    # Get the slice at the sweep
    sweep_slice = radar.get_slice(sweep)

    # grab radar data 
    z = radar.get_field(sweep, field)

    # extract lat lons 
    lon = radar.gate_longitude['data'][sweep_slice, :]
    lat = radar.gate_latitude['data'][sweep_slice, :]

    # get the range and time
    ranges = radar.range["data"]
    time = radar.time["data"]

    # get azimuth
    az = radar.azimuth['data'][sweep_slice]
    # get order of azimuths 
    az_ids = np.argsort(az)

    # reorder azs so they are in order 
    az = az[az_ids]
    z = z[az_ids]
    lon = lon[az_ids]
    lat = lat[az_ids]
    time = time[az_ids]
    
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
    plot = ax1.pcolormesh(lon, lat, z, cmap = cmap, vmin = 0, vmax = 32, zorder = 15)

    # Plot Lakes
    ax1.add_geometries(lakes_gdf.geometry, crs.PlateCarree(), zorder = 150, facecolor='none',
                edgecolor = 'black', linewidth = 2.5)

    # Plot Topography
    contour = ax1.contour(lons,lats, topo, zorder = 100, levels = topo_levels, cmap = 'binary', linewidths = 1.3, )#colors = 'black',)

    # Titles
    ax1.set_title('Valid: ' + valid_time_KMTX, loc = 'right')
    ax1.set_title(radar_location +f' {scan}° Reflectivity', loc = 'left')

    # Set limits corresponding with WRF domain
    ax1.set_xlim(left_longitude, right_longitude)
    ax1.set_ylim(bottom_latitude, top_latitude)

    ################### Ax2 WRF #################
    # Plot composite reflectivity
    # domain 3
    plot = ax2.pcolormesh(xlon, xlat,Zhh_scan, transform = crs.PlateCarree(), cmap = cmap,vmin = 0, vmax = 32,)

    # Plot the terrain in 500 m increments
    terrain = ax2.contour(terr_d03.XLONG, terr_d03.XLAT, terr_d03,transform=crs.PlateCarree(),zorder = 100, levels = topo_levels, cmap = 'binary', linewidths = 1.3, )

    # Add lake
    ax2.contour(lake.XLONG, lake.XLAT, lake,transform=crs.PlateCarree(),zorder = 100, levels = [-0.5,0.5], colors = 'black',)

    # Set titles for WRF
    ax2.set_title('WRF{} Run {}\nInit: {}'.format(path, run, init_time_str), loc = 'left')
    ax2.set_title('d03\nValid: {}'.format(valid_time_str), loc = 'right')

    # Set limits
    ax2.set_xlim(left_longitude, right_longitude)
    ax2.set_ylim(bottom_latitude, top_latitude)


    ################ Colorbar #################
    cax = plt.axes([0.125,0.17, 0.77, 0.04])
    cb = plt.colorbar(plot, cax = cax, pad=0.01, aspect=40, extend='max', shrink=0.5, orientation = 'horizontal')
    cb.ax.tick_params(length=8, width=.25, pad=0.01)
    cb.set_label('Reflectivity (dBZ)', labelpad=8, y=0.5, fontsize = 14)
    cb.ax.set_xticks(np.arange(0, 32+.01, 4).astype(int), labels = np.arange(0, 32+.01, 4).astype(str), fontsize=14)

    # Save figure
    plt.savefig(Fig_dir + 'Radar_validation_{}.png'.format(valid_time_save), dpi = 200, bbox_inches = 'tight')
    #plt.show()
    plt.close()

