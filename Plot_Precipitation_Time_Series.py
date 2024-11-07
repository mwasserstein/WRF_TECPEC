# Michael Wasserstein
# Plot_Precipiation_Time_Series.py
# 11/5/2024

# Script to plot precipitation accumulation from domains 2-3 for the TECPEC simulation 
# For many different stations as well as the observations from synoptic API
# https://docs.synopticdata.com/services/

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_Precipiation_Time_Series.py -r 2 -p 2 -P T
# -r represents the run number you want to plot
# -p represents the path of interest
# -P is a question of if you want to use the Synoptic API precipitation feature or the time series feature. (they are different)

import os, sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from map_script import *
import matplotlib.dates as mdates

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors, ListedColormap, LinearSegmentedColormap
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
parser.add_argument("-P", "--precip", help="Do you want to use the synoptic API precip feature or not? ('T' or 'F')")

args = parser.parse_args()

###########################
##### stuff for WRF  ######
###########################
# Get user inputs
run = str(args.run)
path = int(args.path)
precip = str(args.precip)
print('Plotting data for run', run,)

# Leading zeros
run_number = '{}'.format(run).zfill(2)

# Stuff for times
if path == 1:
    time_start = '202212120000'
    time_end = '202212150000'
elif (path == 2) or (path == 6) or (path == 12):
    time_start = '201903221200'  # verify this is right
    time_end = '201903230600'
    if run in ['14', '16']:
        time_start = '201903221200'  # verify this is right
        time_end = '201903230200'
    if run == '15':
        time_start = '201903220000'  # verify this is right
        time_end = '201903230200'
elif (path == 9) or (path == 8):
    time_start = '202212130000'  # verify this is right
    time_end = '202212140600'
print('Plotting data for run', run,)

run_number = '{}'.format(run).zfill(2)

plot_d04 = False # do you want to plot data for d04
time_offset = False # do you want time offsetting?
time_offset_minutes = 45 # Positive value indicatees that all the wrf times will become 45 minutes later
time_offset_hours = 0

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
if precip == 'T':
    Fig_dir = parent_dir + 'Figures_{}/wrf_{}/stn_precip_time_series_API_precip/'.format(path,run_number)
else:
    Fig_dir = parent_dir + 'Figures_{}/wrf_{}/stn_precip_time_series_API/'.format(path,run_number)

# Make fig dir if it doesn't exist
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

# Settings for plot
color = '#ed3232'
hour_interval = 2

##################################################
#################WRF STUFF for domain 2########################
##################################################
domain = 2

# load in all the wrf output data files
data_files = glob.glob(WRF_path + '*wrfout_d0{}*'.format(domain)) # for the innermost domain
data_files.sort()

# Use the file names to get the start and end time for your plotting output
start_time = datetime.strptime(data_files[0].split('/')[-1][-19:], '%Y-%m-%d_%H:%M:%S')
end_time = datetime.strptime(data_files[-1].split('/')[-1][-19:], '%Y-%m-%d_%H:%M:%S')

# Load in all the wrf files
wrflist_d02 = [Dataset(file) for file in data_files]

# Get the init times
init_time = wrflist_d02[0].SIMULATION_START_DATE
init_time = datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
init_time_str = datetime.strftime(init_time,  '%Y-%m-%d %H:%M:%SZ')

# Extract variables of interest
rain_d02 = wrf.getvar(wrflist_d02, "RAINNC", timeidx=ALL_TIMES, method="cat")
rainc_d02 = wrf.getvar(wrflist_d02, "RAINC", timeidx=ALL_TIMES, method="cat")
terr_d02 = wrf.getvar(wrflist_d02, "ter", timeidx=-1)

# Get the times for the data
times_d02 = pd.to_datetime(rainc_d02.Time.values)

##################################################
#################WRF STUFF for domain 3########################
##################################################
domain = 3

# load in all the wrf output data files
data_files = glob.glob(WRF_path + '*wrfout_d0{}*'.format(domain)) # for the innermost domain
data_files.sort()

# Stuff for boundaries
geobounds = wrf.geo_bounds(wrfin=Dataset(data_files[0]))
bottom_latitude = geobounds.bottom_left.lat
left_longitude = geobounds.bottom_left.lon
top_latitude = geobounds.top_right.lat
right_longitude = geobounds.top_right.lon

# Load in all the wrf files
wrflist_d03 = [Dataset(file) for file in data_files]

# Extract variables of interest
rain_d03 = wrf.getvar(wrflist_d03, "RAINNC", timeidx=ALL_TIMES, method="cat")
rainc_d03 = wrf.getvar(wrflist_d03, "RAINC", timeidx=ALL_TIMES, method="cat")
terr_d03 = wrf.getvar(wrflist_d03, "ter", timeidx=-1)
XLAND_d03 = wrf.getvar(wrflist_d03, 'XLAND', timeidx=-1) # only needed for locator map

# Get the times for the data
times_d03 = pd.to_datetime(rainc_d03.Time.values)


##################################################
#################API STUFF########################
##################################################

# String to use for api request - all stations in requested domain
bbox_string = f'{left_longitude},{bottom_latitude},{right_longitude},{top_latitude}'

# let's get some latest data

if precip == 'T':
    # Get api url for request
    api_request_url = os.path.join(API_ROOT, "stations/precip")            # https://docs.synopticdata.com/services/precipitation
    api_request_url += "?bbox={}&start={}&end={}&pmode=intervals&interval=hour&token={}".format(bbox_string, time_start, time_end, API_TOKEN)

else:
    # Get api url for request
    api_request_url = os.path.join(API_ROOT, "stations/timeseries")
    api_request_url += "?bbox={}&start={}&end={}&vars=precip_accum&token={}".format(bbox_string, time_start, time_end, API_TOKEN)

# Read the api request
response = req.urlopen(api_request_url)
api_text_data = response.read() # read the api request

use_data = json.loads(api_text_data) # Now you can work with use_data because it is a dictionary of the data the API returned.


for i in range(len(use_data['STATION'])): # Loop through all the stations in the bounding box
    
    # Get station data
    stn_data = use_data['STATION'][i]
    
    if precip == 'T':
        precip_data = np.array([stn_data['OBSERVATIONS']['precipitation'][hr]['total'] for hr in range(len(stn_data['OBSERVATIONS']['precipitation']))])
        stn_precip_accum = np.cumsum(precip_data)
        
        # Get the start times and end times in a list of lists as form [[start_time, end_time], [start_time, end_time], ....]
        times = np.array([[stn_data['OBSERVATIONS']['precipitation'][hr]['first_report'], 
            stn_data['OBSERVATIONS']['precipitation'][hr]['last_report']] for hr in range(len(stn_data['OBSERVATIONS']['precipitation']))])
        
        # Convert time strings to datetime objects
        start_times = [datetime.fromisoformat(t[0].replace("Z", "+00:00")) for t in times]
        end_times = [datetime.fromisoformat(t[1].replace("Z", "+00:00")) for t in times]
        
        stn_time = end_times.copy()
        
    else:
        
        # Get station times and precip if using typical time series
        stn_time = pd.to_datetime(np.array(stn_data['OBSERVATIONS']['date_time'])) # Station time
        stn_precip_accum = np.array(stn_data['OBSERVATIONS']['precip_accum_set_1']) # PRecipitation

    # info about station
    stn_ID, stn_name = stn_data['STID'], stn_data['NAME'] # Extract station name
    stn_lat, stn_lon, stn_elev = float(stn_data['LATITUDE']), float(stn_data['LONGITUDE']), float(stn_data['ELEVATION']) / 3.281 # Extract station lat, lon, and elevation (m)

    # Path for saving figure
    save_path = Fig_dir + f'{stn_ID}_precip_accum.png'

    try:
        # Adjust station accumulated precipitation so that it starts at zero always
        stn_precip_accum -= stn_precip_accum[0]

        ############################### WRF data corresponding with stn #####################################
        # Get x_y for the location of interest
        x_y_d02 = wrf.ll_to_xy(wrflist_d02, latitude = stn_lat, longitude = stn_lon) # some interpolations for the xy of the data

        # Compute total modeled preciptiation
        stn_modeled_precip_d02 = rain_d02[:,x_y_d02[1], x_y_d02[0]].values + rainc_d02[:,x_y_d02[1], x_y_d02[0]].values #+ rain[:,x_y[1], x_y[0]].values + graupel[:,x_y[1], x_y[0]].values + hail[:,x_y[1], x_y[0]].values

        # Get elevation for location of interest
        stn_modeled_elev_d02 = str(np.round(terr_d02[x_y_d02[1], x_y_d02[0]].values, 1))

        # Get x_y for the location of interest
        x_y_d03 = wrf.ll_to_xy(wrflist_d03, latitude = stn_lat, longitude = stn_lon) # some interpolations for the xy of the data

        # Compute total modeled preciptiation
        stn_modeled_precip_d03 = rain_d03[:,x_y_d03[1], x_y_d03[0]].values + rainc_d03[:,x_y_d03[1], x_y_d03[0]].values #+ rain[:,x_y[1], x_y[0]].values + graupel[:,x_y[1], x_y[0]].values + hail[:,x_y[1], x_y[0]].values

        # Get elevation for location of interest
        stn_modeled_elev_d03 = str(np.round(terr_d03[x_y_d03[1], x_y_d03[0]].values, 1))

        #################### Plotting ########################
        fig, ax = plt.subplots(1,1,figsize = (12,6),facecolor = 'white', edgecolor = 'k')

        # Plot model results
        ax.plot(times_d02, stn_modeled_precip_d02, color = color, linewidth = 1, linestyle = '-', 
                label = f'{stn_ID} Simulated d02 - Elev {stn_modeled_elev_d02} m')
        ax.plot(times_d03, stn_modeled_precip_d03, color = color, linewidth = 2, linestyle = '-', 
                label = f'{stn_ID} Simulated d03 - Elev {stn_modeled_elev_d03} m')

        # Plot observations
        ax.plot(stn_time, stn_precip_accum, color = color, linewidth = 2, linestyle = '--', label = f'{stn_ID + " Observed - Elev " + str(round(stn_elev, 1)) + " m"} ')

        # Y label
        ax.set_ylabel('Accumulated LPE (mm)', fontsize = 14)

        # For x ticks
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=hour_interval))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d\n%H:%M'))

        # Set y limit
        ax.set_ylim(0,10)

        # Add titles
        ax.set_title('WRF{} run {}\nInit {} UTC'.format(path,run_number, init_time), loc = 'left')
        ax.set_title(f'{stn_ID + "; " + stn_name} accumulated LPE', loc = 'right')

        # Add a legend and a grid
        plt.legend()
        ax.grid()

        ################# Ax2 Locator Map ######################  
        # Set projection and extent of plot
        ax2 = fig.add_axes( [0.12,0.55,0.3,0.3], projection=ShadedReliefESRI().crs) # Need to add a new axis
        ax2.set_extent([left_longitude, right_longitude, bottom_latitude, top_latitude])

        # Add terrain to locator map
        ax2.contourf(terr_d03.XLONG.values, terr_d03.XLAT.values, terr_d03.values, cmap = 'terrain', levels = np.arange(0,3500,250),
                     transform = ccrs.PlateCarree(), zorder = 1)

        # Add lake to locator map
        ax2.contourf(XLAND_d03.XLONG.values, XLAND_d03.XLAT.values, XLAND_d03.values, colors = 'midnightblue', levels = np.arange(1.5,2.6,1),
                     transform = ccrs.PlateCarree(), zorder = 1)

        # Add locator marker
        ax2.scatter(stn_lon, stn_lat, color = 'red', s = 150, edgecolor = 'black', marker = '*', transform=ccrs.PlateCarree(), zorder = 100)

        # Save figure, show and close
        plt.savefig(save_path, dpi = 300, bbox_inches = 'tight')
        #plt.show()
        plt.close()
    except:
        pass