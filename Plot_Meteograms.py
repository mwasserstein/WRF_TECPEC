# Michael Wasserstein
# Plot_Meteograms.py
# 3/15/2024
# Script takes in WRF outputs and plots a meteogram for the site of interst
# DOmain 3 works best

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_Meteograms.py -r 2 -p 2 -d 3
# -r represents the run number you want to plot
# -p represents the path of interest
# -d represents the domain of interest

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os, sys
import wrf
import matplotlib.patheffects as path_effects

import datetime
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim, ALL_TIMES)
from netCDF4 import Dataset
import matplotlib.dates as mdates
import metpy.calc as mpcalc
from metpy.units import units as units
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors, ListedColormap, LinearSegmentedColormap

import urllib.request as req
import os.path
import json
import urllib
import random

sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from map_script import *

######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")
parser.add_argument("-d", "--domain", help="WRF domain")


args = parser.parse_args()


###########################
##### stuff for WRF  ######
###########################
# Get user inputs
run = str(args.run)
path = int(args.path)
domain = int(args.domain)

time_offset = False
time_offset_hours = 0
time_offset_minutes = 45 # Means WRF time will be 45 minutes later (e.g. real data at 12:00, WRF at 12:45)

# Convert hours to minutes and add them (for title)
total_minutes = (time_offset_hours * 60) + time_offset_minutes

import math

def round_down_to_nearest_5(number):
    return 5 * math.floor(number / 5)

def round_up_to_nearest_5(number):
    return 5 * math.ceil(number / 5)

############################
## stuff for synoptic API ##
############################
API_ROOT = "https://api.synopticdata.com/v2/"
API_TOKEN = "XXXXXXXXXXXXX" # synoptic data api token

# User inputs
if path in [2,6,12]:
    init_time_api = '201903221200'
    end_time_api  = '201903230600'
    hour_interval = 3
    if run in ['14','16', '18']:
        init_time_api = '201903221200'
        end_time_api  = '201903230200' 
    if run in ['15']:
        init_time_api = '201903220000'
        end_time_api  = '201903230200'      
if path in [1,5,8,9]:
    init_time_api = '202212120000'
    end_time_api  = '202212150000'
    hour_interval = 6
    if run == '13':
        init_time_api = '202212130000'
        end_time_api  = '202212140600'    
        hour_interval = 3


run_number = '{}'.format(run).zfill(2)

# Base paths
base = '/uufs/chpc.utah.edu/common/home/'
home = base + 'u1371671/'
# paths for data
if path ==1:
    base_path = base + 'steenburgh-group12/michael/wrf/'
else:
    base_path = base + 'steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number)
WPS_path = base_path + 'WPS/'

# paths for saving fig
if time_offset:
    Fig_dir = home + 'WRF/Figures_{}/wrf_{}/Meteograms_d0{}_time_offset_{}h_{}min/'.format(path,run_number, domain,time_offset_hours, time_offset_minutes)
else:
    Fig_dir = home + 'WRF/Figures_{}/wrf_{}/Meteograms_d0{}/'.format(path,run_number, domain,)    

# Make figure directory if it does not exist
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)


# load in all the wrf output data files
data_files = glob.glob(WRF_path  + '*wrfout_d0{}*'.format(domain)) # for the domain of interest
data_files.sort() # sort files

# Get all wrf files into the list
wrflist = [Dataset(file) for file in data_files]

# Information for start time
init_time = wrflist[0].SIMULATION_START_DATE
init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')
init_time_api = datetime.datetime.strftime(init_time, '%Y%m%d%H%M')

# End time information
end_time = pd.to_datetime(end_time_api)
end_time_str = datetime.datetime.strftime(end_time, '%b %-d, %Y %H:%M UTC')
end_time_api = datetime.datetime.strftime(end_time, '%Y%m%d%H%M')

#########################################################################################
########################  WRF Data ######################################################
#########################################################################################
# Stuff for boundaries
geobounds = wrf.geo_bounds(wrfin=wrflist[0])
bottom_latitude = geobounds.bottom_left.lat
left_longitude = geobounds.bottom_left.lon
top_latitude = geobounds.top_right.lat
right_longitude = geobounds.top_right.lon

# Extract the relevant WRF data 
tc = getvar(wrflist, "T2", timeidx=ALL_TIMES, method="cat")
p = getvar(wrflist,"pressure", timeidx=ALL_TIMES, method="cat")
td = getvar(wrflist, 'td2', timeidx=ALL_TIMES, method="cat")
wind10 = wrf.getvar(wrflist, 'wspd_wdir10', timeidx=ALL_TIMES, method = 'cat')
rh = getvar(wrflist, 'rh2', timeidx=ALL_TIMES, method="cat")
WRF_times = tc.Time

if time_offset:      # make the adjustment for a time offset
    WRF_times = pd.to_datetime(WRF_times.values) + datetime.timedelta(minutes = 45)
#######################################################################################

#########################################################################################
########################  API Data ######################################################
#########################################################################################
times_for_event = pd.date_range(init_time, end_time, freq = '1h', tz = 'utc')

# String to use for api request - all stations in requested domain
bbox_string = f'{left_longitude},{bottom_latitude},{right_longitude},{top_latitude}'

# let's get some latest data   Resource: https://docs.synopticdata.com/services/time-series
api_request_url = os.path.join(API_ROOT, "stations/timeseries") # you want time series data
api_request_url += "?token={}&bbox={}&start={}&end={}".format(API_TOKEN,bbox_string,init_time_api, end_time_api) # if you want data for Utah

response = req.urlopen(api_request_url)
api_text_data = response.read() # read the api request

use_data = json.loads(api_text_data) # Now you can work with use_data because it is a dictionary of the data the API returned.

# List of all the stations within the bounding box
list_stations = [use_data['STATION'][i]['STID'] for i in range(len(use_data['STATION']))]

print(list_stations)

# Loop through all the stations
for station in list_stations[:]:
    save_path = Fig_dir + 'Meteogram_{}.png'.format(station)
    try: 

        # let's get some latest data   Resource: https://docs.synopticdata.com/services/time-series
        api_request_url = os.path.join(API_ROOT, "stations/timeseries") # you want time series data
        # string of api request for the station you want
        api_request_url += "?token={}&stid={}&start={}&end={}".format(API_TOKEN, station, init_time_api, end_time_api)

        response = req.urlopen(api_request_url)
        api_text_data = response.read() # read the api request

        use_data = json.loads(api_text_data) # Now you can work with use_data because it is a dictionary of the data the API returned.

        # Latitude and longitude
        stn_lat = float(use_data['STATION'][0]['LATITUDE'])
        stn_lon = float(use_data['STATION'][0]['LONGITUDE'])

        # check the units of the data you're working with
        print('Temperature units:', use_data['UNITS']['air_temp'])
        print('Wind degree units:', use_data['UNITS']['wind_direction'])
        print('Wind speed units: ', use_data['UNITS']['wind_speed'])
        print('Rel Humid units: ', use_data['UNITS']['relative_humidity'])

        # Extract the times where data are available
        api_times = pd.to_datetime(use_data['STATION'][0]['OBSERVATIONS']['date_time'])

        # get the indicies where the hour is zero
        indicies_for_hour_zero = np.where(api_times.minute == 0)[0]

        # Extract the relevant API data
        temperatures_API = np.array(use_data['STATION'][0]['OBSERVATIONS']['air_temp_set_1'])[indicies_for_hour_zero]
        dewpoint_API = np.array(use_data['STATION'][0]['OBSERVATIONS']['dew_point_temperature_set_1d'])[indicies_for_hour_zero] # set_1d is the "Derived" dewpoint but I thnk that is fine
        winddir_API = np.array(use_data['STATION'][0]['OBSERVATIONS']['wind_direction_set_1'])[indicies_for_hour_zero]
        windspeed_API = np.array(use_data['STATION'][0]['OBSERVATIONS']['wind_speed_set_1'])[indicies_for_hour_zero]
        relative_humidity_API = np.array(use_data['STATION'][0]['OBSERVATIONS']['relative_humidity_set_1'])[indicies_for_hour_zero]
        api_times = api_times[indicies_for_hour_zero]

        # Get the x and y gridpoints for the latitude and longitude of the station
        x_y = wrf.ll_to_xy(wrflist[0], float(stn_lat), float(stn_lon)) # some interpolations for the xy of the data

        # GET the wrf data for the station of interest
        # Prepare data for plotting
        tc_WRF = tc[:,x_y[1], x_y[0]] - 273.15
        td_WRF = td[:,x_y[1], x_y[0]]
        rh_WRF = rh[:,x_y[1], x_y[0]]

        # WInd information
        wspd_WRF = wind10.sel(wspd_wdir = 'wspd')[:,x_y[1], x_y[0]]
        wdir_WRF = wind10.sel(wspd_wdir = 'wdir')[:,x_y[1], x_y[0]]

        #########################################################################################
        ########################  Plotting ######################################################
        #########################################################################################
        fig, (ax1, ax1_, ax2) = plt.subplots(3,1,figsize = (12,10), facecolor = 'white', edgecolor = 'k')

        # Plot temperatures
        ax1.plot(WRF_times, tc_WRF,label="Simulated 2 m Temperature", color="tab:red")
        ax1.plot(api_times, temperatures_API,linestyle = '--', label="Observed 2 m Temperature", color="tab:red")

        # Dewpoint
        ax1.plot(WRF_times,td_WRF, label = 'Simulated 2 m Dew point', color = 'tab:blue')
        ax1.plot(api_times,dewpoint_API, linestyle = '--',label = 'Observed 2 m Dew point', color = 'tab:blue')

        # Y label
        ax1.set_ylabel('Temperature (°C)')
        ax1.set_ylim(round_down_to_nearest_5(td_WRF.min())-5, round_up_to_nearest_5(tc_WRF.max())+5)

        # Y Label and limits
        ax1.xaxis.set_major_locator(mdates.HourLocator(interval=hour_interval))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b %d\n%H:%M'))
        #ax1.set_xlim(api_times[0], pd.to_datetime(WRF_times.values)[-1])

        # Legend
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        # Titles
        if time_offset == True:
            ax1.set_title('{} Meteogram\nOffset {} minutes'.format(station.upper(), total_minutes))
        else:
            ax1.set_title('{} Meteogram'.format(station.upper()))
        ax1.set_title("Init "+init_time_str, loc = 'left')
        ax1.set_title('WRF{} run {} - d0{}'.format(path, run_number, domain), loc = 'right')

        # Grid
        ax1.grid()

        ################# Ax1_ RH data ######################

        # Plot relative humidity
        ax1_.plot(WRF_times, rh_WRF, color = "tab:green", label = 'Simulated 10 m Relative Humidity', zorder = 0)
        ax1_.plot(api_times, relative_humidity_API, linestyle = '--', color = "tab:green", label = 'Observed 10 m Relative Humidity', zorder = 0)

        # For the legend
        ax1_.plot([], [], color = "tab:green", label = 'Simulated 2 m Relative Humidity')
        ax1_.plot([], [], linestyle = '--', color = "tab:green", label = 'Observed 2 m Relative Humidity', zorder = 0)

        # Y ticks and limits and labels
        ax1_.set_yticks(np.arange(round_down_to_nearest_5(rh_WRF.min()),101,5))
        ax1_.set_ylim(round_down_to_nearest_5(rh_WRF.min()),100)
        ax1_.set_ylabel('Relative Humidity (%)')

        # Xlimits and labels
        #ax1_.set_xlim(api_times[0], pd.to_datetime(WRF_times.values)[-1])
        ax1_.xaxis.set_major_locator(mdates.HourLocator(interval=hour_interval))
        ax1_.xaxis.set_major_formatter(mdates.DateFormatter('%b %d\n%H:%M'))

        # Legend
        ax1_.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        # Grid
        ax1_.grid()

        ################# Ax2 wind data ######################  
        # Make a twin axis for the wind directioon
        ax2_ = ax2.twinx()

        # Plot wind speed
        ax2.plot(WRF_times, wspd_WRF, color="tab:red", label = 'Simulated 10 m Wind Speed')
        ax2.plot(api_times, windspeed_API,linestyle = '--', label="Observed 10 m Wind Speed", color="tab:red")

        # WInd direction    
        ax2_.scatter(pd.to_datetime(WRF_times.values), wdir_WRF.values, color="tab:blue",) # not sure why I need to convert to datetime, but it works
        ax2_.scatter(api_times, winddir_API, color="tab:blue",  marker = '*')

        # Stuff for legend
        ax2.scatter([], [], color="tab:blue", label = 'Simulated 10 m Wind Direction') # just for the legend
        ax2.scatter([],[],color="tab:blue", label = 'Observed 10 m Wind Direction', marker = '*')
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        # Y limits, labels, and ticks
        ax2.set_ylim(0,)
        ax2.set_ylabel('Wind Speed (m s$^{-1}$)')

        ax2_.set_ylabel('Wind Direction')
        ax2_.set_yticks(np.arange(0,361,45))
        ax2_.set_yticklabels(['N', 'NE', 'E', 'SE', 'S', "SW", 'W', 'NW', 'N'])

        # X limits, labels, and ticks
        ax2.xaxis.set_major_locator(mdates.HourLocator(interval=hour_interval))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b %d\n%H:%M'))
        #ax2.set_xlim(api_times[0], pd.to_datetime(WRF_times.values)[-1])
        ax2.set_xlabel('Time (UTC)')

        # Grid
        ax2.grid()

        ################# Ax3 Locator Map ######################  
        # Set projection and extent of plot
        ax3 = fig.add_axes( [0.94,0.34,0.2,0.18], projection=ShadedReliefESRI().crs) # Need to add a new axis
        ax3.set_extent([left_longitude, right_longitude, bottom_latitude, top_latitude])

        # Add the map tiles as level 1
        ax3.add_image(ShadedReliefESRI(), map_tile_detail)#, zorder=1)

        # Add locator marker
        ax3.scatter(float(stn_lon), float(stn_lat), color = 'red', s = 150, edgecolor = 'black', marker = '*', transform=ccrs.PlateCarree(), zorder = 100)

        print('Saving', Fig_dir + 'Meteogram_{}.png'.format(station))
        plt.savefig(save_path, bbox_inches = 'tight', dpi = 300)
        #plt.show()
        plt.close()
    except:
        pass
