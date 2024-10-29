# Michael Wasserstein
# Plot_Mesonet_Bias.py
# 10/10/2024
# Script takes in WRF outputs and plots them, as well as real mesonet data,
# downloaded from the synoptic API
# Script plots surface winds and temperatures

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_Mesonet_Bias.py -r 2 -p 2
# -r represents the run number you want to plot
# -p is the wrf path

import os, sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from map_script import *

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors, ListedColormap, LinearSegmentedColormap
import cartopy
import cartopy.crs as ccrs
from cartopy.io.img_tiles import GoogleTiles
import cartopy.feature as cf
from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
import pyart
import metpy.calc as mpcalc
from metpy.units import units

import numpy as np
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair)
import wrf
import glob
import pandas as pd
import datetime

import urllib.request as req
import os.path
import json
import urllib

######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")

args = parser.parse_args()

# # Get user inputs
run = str(args.run)
path = int(args.path)

run_number = '{}'.format(run).zfill(2) # leading zeros

# Time offset - do you want it or not?
time_offset = False
time_offset_minutes = -45 # This will add 45 minutes to the synoptic api time
time_offset_hours = 0

# ###########################
# ##### stuff for WRF  ######
# ###########################
if (path == 2) or (path == 6) or (path == 12): #TECPEC
    time_of_interest_list = pd.date_range(datetime.datetime(2019,3,22,12), datetime.datetime(2019,3,23,6), freq = '15min')
    if run in ['14', '16']:
        time_of_interest_list = pd.date_range(datetime.datetime(2019,3,22,12), datetime.datetime(2019,3,23,2), freq = '15min')
if path in [1,5,8,9]:
    time_of_interest_list = pd.date_range(datetime.datetime(2022,12,12,0), datetime.datetime(2022,12,15,0), freq = '15min')
    if run == '13':
        time_of_interest_list = pd.date_range(datetime.datetime(2022,12,13,0), datetime.datetime(2022,12,14,6), freq = '15min')

############################
## stuff for synoptic API ##
############################
API_ROOT = "https://api.synopticdata.com/v2/"
API_TOKEN = "0d0f87d395244930af38c1460e0c1a0f"
within = 30 # how many minutes within the time of interest do you want to get data for

# Base paths
base = '/uufs/chpc.utah.edu/common/home/'
home = base + 'u1371671/'
# paths for data
if path ==1:
    base_path = base + 'steenburgh-group12/michael/wrf/'
else:
    base_path = base + 'steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number) # Path of wrf data

# paths for saving fig
Fig_dir = home + 'WRF/Figures_{}/wrf_{}/Surface_Temperature_Bias/'.format(path,run_number)
if time_offset:
    Fig_dir = home + 'WRF/Figures_{}/wrf_{}/Surface_Temperature_Bias_time_offset_{}h_{}min/'.format(path,run_number,
                                                                                                   time_offset_hours, time_offset_minutes)
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)
    
for time_of_interest in time_of_interest_list:
    try:
        ###### Stuff for times ########
        time_of_interest_str = datetime.datetime.strftime(time_of_interest, '%Y-%m-%d_%H:%M:%S')
        time_of_interest_api = datetime.datetime.strftime(time_of_interest, '%Y%m%d%H%M') # Also to be used for saving
        time_of_interest_save = time_of_interest_api # make the variable names more clear
        time_of_interest_title = datetime.datetime.strftime(time_of_interest, '%b %-d, %Y %H:%M UTC')

        time_API_title = datetime.datetime.strftime(time_of_interest, '%b %-d, %Y %H:%M UTC')

        if time_offset:
            time_of_interest += datetime.timedelta(hours = time_offset_hours, minutes = time_offset_minutes)
            time_of_interest_str = datetime.datetime.strftime(time_of_interest, '%Y-%m-%d_%H:%M:%S')
            time_of_interest_save = datetime.datetime.strftime(time_of_interest, '%Y%m%d%H%M') # Also to be used for saving
            time_of_interest_title = datetime.datetime.strftime(time_of_interest, '%b %-d, %Y %H:%M UTC')


        ###########################
        ####### Data Access #######
        ###########################
        ####################################### WRF ##########################################
        f = WRF_path + 'wrfout_d03_{}'.format(time_of_interest_str) # get file name for the time of interest

        # Load in file
        wrf_file = Dataset(f)

        # Stuff for boundaries
        geobounds = wrf.geo_bounds(wrfin=wrf_file)
        bottom_latitude = geobounds.bottom_left.lat
        left_longitude = geobounds.bottom_left.lon
        top_latitude = geobounds.top_right.lat
        right_longitude = geobounds.top_right.lon

        # Get the initialization times
        init_time = wrf_file.SIMULATION_START_DATE
        init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
        init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')

        # Get the necessary wrf data
        ter = getvar(wrf_file, "ter", timeidx=-1)
        lati, long = latlon_coords(ter)
        T2 = getvar(wrf_file, 'T2') - 273.15
        u10 = getvar(wrf_file, 'uvmet10').sel(u_v = 'u')
        v10 = getvar(wrf_file, 'uvmet10').sel(u_v = 'v')
        #####################################################################################

        ####################################### Synoptic Data ##########################################
        # String to use for api request - all stations in requested domain
        bbox_string = f'{left_longitude},{bottom_latitude},{right_longitude},{top_latitude}'

        # Get data for the nearest time to a specified time
        api_request_url = os.path.join(API_ROOT, "stations/nearesttime")

        # Get data for the state of utah
        api_request_url += "?token={}&bbox={}&attime={}&within={}".format(API_TOKEN, bbox_string, time_of_interest_api,within) 

        response = req.urlopen(api_request_url) # Make request
        api_text_data = response.read() # read the api request

        # Load data as a dictionary
        use_data = json.loads(api_text_data) 

        # check the units of the data you're working with
        print('Temperature units:', use_data['UNITS']['air_temp'])
        print('Wind degree units:', use_data['UNITS']['wind_direction'])
        print('Wind speed units: ', use_data['UNITS']['wind_speed'])
        print('Rel Humid units: ', use_data['UNITS']['relative_humidity'])
        #################################################################################


        ###########################
        ####### Plotting ##########
        ###########################
        # Plot settings
        c_WRF = 'slategrey'
        c_obs = 'black'

        # COlormap for the temperature bias plot
        cmap = get_cmap('pyart_HomeyerRainbow')

        # Create a ListedColormap with 20 colors
        colors = [cmap(i) for i in np.linspace(0, 1, 6)]

        # change colors and make middle colors white
        colors.insert(3, (1,1,1,1))
        colors.insert(3, (1,1,1,1))
        new_cmap = ListedColormap(colors)

        ################## plot temperature #######################
        # Make Figure
        fig = plt.figure(figsize=(12,10))

        # Set projection and extent of plot
        ax = plt.axes(projection=ShadedReliefESRI().crs)
        ax.set_extent([left_longitude, right_longitude, bottom_latitude, top_latitude])

        # Add the map tiles as level 1
        ax.add_image(ShadedReliefESRI(), map_tile_detail)#, zorder=1)

        # Add north arrow
        ax.text(.05, .98 ,u'\u25B2 \nN ', ha='center', va='top', fontsize=16, family='Arial', rotation = 0, zorder=50, transform=ax.transAxes)

        # Empty lists to store data
        stn_longs = []
        stn_lats  = []
        stn_longs_T = []
        stn_lats_T  = []
        stn_temps = []
        stn_us     = []
        stn_vs     = []
        stn_temps_WRF = []

        # Loop through all the stations
        for i in range(len(use_data['STATION'])):
            # Extract station lat and lon
            stn_long = float(use_data['STATION'][i]['LONGITUDE'])
            stn_lati = float(use_data['STATION'][i]['LATITUDE'])

            # Extract WRF x and y
            x, y = wrf.ll_to_xy(wrf_file, stn_lati, stn_long)

            try:
                # Stuff for temps
                try: 
                    stn_temp = float(use_data['STATION'][i]['OBSERVATIONS']['air_temp_value_1']['value']) # The observed temperature
                    stn_temp_WRF = T2[y,x]                      # The WRF temperature corresponding to the location of the gridpoint

                    # Add the data to your lists
                    stn_longs_T.append(stn_long)        
                    stn_lats_T.append(stn_lati)
                    stn_temps.append(stn_temp)
                    stn_temps_WRF.append(stn_temp_WRF)

                    # Calculate the temperature difference from the lists, converting them to arrays
                    temperature_difference = np.array(stn_temps_WRF) - np.array(stn_temps)
                except:
                    pass

                try:
                    # Stuff for winds
                    stn_dir = use_data['STATION'][i]['OBSERVATIONS']['wind_direction_value_1']['value'] * units('degree') # Direction
                    stn_spd = use_data['STATION'][i]['OBSERVATIONS']['wind_speed_value_1']['value'] * units('meter_per_second') # Speed
                    stn_u, stn_v = mpcalc.wind_components(stn_spd, stn_dir) # Get the u and v components from the direction and speed

                    # Append data to lists to be used in plt.barbs
                    stn_longs.append(stn_long)
                    stn_lats.append(stn_lati)
                    stn_us.append(stn_u.magnitude)
                    stn_vs.append(stn_v.magnitude)

                except:
                    pass

            except:
                pass

        # Temperature settings
        min_temp = -8
        max_temp = 8

        # Plot temperature difference scatter from the wrf data
        sc = ax.scatter(stn_longs_T, stn_lats_T, c=temperature_difference,transform=ccrs.PlateCarree(), cmap =new_cmap, 
                          vmin = min_temp, vmax = max_temp,  alpha = 1, zorder = 30, edgecolors = c_obs)

        skip = 12 # Use this to display frequency of barbs
        # Plot wind barb from wrf data
        barb = ax.barbs(long.values[::skip,::skip], lati.values[::skip,::skip], u10.values[::skip,::skip], v10.values[::skip,::skip],
                        transform=ccrs.PlateCarree(), length=5, linewidth=0.5, zorder=25, color=c_WRF, alpha=1)

        # Plot wind barb from obs data
        ax.barbs(np.array(stn_longs), np.array(stn_lats), np.array(stn_us), np.array(stn_vs), 
                 transform=ccrs.PlateCarree(), length=7, linewidth=0.5, zorder=25, color=c_obs, alpha=1.0)

        # Add colorbar
        cax = plt.axes([0.93,0.13, 0.035, 0.723])
        plt.colorbar(sc, cax = cax, orientation = 'vertical', label = 'Model - Observation Temperature (°C)', ax = ax, fraction = 0.05,
                     ticks = np.arange(min_temp,max_temp+1,2), drawedges=True)



        # Add a title
        ax.set_title("WRF{} Run {}\nInit {}".format(path, run, init_time_str), fontsize = 10, loc = 'left',)
        ax.set_title("WRF Model Valid {}\nSynoptic API Valid {}".format(time_of_interest_title,time_API_title), fontsize = 10, loc = 'right',)


        # Save figure
        plt.savefig(Fig_dir + time_of_interest_save + '_WRF_obs_temperature_bias_wind.png', dpi = 300, bbox_inches = 'tight')
        #plt.show()
        plt.close()

    except:
        pass