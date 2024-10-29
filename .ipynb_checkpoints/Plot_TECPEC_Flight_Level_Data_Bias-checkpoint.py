# Michael Wasserstein
# Plot_TECPEC_Flight_Level_Data_Bias.py
# 10/18/2024
# Script takes TECPEC flight level data
# And WRF data and plots a flight level bias for various variables for the data

####### Usage #########
# Conda environment - Radar_env_3
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_TECPEC_Flight_Level_Data_Bias.py -r 2 -p 2
# -r represents the run number you want to plot
# -p is the wrf path (wrf1 or wrf2)

import numpy as np
import pandas as pd
import os, sys
import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import SkewT
import matplotlib.pyplot as plt
#import datetime
from datetime import datetime, timedelta
import datetime

import glob
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair, ALL_TIMES)
import wrf
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from Useful_Python_Functions import truncate_colormap
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors
import matplotlib.patheffects as path_effects
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as mimage
from matplotlib.gridspec import GridSpec
import cartopy
import cartopy.crs as ccrs
from cartopy.io.img_tiles import GoogleTiles
import cartopy.feature as cf
import pandas as pd
import xarray as xr
import pyproj

######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")

args = parser.parse_args()

##########################
#### stuff for WRF  ######
##########################
#Get user inputs
run = str(args.run)
path = int(args.path)
print('Plotting data for run', run,)

run_number = '{}'.format(run).zfill(2) # leading zeros

time_offset = False # DO you want a time offset
time_offset_minutes = -45 # This will be added to the TECPEC flight level data time
time_offset_hours = 0
    
############################## Plotting Stuff #########################
# Extent of plot
left_lon = -112.35
right_lon = -111.4
bot_lat = 40.3
top_lat = 40.8

# Get colormap for terrain and truncate it
cmap = plt.get_cmap('terrain')
new_cmap = truncate_colormap(cmap, 0.23, 1)

# Path where you have TECPEC flight level data
data_path = '/uufs/chpc.utah.edu/common/home/u1371671/steenburgh-group12/michael/TECPEC_Flight_level_data/'

# Load in TECPEC flight leg data
flight_df = pd.read_csv(data_path + 'TECPEC_leg_data.csv')
legs = flight_df['Leg'].values
start_lats = flight_df['Start_lat'].values
end_lats = flight_df['End_lat'].values
start_lons = flight_df['Start_lon'].values
end_lons = flight_df['End_lon'].values

# Extract flight latitudes and longitudes (organized by leg and 200 lats and lons)
latitudes = np.load(data_path + 'latitudes.npy')
longitudes = np.load(data_path + 'longitudes.npy')


# From daves email:
# dayt                    time in Matlab datenumber format (days since the year 0; dayt=1 is 1 January 0000)
# latc                      latitutde in deg
# lonc                     longitude in deg
# zmsl                     altitude above mean sea level in meters
# tamb                   temperature in deg C
# licor_mr              water vapor mixing ratio in g/kg from Licor probe (generally best option)
# licor_tdp             dewpoint temperature in deg C from Licor probe (generally best option)
# pmb                    static pressure in mb
# uwind                 zonal wind component in m/s
# vwind                  meridional wind component in m/s
# wwind                 vertical wind component in m/s
# head                    aircraft heading in deg from N
# track                    aircraft track in deg from N

# Note for michael
#datearray is the times you want to use # This was converted in matlab

# Load in all arrays for data
# We flatten them because they are arrays of shape (1,XXX)
#dayt = pd.read_csv(data_path+'dayt.csv',header = None).values.flatten()
#time = pd.read_csv(data_path+'time.csv',header = None).values.flatten()
latc = pd.read_csv(data_path+'latc.csv',header = None).values.flatten()
lonc = pd.read_csv(data_path+'lonc.csv',header = None).values.flatten()
zmsl = pd.read_csv(data_path+'zmsl.csv',header = None).values.flatten()
tamb = pd.read_csv(data_path+'tamb.csv',header = None).values.flatten()
licor_mr = pd.read_csv(data_path+'licor_mr.csv',header = None).values.flatten()
licor_rh = pd.read_csv(data_path+'licor_rh.csv',header = None).values.flatten()
licor_tdp = pd.read_csv(data_path+'licor_tdp.csv',header = None).values.flatten()
rh = pd.read_csv(data_path+'rh.csv',header = None).values.flatten()
mr = pd.read_csv(data_path+'rh.csv',header = None).values.flatten()
tdp = pd.read_csv(data_path+'tdp.csv',header = None).values.flatten()
pmb = pd.read_csv(data_path+'pmb.csv',header = None).values.flatten()
uwind = pd.read_csv(data_path+'uwind.csv',header = None).values.flatten()
vwind = pd.read_csv(data_path+'vwind.csv',header = None).values.flatten()
wwind = pd.read_csv(data_path+'wwind.csv',header = None).values.flatten()
head = pd.read_csv(data_path+'head.csv',header = None).values.flatten()
track = pd.read_csv(data_path+'track.csv',header = None).values.flatten()
datearray = pd.read_csv(data_path+'datearray.csv',header = None).values.flatten() # This was converted in matlab

# Convert times to pd datetime
datearray = pd.to_datetime(datearray)

# Generate a dataframe
flight_data = pd.DataFrame({
    'datearray' : datearray,
    'latc' : latc,
    'lonc' : lonc,
    'zmsl' : zmsl,
    'tamb' : tamb,
    'licor_mr' : licor_mr,
    'licor_rh' : licor_rh,
    'licor_tdp' : licor_tdp,
    'rh' : rh, 
    'mr' : mr,
    'tdp' : tdp,
    'pmb' : pmb,
    'uwind' : uwind,
    'vwind' : vwind,
    'wwind' : wwind,
    'head' : head,
    'track' : track,
             })

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
Fig_dir = home + 'WRF' + '/Figures_{}/wrf_{}/Flight_Level/'.format(path,run_number)
if time_offset:
    Fig_dir = home + 'WRF' + '/Figures_{}/wrf_{}/Flight_Level_time_offset/'.format(path,run_number)

# Make fig dir if it doesn't already exist
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

###########################
####### Data Access #######
###########################
# list of the wrf output data files
data_files_d03 = glob.glob(WRF_path  + '*wrfout_d03*') 
data_files_d03.sort()

# Load in the datsets
wrflist = [Dataset(file) for file in data_files_d03]

# Get the initialization times
init_time = wrflist[0].SIMULATION_START_DATE
init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')

# Boundaries
geobounds = wrf.geo_bounds(wrfin=wrflist[0])
min_lat = geobounds.bottom_left.lat
min_lon = geobounds.bottom_left.lon
max_lat = geobounds.top_right.lat
max_lon = geobounds.top_right.lon

# Get all the WRF data for all times you need
tamb_WRF = getvar(wrflist, "temp", timeidx=ALL_TIMES, method="cat", units = 'degC')
zmsl_WRF = wrf.g_geoht.get_height(wrflist, timeidx=ALL_TIMES,method="cat")
tdp_WRF = getvar(wrflist,"td", timeidx=ALL_TIMES, method="cat")
mr_WRF = getvar(wrflist,"QVAPOR", timeidx=ALL_TIMES, method="cat") * 1000 # Convert to g/kg
u_WRF = getvar(wrflist,'ua' , timeidx=ALL_TIMES, method="cat")
v_WRF = getvar(wrflist,'va',timeidx=ALL_TIMES, method="cat")
p_WRF = getvar(wrflist,'p',timeidx=ALL_TIMES, method="cat")
ter_WRF = getvar(wrflist, 'ter', timeidx = 0,)
lake_WRF = getvar(wrflist, 'LANDMASK', timeidx = 0)

# Loop through all the flight legs
for leg in legs[:]:
    
    # Extract flight data
    leg_data = flight_df[flight_df['Leg'] == leg]
    
    # Shorten the line to fall within geobounds (not totally sure this is ncecssary, but it can help
    indices_within_bounds = np.where((latitudes[leg-1] >= min_lat) & (latitudes[leg-1] <= max_lat) & (longitudes[leg-1] >= min_lon) & (longitudes[leg-1] <= max_lon))[0]

    # Get all the shortened list of longitudes and latitudes for the leg of interest
    shortened_longitudes = longitudes[leg-1][indices_within_bounds]
    shortened_latitudes = latitudes[leg-1][indices_within_bounds]
    
    # Change the dataframe based on the shorted latitudes and longitudes
    leg_data['End_lat'] = shortened_latitudes[-1]
    leg_data['Start_lat'] = shortened_latitudes[0]
    leg_data['End_lon'] = shortened_longitudes[-1]
    leg_data['Start_lon'] = shortened_longitudes[0]
    
    # Extract the start and end times for the flight leg
    start_time, end_time = pd.to_datetime(leg_data['Start_dt'].values), pd.to_datetime(leg_data['End_dt'].values)
    
    print(len(start_time))
    print(len(flight_data['datearray'].values))
    
    # Subset the leg data for the correct times
    leg_data = flight_data[(flight_data['datearray'].values >= start_time.values) & (flight_data['datearray'].values <= end_time.values)].reset_index(drop = True)
    
    if leg % 2 != 0: # Check if the leg number is odd If the remainder in the division is not zero
        # Subset the data so that all latitude and longitude points in the data are within your wrf model grid
        leg_data = leg_data[(leg_data['lonc'] >= shortened_longitudes[0]) & (leg_data['lonc'] <= shortened_longitudes[-1]) &
             (leg_data['latc'] >= shortened_latitudes[0]) & (leg_data['latc'] <= shortened_latitudes[-1])].reset_index(drop = True)
    else: # If the leg is even
        leg_data = leg_data[(leg_data['lonc'] <= shortened_longitudes[0]) & (leg_data['lonc'] >= shortened_longitudes[-1]) &
             (leg_data['latc'] <= shortened_latitudes[0]) & (leg_data['latc'] >= shortened_latitudes[-1])].reset_index(drop = True)
    
    # Get the data
    leg_lons = leg_data.lonc.values
    leg_lats = leg_data.latc.values
    leg_tamb = leg_data.tamb.values
    
    # the flight leg start and end time
    start_time_str = pd.to_datetime(start_time).strftime('%b %-d, %Y %H:%M').values[0]
    end_time_str = pd.to_datetime(end_time).strftime('%b %-d, %Y %H:%M').values[0]
    
    if time_offset:
        # offset times for the wrf
        start_time_str = (pd.to_datetime(start_time) + datetime.timedelta(minutes=time_offset_minutes,
                                                                            hours = time_offset_hours)).strftime('%b %-d, %Y %H:%M').values[0]
        end_time_str = (pd.to_datetime(end_time) + datetime.timedelta(minutes=time_offset_minutes,
                                                                        hours = time_offset_hours)).strftime('%b %-d, %Y %H:%M').values[0]

    print(f'Flight leg {leg} starting {start_time_str} ending {end_time_str}')

    # Initialize an empty array to store the data
    temperature_WRF = np.empty(len(leg_data))

    # Loop through all the flight leg data
    for index, row in leg_data.iterrows():
        # Find the time of the flight for the given row. THis will enable you to select that time in the wrf data
        if time_offset:
            # Adjust time for time offset
            time = pd.to_datetime(row.datearray) + datetime.timedelta(minutes=time_offset_minutes,hours = time_offset_hours)
        else:
            # Extract time from leg data
            time = pd.to_datetime(row.datearray)
        
        # latitude, longitude, and height of the flight at the given time
        latc = row.latc
        lonc = row.lonc
        zmsl = row.zmsl

        # Get the data at the time of interest
        tamb_time = tamb_WRF.sel(Time = time, method = 'nearest')
        zmsl_time = zmsl_WRF.sel(Time = time, method = 'nearest')

        # Get the wrf x and y from the lon and lat
        x_y = wrf.ll_to_xy(wrfin=wrflist,latitude=latc, longitude=lonc) # zeroth index give west_east, first index gives south_norht
        south_north = x_y[1].values
        west_east = x_y[0].values

        # In the first line, you interpolate the variable of interest to the level of the flight
        # In the next, you select the value
        tamb_interpolated_to_level = wrf.interplevel(field3d=tamb_time, vert = zmsl_time, desiredlev = zmsl) # Get variable at the level of the flight
        tamb_at_location = tamb_interpolated_to_level.sel(south_north = south_north, west_east = west_east,).values # Get variable at the nearest x and y gridpoints
        temperature_WRF[index] = tamb_at_location        
    
    # Calculate the temperature bias (model - obs)
    temperature_bias = temperature_WRF - leg_tamb

    
    ################################# Ambient temperature Plot ############################################
    # Create figure
    fig = plt.figure(constrained_layout=True, figsize = (8, 13), 
                 facecolor = 'white', edgecolor = 'k', )
    gs = GridSpec(1, 1, figure=fig,)

    # create sub plots as grid
    ax = fig.add_subplot(gs[0,:], projection = ccrs.PlateCarree())

    # Set extent
    ax.set_extent([left_lon, right_lon, bot_lat, top_lat])

    # Plot temperature bias data
    path = ax.scatter(leg_lons, leg_lats, c = temperature_bias, cmap = 'bwr', 
                    s = 12, transform=ccrs.PlateCarree(), zorder = 100, vmin = -2, vmax = 2)

    # Plot terrain data
    topo_fill_zoom = ax.contourf(ter_WRF.XLONG,ter_WRF.XLAT,ter_WRF.values,levels = np.arange(1000, 3001, 100),
                extend='max', alpha=1,antialiased=True,cmap = new_cmap,zorder=2,transform=ccrs.PlateCarree())

    # Add lakes
    lake = ax.contourf(lake_WRF.XLONG, lake_WRF.XLAT, lake_WRF.values, levels = np.arange(-0.1,1.0), colors = 'midnightblue',
                        alpha = 1, zorder = 3, transform = ccrs.PlateCarree())

    # Plot terrain colorbar
    cax = plt.axes([0.12,0.23, 0.8, 0.012])
    plt.colorbar(topo_fill_zoom, cax = cax, fraction = 0.5, label = 'Elevation (m)', orientation = 'horizontal')

    # Plot colorbar for temperature
    cax = plt.axes([0.12,0.13, 0.8, 0.012])
    plt.colorbar(path, cax = cax, fraction = 0.5, label = 'Temperature Bias (WRF - Observation °C)', orientation = 'horizontal')

    # Titles
    ax.set_title('TECPEC leg {}\n{}-{} UTC'.format(leg, start_time_str, end_time_str), loc = 'left')
    ax.set_title('WRF run {} - d03'.format(run_number) + "\nInit "+init_time_str + '\n{}-{} UTC'.format(start_time_str,
                                                                                                   end_time_str), loc = 'right')
    # Save, show and close
    plt.savefig(Fig_dir + f'Temperature_Bias_leg_{leg}.png', dpi = 300, bbox_inches = 'tight')
    #plt.show()
    plt.close()