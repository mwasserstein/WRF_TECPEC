# Michael Wasserstein
# 2_Domain_2_Panel_Temperature_Wind_WRF_ERA5.py
# 10/18/2024
# Script takes in WRF outputs and ERA5 data
# and plots 700 mbar temperature and winds

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/2_Domain_2_Panel_Temperature_Wind_WRF_ERA5.py -r 2 -p 2 -l 700
# -r represents the run number you want to plot
# -p is the wrf path (wrf1 or wrf2)
# -l is the level at which to plot (mbar)

import os, sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from Useful_Python_Functions import round_to_nearest_hour
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.gridspec import GridSpec
import matplotlib
from matplotlib.colors import from_levels_and_colors, ListedColormap, LinearSegmentedColormap
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature, COLORS
import pyart
import metpy.calc as mpcalc
from metpy.units import units
import math
from metpy.plots import USCOUNTIES
import xarray as xr
import numpy as np
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair,ALL_TIMES)
import wrf
import glob
import pandas as pd
import datetime
from scipy.ndimage import gaussian_filter
from datetime import datetime, timedelta

def find_time_idx(time, files, domain):
    '''
    Function will find the index in a list of files corresponding to the time you want
    '''
    
    # Loop through the file paths and extract the datetime part
    for idx, file_path in enumerate(data_files_d03):
        # Extract the datetime part from the file path
        date_str = file_path.split('_d0{}_'.format(domain))[1]
        file_datetime = datetime.strptime(date_str, '%Y-%m-%d_%H:%M:%S')

        # Compare with the target datetime
        if file_datetime == time:
            break # end the loop
            
    return idx

# ######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")
parser.add_argument("-l", "--level", help="Level at which to interpolate")

args = parser.parse_args()

###########################
##### stuff for WRF  ######
###########################
# Get user inputs
run = str(args.run)
path = int(args.path)
level_to_interp = int(args.level)
print('Plotting data for run', run,)

run_number = '{}'.format(run).zfill(2) # leading zeros

plot_d04 = False # plot data for domain 4?

if (path == 2) or (path == 6) or (path == 12): #TECPEC
    time_of_interest_list = pd.date_range(datetime(2019,3,22,12), datetime(2019,3,23,6), freq = '15min')
    ERA5_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group8/era5/iso/201903/' # path containing era5 data
    # Start and end time for analysis period
    start_time_analysis = datetime(2019,3,22,19,30)
    end_time_analysis = datetime(2019,3,23,0,15) 
    if run in ['14']:
        time_of_interest_list = pd.date_range(datetime(2019,3,22,12), datetime(2019,3,23,2), freq = '15min')
    if run in ['15']:
        time_of_interest_list = pd.date_range(datetime(2019,3,22,0), datetime(2019,3,23,2), freq = '15min')
if (path == 1) or (path == 5) or (path == 8) or (path == 9):
    time_of_interest_list = pd.date_range(datetime(2022,12,13,0), datetime(2022,12,14,6), freq = '1H')
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
if plot_d04:
    Fig_dir = home + 'WRF/Figures_{}/wrf_{}/2_Panel_ERA5_700_mb_d03_d04_Temp_Wind/'.format(path,run_number)
else:
    Fig_dir = home + 'WRF/Figures_{}/wrf_{}/2_Panel_ERA5_700_mb_d03_Temp_Wind/'.format(path,run_number)

# Make fig dir if it doesn't exist
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

# load in all the wrf output data files
data_files_d04 = glob.glob(WRF_path  + '*wrfout_d04*') # for the outermost domain
data_files_d04.sort()

# load in all the wrf output data files
data_files_d03 = glob.glob(WRF_path  + '*wrfout_d03*') # for the outermost domain
data_files_d03.sort()

# Find the start and end indicies for your period of analysis - only necessary for d03 cause it has same timing as d04
start_ind = find_time_idx(time=start_time_analysis, files=data_files_d03, domain=3)
end_ind = find_time_idx(time=end_time_analysis, files=data_files_d03, domain=3)

# Hours for time offset from the wrf time of intereset
time_delta_hours = [-2,-1,0,1,2,]

# Levels for plotting temperature
T_levels = np.arange(-8,2.1,0.25)

# Colormap to use
cmap = plt.get_cmap('pyart_HomeyerRainbow')

# Colors corresponding to the number of levels in your temperature
colors = []
for v in np.linspace(0,1,len(T_levels)):
    colors.append(cmap(v))

# Loop through times of interest, only extracting from the start index to the end
for i, time_of_interest in enumerate(time_of_interest_list[start_ind:end_ind]):
    i += start_ind # adjustment to i based on what your start index is

    # Time stuff
    time_nearest_hour = round_to_nearest_hour(time_of_interest) # round your time of interest to the nearest hour (for ERA5 since its on the hour)
    time_of_interest_WRF = time_of_interest.strftime('%Y-%m-%d_%H:%M:%S')
    time_of_interest_title = datetime.strftime(time_of_interest, '%b %-d, %Y %H:%M UTC')
    time_of_interest_save = time_of_interest.strftime('%Y%m%d%H%M%S')
    
    # WRF data files
    ds_d03 = Dataset(data_files_d03[i])
    if plot_d04:
        ds_d04 = Dataset(data_files_d04[i])

    #################################
    ########### DOmain 3 ############
    #################################
    # Extract data
    ht_d03 = getvar(ds_d03, "z", timeidx=-1)
    u_d03 =  getvar(ds_d03, "ua", timeidx = -1)
    v_d03 =  getvar(ds_d03, "va",  timeidx = -1)
    lati_d03, long_d03 = latlon_coords(ht_d03)
    p_d03 = getvar(ds_d03, 'pressure',  timeidx = -1)
    tc_d03 = getvar(ds_d03, 'tc',  timeidx = -1)
    
    # Stuff for boundaries
    geobounds = wrf.geo_bounds(wrfin=ds_d03)
    bottom_latitude_d03 = geobounds.bottom_left.lat
    left_longitude_d03 = geobounds.bottom_left.lon
    top_latitude_d03 = geobounds.top_right.lat
    right_longitude_d03 = geobounds.top_right.lon
    
    # Interpolate to level of interest 
    u_lev_d03 = wrf.interplevel(field3d=u_d03, vert=p_d03, desiredlev=level_to_interp,)
    v_lev_d03 = wrf.interplevel(field3d=v_d03, vert=p_d03, desiredlev=level_to_interp,)
    tc_lev_d03 = wrf.interplevel(field3d=tc_d03, vert=p_d03, desiredlev=level_to_interp,)

    if plot_d04:
        #################################
        ########### DOmain 4 ############
        #################################
        # Extract data
        ht_d04 = getvar(ds_d04, "z", timeidx=-1)
        u_d04 =  getvar(ds_d04, "ua", timeidx = -1)
        v_d04 =  getvar(ds_d04, "va",  timeidx = -1)
        lati_d04, long_d04 = latlon_coords(ht_d04)
        p_d04 = getvar(ds_d04, 'pressure',  timeidx = -1)
        tc_d04 = getvar(ds_d04, 'tc',  timeidx = -1)

        #Stuff for boundaries
        geobounds = wrf.geo_bounds(wrfin=ds_d04)
        bottom_latitude_d04 = geobounds.bottom_left.lat
        left_longitude_d04 = geobounds.bottom_left.lon
        top_latitude_d04 = geobounds.top_right.lat
        right_longitude_d04 = geobounds.top_right.lon

        # Interpolate to level of interest 
        u_lev_d04 = wrf.interplevel(field3d=u_d04, vert=p_d04, desiredlev=level_to_interp,)
        v_lev_d04 = wrf.interplevel(field3d=v_d04, vert=p_d04, desiredlev=level_to_interp,)
        tc_lev_d04 = wrf.interplevel(field3d=tc_d04, vert=p_d04, desiredlev=level_to_interp,)

    for delta_hour in time_delta_hours: # loop through all time offset hours
        
        # Stuff for era5 times
        time_of_interest_ERA5_dt = time_nearest_hour + timedelta(hours=delta_hour)
        time_of_interest_ERA5 = time_of_interest_ERA5_dt.strftime('%Y%m%d00_%Y%m%d23')
        time_of_interest_ERA5_title = datetime.strftime(time_of_interest_ERA5_dt, '%b %-d, %Y %H:%M UTC')
        time_of_interest_ERA5_save = time_of_interest_ERA5_dt.strftime('%Y%m%d%H%M%S')
        
        # Print times to help user
        print('ERA5 time:', time_of_interest_ERA5_title)
        print('WRF time: ', time_of_interest_title,)
        
        ################################ ERA5 ##################################
        # data paths for t, u, v
        ERA5_t_file = ERA5_path + 'e5.oper.an.pl.128_130_t.ll025sc.' + time_of_interest_ERA5 + '.WE.nc'
        ERA5_v_file = ERA5_path + 'e5.oper.an.pl.128_132_v.ll025uv.' + time_of_interest_ERA5 + '.WE.nc'
        ERA5_u_file = ERA5_path + 'e5.oper.an.pl.128_131_u.ll025uv.' + time_of_interest_ERA5 + '.WE.nc'

        # combine files and open using xarray
        ERA5_files = [ERA5_u_file, ERA5_v_file, ERA5_t_file]
        ERA5_ds = xr.open_mfdataset(ERA5_files)

        # Extract the data
        ERA5_u = ERA5_ds.sel(level = level_to_interp, time = time_of_interest_ERA5_dt).u.values
        ERA5_v = ERA5_ds.sel(level = level_to_interp, time = time_of_interest_ERA5_dt).v.values
        ERA5_t = ERA5_ds.sel(level = level_to_interp, time = time_of_interest_ERA5_dt).t.values - 273.15 # Convert to celcius
        ERA5_longitudes = ERA5_ds.longitude.values
        ERA5_latitudes = ERA5_ds.latitude.values

        #########################################################################
        ######################### Plotting ######################################
        #########################################################################
        # Create figure
        fig, ((ax1, ax2)) = plt.subplots(1,2,figsize = (14, 14), facecolor = 'white', edgecolor = 'k',
                                  subplot_kw = {'projection' : ccrs.PlateCarree()})

        ####################################################################################
        ############################### ax1 ERA5 - 700 mbar ################################
        ####################################################################################
        # Set axis extent to correspond with WRF d03
        ax1.set_extent([left_longitude_d03, right_longitude_d03, bottom_latitude_d03, top_latitude_d03], crs=ccrs.PlateCarree())

        # Plot Temperature
        T_Plot = ax1.contourf(ERA5_longitudes, ERA5_latitudes, ERA5_t, transform = ccrs.PlateCarree(),levels = T_levels,
                colors = colors, extend = 'both', zorder=1)

        # Add the wind barbs
        ax1.barbs(ERA5_longitudes, ERA5_latitudes, ERA5_u, ERA5_v,
            length=4, zorder=3, color = 'black',regrid_shape = 15)

        if plot_d04:
            # Add a bounding box indicating where d04 from WRF is 
            ax1.add_patch(matplotlib.patches.Rectangle((left_longitude_d04, bottom_latitude_d04), right_longitude_d04-left_longitude_d04, top_latitude_d04-bottom_latitude_d04,
                     fill=None, lw=2, edgecolor='grey', zorder=50,transform = ccrs.PlateCarree(), alpha = 0.5))

        # Add U.S. state borders and counties
        ax1.add_feature(cfeature.STATES, edgecolor='white', linewidth=0.5, zorder = 10)
        ax1.add_feature(USCOUNTIES.with_scale('500k'), edgecolor = 'white', zorder = 51)

        # Add titles
        ax1.set_title('700-mb Temperature\nand Winds', loc = 'right')
        ax1.set_title('ERA5 Valid {}'.format(time_of_interest_ERA5_title), loc = 'left')

        ####################################################################################
        ############################### ax2 WRF - 700 mbar #################################
        ####################################################################################
        ####### DOmain 3
        
        # Set extent for d03 boundaries
        ax2.set_extent([left_longitude_d03, right_longitude_d03, bottom_latitude_d03, top_latitude_d03], crs=ccrs.PlateCarree())
        
        # Plot temperature for 700 mbar
        T_Plot = ax2.contourf(long_d03, lati_d03, tc_lev_d03, transform = ccrs.PlateCarree(),levels = T_levels,
                    colors = colors, extend = 'both', zorder=1)

        # Add the wind barbs
        ax2.barbs(to_np(long_d03[::,::]), to_np(lati_d03[::,::]),
              to_np(u_lev_d03[::, ::]), to_np(v_lev_d03[::, ::]),
            length=4, zorder=3, color = 'black',regrid_shape = 15)

        ###### DOmain 4 
        if plot_d04:
            # Plot temperature
            T_Plot = ax2.contourf(long_d04, lati_d04, tc_lev_d04, transform = ccrs.PlateCarree(),levels = T_levels,
                        colors = colors, extend = 'both', zorder = 9)

            # Add the wind barbs
            ax2.barbs(to_np(long_d04[::,::]), to_np(lati_d04[::,::]),
                  to_np(u_lev_d04[::, ::]), to_np(v_lev_d04[::, ::]),
                length=4,  color = 'black', zorder = 11, regrid_shape = 15)

            # bounding box for d04
            ax2.add_patch(matplotlib.patches.Rectangle((left_longitude_d04, bottom_latitude_d04), right_longitude_d04-left_longitude_d04, top_latitude_d04-bottom_latitude_d04,
                     fill=None, lw=2, edgecolor='grey', zorder=50,transform = ccrs.PlateCarree(), alpha = 0.5))

        # Add U.S. state borders and counties
        ax2.add_feature(cfeature.STATES, edgecolor='white', linewidth=0.5, zorder = 51)
        ax2.add_feature(USCOUNTIES.with_scale('500k'), edgecolor = 'white', zorder = 51)

        # Set titles
        if plot_d04:
            ax2.set_title(f'Domains 3-4\n{level_to_interp}-mb Temperature\nand Winds', loc = 'right')
        else:
            ax2.set_title(f'Domain 3\n{level_to_interp}-mb Temperature\nand Winds', loc = 'right')
        ax2.set_title('WRF{} Run {}\nValid {}'.format(path, run, time_of_interest_title), loc = 'left')

        # Set a background color
        ax2.set_facecolor('gainsboro')

        # Add the color bar
        cax = plt.axes([0.16,0.33, 0.7, 0.02])
        cb_wa = fig.colorbar(T_Plot, cax = cax, orientation = 'horizontal', label = 'Temperature (°C)', fraction = 0.05,
                     ticks = T_levels[::4])
        cb_wa.ax.tick_params(labelsize=10)

        # Adjust subplot positioning
        plt.subplots_adjust(wspace = 0.1)

        # Save, show, and close
        plt.savefig(Fig_dir + 'WRF_{}_E5_{}_{}_Wind_Temp.png'.format(time_of_interest_save, time_of_interest_ERA5_save, level_to_interp), dpi = 300, bbox_inches = 'tight')
        #plt.show()
        plt.close()