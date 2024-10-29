# Michael Wasserstein
# 2_Domain_2_Panel_Temperature_Wind_HRRR_WRF.ipynb
# 10/18/2024
# Script takes in WRF outputs as well as HRRR data
# And plots a two panel of 700 mbar temperature and wind

####### Usage #########
# Conda environment - xesmf_env
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/2_Domain_2_Panel_Temperature_Wind_HRRR_WRF.py -r 2 -p 2
# -r represents the run number you want to plot
# -p is the wrf path (wrf1 or wrf2)

import os, sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from Useful_Python_Functions import round_to_nearest_hour
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (NullFormatter, ScalarFormatter)
import matplotlib
from matplotlib.colors import from_levels_and_colors, ListedColormap, LinearSegmentedColormap
import cartopy
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, COLORS
import pyart
import metpy.calc as mpcalc
from metpy.units import units
import math
from metpy.plots import USCOUNTIES
import numpy as np
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair,ALL_TIMES)
import wrf
import glob
import pandas as pd
import datetime
from scipy.ndimage import gaussian_filter

from herbie import Herbie # Package for HRRR


def find_time_idx(time, files, domain):
    '''
    Function will find the index in a list of files corresponding to the time you want
    '''
    
    # Loop through the file paths and extract the datetime part
    for idx, file_path in enumerate(data_files_d03):
        # Extract the datetime part from the file path
        date_str = file_path.split('_d0{}_'.format(domain))[1]
        file_datetime = datetime.datetime.strptime(date_str, '%Y-%m-%d_%H:%M:%S')

        # Compare with the target datetime
        if file_datetime == time:
            break # if the time is the same, end the loop
            
    return idx



######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")
parser.add_argument("-l", "--level", help="Level to plot in mbar")

args = parser.parse_args()

###########################
##### stuff for WRF  ######
###########################
# Get user inputs
run = str(args.run)
path = int(args.path)
level_to_interp = int(args.level) #mbar

run_number = '{}'.format(run).zfill(2) # Leading zeros

# Do you want to plot d04 data
plot_d04 = False

if (path == 2) or (path == 6) or (path == 12): #TECPEC
    # TImes for entire run
    time_of_interest_list = pd.date_range(datetime.datetime(2019,3,22,12), datetime.datetime(2019,3,23,6), freq = '15min')
    # Start and end time for the period you want to analyze
    start_time_analysis = datetime.datetime(2019,3,22,19,30)
    end_time_analysis = datetime.datetime(2019,3,23,0,15) 
    
    if run in ['14']:
        # TImes for entire run
        time_of_interest_list = pd.date_range(datetime.datetime(2019,3,22,12), datetime.datetime(2019,3,23,2), freq = '15min')
        
    if run in ['15']:
        # TImes for entire run
        time_of_interest_list = pd.date_range(datetime.datetime(2019,3,22,0), datetime.datetime(2019,3,23,2), freq = '15min')   
    
if (path == 1) or (path == 5) or (path == 8) or (path == 9):
    time_of_interest_list = pd.date_range(datetime.datetime(2022,12,13,0), datetime.datetime(2022,12,14,6), freq = '1H')

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
    Fig_dir = home + 'WRF/Figures_{}/wrf_{}/2_Panel_HRRR_700_mb_d03_d04_Temp_Wind/'.format(path,run_number)
else:
    Fig_dir = home + 'WRF/Figures_{}/wrf_{}/2_Panel_HRRR_700_mb_d03_Temp_Wind/'.format(path,run_number)

if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)
    
# load in all the wrf output data files
data_files_d04 = glob.glob(WRF_path  + '*wrfout_d04*') # for the outermost domain
data_files_d04.sort()

# load in all the wrf output data files
data_files_d03 = glob.glob(WRF_path  + '*wrfout_d03*') # for the outermost domain
data_files_d03.sort()

# Find start and end indicies in the files corresponding to your times you want to analyze
start_ind = find_time_idx(time=start_time_analysis, files=data_files_d03, domain=3)
end_ind = find_time_idx(time=end_time_analysis, files=data_files_d03, domain=3)

# Offset times from time for the HRRR model
time_delta_hours = [-2,-1,0,1,2,]

################################### Plotting ###################################
# Levels for temperature plot
T_levels = np.arange(-8,2.1,0.25)

# COlormap for temperature
cmap = plt.get_cmap('pyart_HomeyerRainbow')

# Get a list of colors based on the number of levels you're using
colors = []
for v in np.linspace(0,1,len(T_levels)):
    colors.append(cmap(v))
    
# Loop through times of interest, only extracting from the start index to the end
for i, time_of_interest in enumerate(time_of_interest_list[start_ind:end_ind]):
    i += start_ind # adjustment to i based on what your start index is

    # Stuff for times
    time_nearest_hour = round_to_nearest_hour(time_of_interest) # round your time of interest to the nearest hour (for HRRR since its on the hour)
    time_of_interest_WRF = time_of_interest.strftime('%Y-%m-%d_%H:%M:%S')
    time_of_interest_title = datetime.datetime.strftime(time_of_interest, '%b %-d, %Y %H:%M UTC')
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
    
    ###############################################################
    ############################# HRRR ############################
    ###############################################################
    for delta_hour in time_delta_hours: # loop through the time offsets that you've specified
        # Times for the HRRR data
        time_of_interest_HRRR_dt = time_nearest_hour + datetime.timedelta(hours=delta_hour) # TIme based on time offset
        time_of_interest_HRRR = time_of_interest_HRRR_dt.strftime('%Y-%m-%d %H:%M')
        time_of_interest_HRRR_title = datetime.datetime.strftime(time_of_interest_HRRR_dt, '%b %-d, %Y %H:%M UTC')
        time_of_interest_HRRR_save = time_of_interest_HRRR_dt.strftime('%Y%m%d%H%M%S')
        
        # Print times to help user
        print('HRRR time:', time_of_interest_HRRR_title)
        print('WRF time :', time_of_interest_title,)

        # Create Herbie object for the HRRR model 6-hr surface forecast product
        H = Herbie(
          time_of_interest_HRRR,
          model='hrrr',
          product='prs', # prs or sfc (sometimes sfc works even for pressure data, but prs is better I think)
          fxx=0,
          overwrite = False # Use local copy of file if it already exists if set to false
        )

        # Get a u, v wind, and temperature dataset for 700 mbar
        ds_u_v = H.xarray(":[U\|V]GRD:700 mb")
        ds_T = H.xarray(":TMP:700 mb")


        ############################### Create Figure ################################
        fig = plt.figure(figsize = (14, 14), facecolor = 'white', edgecolor = 'k',)

        nrows = 20
        ncols = 24

        ############################### ax1 HRRR ################################
        # Create axis using subplot2grid
        ax1 = plt.subplot2grid((nrows, ncols), (0, 0), rowspan = 20, colspan = 10, projection = ccrs.PlateCarree())

        # Set extent and make it the same as the wrf d03 boundaries
        ax1.set_extent([left_longitude_d03, right_longitude_d03, bottom_latitude_d03, top_latitude_d03], crs=ccrs.PlateCarree())

        # Make longitudes go from -180 to 180
        longitude_u_v = np.where(ds_u_v.longitude.values > 180, ds_u_v.longitude.values - 360, ds_u_v.longitude.values)

        # Add the wind barbs
        ax1.barbs(longitude_u_v, ds_u_v.latitude.values,
              ds_u_v.u.values[::,::], ds_u_v.v.values[::,::],
            length=4, zorder=100, color = 'black',regrid_shape = 15)

        # Plot temperature, adjusting longitude to be -180 to 180 and adjusting temeprature
        ax1.contourf(ds_T.longitude-360, ds_T.latitude, ds_T.t-273.15, transform = ccrs.PlateCarree(),levels = T_levels,
                    colors = colors, zorder = 98, extend = 'both')


        # Add U.S. state borders
        ax1.add_feature(cfeature.STATES, edgecolor='white', linewidth=0.5, zorder = 100)
        ax1.add_feature(USCOUNTIES.with_scale('500k'), edgecolor = 'white', zorder = 100)

        # Plot title
        ax1.set_title('HRRR FCST h 0 \nValid {}'.format( time_of_interest_HRRR_title), loc = 'left')

        if plot_d04:
            # d04 box
            ax1.add_patch(matplotlib.patches.Rectangle((left_longitude_d04, bottom_latitude_d04), right_longitude_d04-left_longitude_d04, top_latitude_d04-bottom_latitude_d04,
                     fill=None, lw=2, edgecolor='grey', zorder=2000,transform = ccrs.PlateCarree(), alpha = 0.5))

        ############################### ax2 WRF ################################
        # Create axis using subplot2grid
        ax2 = plt.subplot2grid((nrows, ncols), (1, 14), rowspan = 18, colspan = 12, projection = ccrs.PlateCarree())

        ##### DOmain 3
        T_Plot = ax2.contourf(long_d03, lati_d03, tc_lev_d03, transform = ccrs.PlateCarree(),levels = T_levels,
                    colors = colors, extend = 'both', zorder=1)

        # Add the wind barbs
        ax2.barbs(to_np(long_d03[::,::]), to_np(lati_d03[::,::]),
              to_np(u_lev_d03[::, ::]), to_np(v_lev_d03[::, ::]),
            length=4, zorder=3, color = 'black',regrid_shape = 15)

        if plot_d04:
            ######### DOmain 4
            T_Plot = ax2.contourf(long_d04, lati_d04, tc_lev_d04, transform = ccrs.PlateCarree(),levels = T_levels,
                        colors = colors, extend = 'both', zorder = 9)

            # Add the wind barbs
            ax2.barbs(to_np(long_d04[::,::]), to_np(lati_d04[::,::]),
                  to_np(u_lev_d04[::, ::]), to_np(v_lev_d04[::, ::]),
                length=4,  color = 'black', zorder = 11, regrid_shape = 15)

            # d04 box to indicate domain 4 boundary
            ax2.add_patch(matplotlib.patches.Rectangle((left_longitude_d04, bottom_latitude_d04), right_longitude_d04-left_longitude_d04, top_latitude_d04-bottom_latitude_d04,
                     fill=None, lw=2, edgecolor='grey', zorder=50,transform = ccrs.PlateCarree(), alpha = 0.5))

            # Set title
            ax2.set_title('Domains 3-4\n700-mb Temperature\nand Winds', loc = 'right')
            
        else:
            # Set title
            ax2.set_title('Domain 3\n700-mb Temperature\nand Winds', loc = 'right')
        ax2.set_title('WRF{} Run {}\nValid {}'.format(path, run, time_of_interest_title), loc = 'left')
        
        # Add U.S. state borders and counties
        ax2.add_feature(cfeature.STATES, edgecolor='white', linewidth=0.5, zorder = 51)
        ax2.add_feature(USCOUNTIES.with_scale('500k'), edgecolor = 'white', zorder = 51)

        # Set the background color
        ax2.set_facecolor('gainsboro')

        # Add the color bar
        cax = plt.axes([0.13,0.31, 0.7, 0.02])
        cb_wa = fig.colorbar(T_Plot, cax = cax, orientation = 'horizontal', label = 'Temperature (°C)', fraction = 0.05,
                     ticks = T_levels[::4])
        cb_wa.ax.tick_params(labelsize=10)

        # Adjust spacing
        plt.subplots_adjust(hspace = -0.6, wspace = -0.7)

        # Save figure, show, and close
        plt.savefig(Fig_dir + 'WRF_{}_HRRR_{}_{}_Wind_Temp.png'.format(time_of_interest_save, 
                                                            time_of_interest_HRRR_save, level_to_interp), dpi = 300, bbox_inches = 'tight')
        #plt.show()
        plt.close()