# Michael Wasserstein
# Plot_Time_Height.py
# 11/25/2024
# Script takes in WRF outputs and plots a time height for the site of interest

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_Time_Height.py -r 2 -p 2 
# -r represents the run number you want to plot
# -p represents the path of interest

# Imports
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

sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from map_script import *

######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")

args = parser.parse_args()

###########################
##### stuff for WRF  ######
###########################
# Get user inputs
run = str(args.run)
path = int(args.path)

# leading zeros
run_number = '{}'.format(run).zfill(2)

# paths for data
if path == 1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
else:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number)
WPS_path = base_path + 'WPS/'

# paths for saving fig
parent_dir = '/uufs/chpc.utah.edu/common/home/u1371671/WRF'
Fig_dir = parent_dir + '/Figures_{}/wrf_{}/Time_Heights/'.format(path, run_number)

# If it doesn't already exist, make directory
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

########### End of user inputs ############
# stuff for latitude longitude
station_dicts = {
                'CLN' : [40.5763,-111.6383], # Alta collins
                'SLC' : [40.77,-111.95],   # SLC airport
                'SVD' : [40.52, -111.87], # South part of valley near Draper
                'WVH' : [40.52, -112.07], # West part of valley near Herriman
                'KSY' : [40.45, -112.21], # Kelsey peak in oquirrhs
                'DST' : [40.459, -112.627], # Deseret peak in stansburys
                'TDT' : [40.53, - 112.285] # Tooele downtown    
                }   


# Load in the data files and sort
data_files_d03 = glob.glob(WRF_path + '*wrfout_d03*') # for the innermost domain
data_files_d03.sort()

data_files_d03 = data_files_d03[::4] # for now, only get every 4th filen (this speeds up some of the wrf python code)

# Load in all files and make into a list
wrflist = [Dataset(file) for file in data_files_d03]

# Use the file names to get the start and end time for your plotting output
start_time = datetime.datetime.strptime(data_files_d03[0].split('/')[-1][-19:], '%Y-%m-%d_%H:%M:%S')
end_time = datetime.datetime.strptime(data_files_d03[-1].split('/')[-1][-19:], '%Y-%m-%d_%H:%M:%S')

valid_time_str = datetime.datetime.strftime(end_time, '%b %-d, %Y %H:%M UTC') # String for title for end time

# Intitialization time
init_time = wrflist[0].SIMULATION_START_DATE
init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')

# Stuff for boundaries
geobounds = wrf.geo_bounds(wrfin=wrflist[0])
bottom_latitude = geobounds.bottom_left.lat
left_longitude = geobounds.bottom_left.lon
top_latitude = geobounds.top_right.lat
right_longitude = geobounds.top_right.lon

# Extract data for all times of intereset
tc = getvar(wrflist, "tc", timeidx=ALL_TIMES, method="cat")
p = getvar(wrflist,"pressure", timeidx=ALL_TIMES, method="cat")
rh = getvar(wrflist,"rh", timeidx=ALL_TIMES, method="cat")
ua = getvar(wrflist,"ua", timeidx=ALL_TIMES, method="cat")
va = getvar(wrflist,"va", timeidx=ALL_TIMES, method="cat")
wa = getvar(wrflist,"wa", timeidx=ALL_TIMES, method="cat")
thte = getvar(wrflist,"eth", timeidx=ALL_TIMES, method="cat")
terr_d03 = getvar(wrflist,"HGT", timeidx=0, method="cat")
LANDMASK_d03 = getvar(wrflist, 'LANDMASK', timeidx=-1) # only needed for locator map

# Determine duration of run / time period to plot time height in nanoseconds
timedelta_run = ua.Time[-1].values - ua.Time[0].values

# Convert to hours
hours = timedelta_run / np.timedelta64(1, 'h')

# Setting for plot
rh_fill_colors = ['#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30']
fsize = 10
sampint = 1 # interval at which you want to sample the data
hbarbint = 1 # horizontal barb interval
vbarbint = 2 # vertical barb inteval
timeint = -2 # TIme interval (use the minus sign because plot is right to left)

for station in station_dicts.keys():
    
    # extract latitude and longitude from dictionary
    lat_lon = station_dicts[station] 

    # some interpolations for the xy of the data based on the station of use
    x_y = wrf.ll_to_xy(wrflist, lat_lon[0], lat_lon[1]) 

    # extract the values at the location of the station
    tc_plot = tc[:,:,x_y[1], x_y[0]] # deg C
    p_plot = p[:,:,x_y[1], x_y[0]]   # hPa
    rh_plot = rh[:,:,x_y[1], x_y[0]] # %
    ua_plot = ua[:,:,x_y[1], x_y[0]] # m / s
    va_plot = va[:,:,x_y[1], x_y[0]] # m / s
    wa_plot = wa[:,:,x_y[1], x_y[0]] # m / s
    thte_plot = thte[:,:,x_y[1], x_y[0]] # Kelvin

    # Determine how many hours are in plot
    numtimes = len(thte_plot.Time)

    # Determine number of levels
    numlevels = len(thte_plot.bottom_top)

    # Set timeindex (should have shape (numtimes, numlevels))
    timeindex = np.full((numtimes,numlevels),np.nan) 
    for i in range (0,numtimes):
        timeindex[i] = float(i)

    # Intiate and set figure size
    fig = plt.figure(figsize=(9,6), facecolor = 'white', edgecolor = 'k')

    ################# ax main plot ######################  
    ax = fig.add_subplot(1,1,1)

    # Invert x axis to time increases to left per convention
    ax.invert_xaxis()

    # Plot rh as colorfil
    plot_rh = ax.contourf(timeindex[:numtimes+1:sampint,:],p_plot[:numtimes+1:sampint,:],rh_plot[:numtimes+1:sampint,:],
                           np.arange(0,101,10),colors = rh_fill_colors, alpha=0.5, extend='max', antialiased=True)

    # Plot vertical velocity contours
    # Positive contours red
    contours = [.05,.25,.5] # use those because it effectively converts to cm/s
    ax.contour(timeindex[:numtimes+1:sampint,:],p_plot[:numtimes+1:sampint,:],wa_plot[:numtimes+1:sampint,:],contours,colors=['#fb6a4a','#de2d26','#a50f15'],linewidths=.5)

    # Negative contours blue
    contours = [-.5, -.25, -.05] # use those because it effectively converts to cm/s
    ax.contour(timeindex[:numtimes+1:sampint,:],p_plot[:numtimes+1:sampint,:],wa_plot[:numtimes+1:sampint,:],contours,colors=['#08519c','#3182bd','#6baed6'],linewidths=.5)

    # Plot thetae as contour plot
    thetaecont = ax.contour(timeindex[:numtimes+1:sampint,:],p_plot[:numtimes+1:sampint,:],thte_plot[:numtimes+1:sampint,:],np.arange(200,500,2),colors='k',linewidths=1)
    ax.clabel(thetaecont, fmt='%.0f', inline=1, fontsize=fsize)

    # Plot 0˚C isotherm
    ax.contour(timeindex[:numtimes+1:sampint,:],p_plot[:numtimes+1:sampint,:],tc_plot[:numtimes+1:sampint,:],np.arange(0,500,500),colors='b',linewidths=2)

    # Plot barbs
    ax.barbs(timeindex[:numtimes+1:hbarbint,::vbarbint],p_plot[:numtimes+1:hbarbint,::vbarbint],
              ua_plot[:numtimes+1:hbarbint,::vbarbint],va_plot[:numtimes+1:hbarbint,::vbarbint], length=5, linewidth=0.25) 

    # x- ticks
    ticklabels=pd.to_datetime(ua.Time.values)[::timeint].strftime('%b %d\n%H:%M') 
    ax.set_xticks(timeindex[:,0][::timeint])   
    ax.set_xticklabels(ticklabels) 

    # y-axis ticks 
    pbot = np.nanmax(p_plot)
    yticks=np.arange(300,pbot+1,50) 
    ax.set_yticks(yticks) 
    ax.set_ylim(pbot,300)    # Set y limits (will plot from highest pressure level in forecast period to 300 mb

    # y-axis labels
    ax.set_ylabel('Pressure (hPa)')

    # Titles
    ax.set_title(f'WRF{run} run {path}\n' +station+" Time-Height Section\n"+"Init "+init_time_str,
                  loc='left', fontsize=fsize)
    ax.set_title(f'{hours}-hr forecast period\n' +'through ' + valid_time_str , loc='right', fontsize=fsize)

    # Legend for thetae
    thetaetext = fig.text(0.072,0.84,'THETA-E',color='k',horizontalalignment='center',fontsize=fsize,transform=ax.transAxes)
    thetaetext.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
                       path_effects.Normal()])

    # Legend for 0˚C isotherm
    tmptext = fig.text(0.072,0.805,'TEMP=0˚C',color='b',horizontalalignment='center',fontsize=fsize,transform=ax.transAxes)
    tmptext.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
                       path_effects.Normal()])

    # Legend for upward vertical velocity
    tmptext = fig.text(0.072,0.77,'UVV',color='#fb6a4a',horizontalalignment='center',fontsize=fsize,transform=ax.transAxes)
    tmptext.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
                       path_effects.Normal()])

    # Legend for downward vertical velocity
    tmptext = fig.text(0.072,0.735,'DVV',color='#6baed6',horizontalalignment='center',fontsize=fsize,transform=ax.transAxes)
    tmptext.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
                       path_effects.Normal()])

    # Add rh colorbar
    cax = plt.axes([0.92,0.12, 0.02, 0.4])
    cb = plt.colorbar(plot_rh, cax=cax, shrink=.6, pad=0.02, drawedges=True, orientation='vertical')
    cb.ax.tick_params(labelsize=fsize)
    cb.set_label("Relative Humidity (%)", fontsize=fsize)

    ################# Ax2 Locator Map ######################  
    # Set projection and extent of plot
    ax2 = fig.add_axes( [0.87,0.57,0.3,0.3], projection=ShadedReliefESRI().crs) # Need to add a new axis
    ax2.set_extent([left_longitude, right_longitude, bottom_latitude, top_latitude])

    # Add terrain to locator map
    ax2.contourf(terr_d03.XLONG.values, terr_d03.XLAT.values, terr_d03.values, cmap = 'terrain', levels = np.arange(0,3500,250),
                 transform = ccrs.PlateCarree(), zorder = 1)

    # Add lake to locator map
    ax2.contourf(LANDMASK_d03.XLONG.values, LANDMASK_d03.XLAT.values, LANDMASK_d03.values, colors = 'midnightblue', levels = np.arange(0,0.6,0.5),
                 transform = ccrs.PlateCarree(), zorder = 1)

    # Add locator marker
    ax2.scatter(lat_lon[1], lat_lon[0], color = 'red', s = 150, edgecolor = 'black', marker = '*', transform=ccrs.PlateCarree(), zorder = 100)

    # Save show and close
    plt.savefig(Fig_dir + 'Time_Height_{}.png'.format(station), dpi = 300, bbox_inches = 'tight')
    #plt.show()
    plt.close()