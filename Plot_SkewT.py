# Michael Wasserstein
# Plot_SkewT_time_offset.py
# 10/9/2024
# Script takes in WRF outputs and plots a skewT , offset by a specified amount of time

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts/Plot_SkewT.py -r 2 -s SLC -p 2 -i 2019032212 -e 2019032400
# -r represents the run number you want to plot
# -s represents the station/location of interest
# -p represents the path of the wrf run
# -i represents the start time for the model run YYYYmmddHH
# -e represents the end time for the model run YYYYmmddHH

import wrf
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#import proplot as pplt
import glob
import xarray
import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units
import datetime
import pandas as pd
from siphon.simplewebservice.wyoming import WyomingUpperAir
import os

######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-s", "--station", help="Station of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")
parser.add_argument("-i", "--init", help="Start time for skew t's")
parser.add_argument("-e", "--end", help="end time for skew t's")


args = parser.parse_args()

# Get user inputs
run = int(args.run)
station = str(args.station)
path = int(args.path)
start = str(args.init)
end = str(args.end)
print('Plotting SkewT for run', run, 'at', station)

run_number = '{}'.format(run).zfill(2)

# Time offset - do you want it or not?
time_offset = False
time_offset_minutes = 45
time_offset_hours = 0


# paths for data
if path ==1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
else:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/'.format(run_number)
WPS_path = base_path + 'WPS/'

# paths for saving fig
parent_dir = '/uufs/chpc.utah.edu/common/home/u1371671/WRF'
if time_offset == True:
    Fig_dir = parent_dir + '/Figures_{}/wrf_{}/Soundings_time_offset/'.format(path,run_number)
else:
    Fig_dir = parent_dir + '/Figures_{}/wrf_{}/Soundings/'.format(path,run_number)
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

# define the start and end time for the skew t plots
start_time = datetime.datetime.strptime(start, '%Y%m%d%H')
end_time   = datetime.datetime.strptime(end, '%Y%m%d%H')

time_of_interest_list = pd.date_range(start = start_time, end = end_time, freq = '12H')

lat_lon_dict = {'SLC' : [40.77, -111.95]}
lat_lon = lat_lon_dict[station] # lat and lon of slc

for time_of_interest in time_of_interest_list:
    time_of_interest_str = datetime.datetime.strftime(time_of_interest, '%Y-%m-%d_%H:%M:%S')
    if time_offset:
        time_of_interest_offset = time_of_interest + datetime.timedelta(hours = time_offset_hours, minutes = time_offset_minutes)
        time_of_interest_offset_str = datetime.datetime.strftime(time_of_interest_offset, '%Y-%m-%d_%H:%M:%S')
        file = WRF_path + 'run/' + 'wrfout_d03_' + time_of_interest_offset_str

    else:
        file = WRF_path + 'run/' + 'wrfout_d03_' + time_of_interest_str

    wrfin = Dataset(file) # Extract the actual file

    x_y = wrf.ll_to_xy(wrfin, lat_lon[0], lat_lon[1]) # some interpolations for the xy of the data

    # Load in needed wrf variables for the script
    p1 = wrf.getvar(wrfin,"pressure",timeidx=0)
    T1 = wrf.getvar(wrfin,"tc",timeidx=0)
    Td1 = wrf.getvar(wrfin,"td",timeidx=0)
    u1 = wrf.getvar(wrfin,"ua",timeidx=0)
    v1 = wrf.getvar(wrfin,"va",timeidx=0)

    # Get the wrf variables for the correct location and add units
    p = p1[:,x_y[1],x_y[0]] * units.hPa
    T = T1[:,x_y[1],x_y[0]] * units.degC
    Td = Td1[:,x_y[1],x_y[0]] * units.degC
    u = u1[:,x_y[1],x_y[0]] * units('m/s')
    v = v1[:,x_y[1],x_y[0]] * units('m/s')

    # Times for the model run
    valid_time = pd.to_datetime(p1.Time.values)
    valid_time_str = datetime.datetime.strftime(valid_time, '%b %-d, %Y %H:%M UTC')
    valid_time_save = datetime.datetime.strftime(valid_time, '%Y%m%d%H%M')  # for saving figure

    # Stuff for the initialization time
    init_time = wrfin.SIMULATION_START_DATE
    init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
    init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')
    
    # Get the real observed sounding data
    try:
        # Load in the dataset
        df = WyomingUpperAir.request_data(time_of_interest, station)
        
        # Get the times for the observations
        time_obs = pd.to_datetime(df.time.values[0])
        time_obs_str = datetime.datetime.strftime(time_obs, '%b %-d, %Y %H:%M UTC')
        
        # prepare the data in the units that you want
        p_obs = df['pressure'].values * units(df.units['pressure'])
        T_obs = df['temperature'].values * units(df.units['temperature'])
        Td_obs = df['dewpoint'].values * units(df.units['dewpoint'])
        u_obs = (df['u_wind'].values * units(df.units['u_wind'])).to(units.meter_per_second)
        v_obs = (df['v_wind'].values * units(df.units['v_wind'])).to(units.meter_per_second)
        
        
        #################################### Plot Figure ###################################
        fig = plt.figure(figsize=(8, 12), facecolor = 'white', edgecolor = 'k')

        skew = SkewT(fig, rotation=45)

        # Plot the data using normal plotting functions, in this case using
        # log scaling in Y, as dictated by the typical meteorological plot
        skew.plot(p, T, 'm', linewidth = 2.5, label = 'WRF Simulation ' + valid_time_str)
        skew.plot(p, Td, 'm', linewidth = 2.5)
        skew.plot_barbs(p, u, v, color = 'm', xloc = 1.05)

        skew.plot(p_obs, T_obs, 'c', linestyle = '-', linewidth = 2.5, label = 'KSLC Sounding ' + time_obs_str)
        skew.plot(p_obs, Td_obs, 'c', linestyle = '-', linewidth = 2.5)
        skew.plot_barbs(p_obs[::2], u_obs[::2], v_obs[::2], color = 'c', xloc = 0.95)

        # Add the relevant special lines
        skew.plot_dry_adiabats()
        skew.plot_moist_adiabats()
        skew.plot_mixing_lines()
        skew.ax.set_xlim(-60, 40)
        skew.ax.set_xlabel('Temperature ($^\circ$C)')
        skew.ax.set_ylabel('Pressure (hPa)')

        # Add some titles
        plt.title('WRF {} Run {}\nInit: {}'.format(path, run, init_time_str), loc='left')
        plt.title('Valid: {}'.format(valid_time_str), loc='right')
        
        # TIcks and labels
        plt.xticks(fontsize = 16)
        plt.xlabel('Temperature (°C)', fontsize = 18)

        plt.yticks(fontsize = 16)
        plt.ylabel('Pressure (hPa)', fontsize = 18)

        # Moist adiabats and mixing lines
        skew.plot_moist_adiabats(alpha = 0.2)
        skew.plot_mixing_lines(alpha = 0.2)

        skew.ax.set_xlim(-30, 40)

        plt.legend(loc = 'upper left', prop={'size': 14})


        plt.savefig(Fig_dir + '/WRF_SND_{}_{}.png'.format(station, valid_time_save), bbox_inches='tight')
        plt.close()
        #plt.show()


    except:
        print('Real sounding data for ' + time_of_interest_str + ' is not available')
