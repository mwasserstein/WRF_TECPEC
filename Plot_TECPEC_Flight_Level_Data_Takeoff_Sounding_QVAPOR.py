# Michael Wasserstein
# Plot_TECPEC_Flight_Level_Data_Takeoff_Sounding_QVAPOR.ipynb
# 10/17/2024
# Script plots the TECPEC takeoff sounding following the flight path
# as well as the q vapor with height

# From chatGPT
# To determine the presence of a cloud using water vapor mixing ratio, you need to compare the current mixing ratio to the saturation mixing ratio at the given temperature and pressure; if the current mixing ratio is equal to or greater than the saturation mixing ratio, then conditions are suitable for cloud formation, indicating the potential for a cloud to exist at that location. 
# Key points to remember: 
# Mixing ratio:
# Represents the amount of water vapor present in a given mass of air, usually expressed in grams of water vapor per kilogram of dry air (g/kg). 
# Saturation mixing ratio:
# The maximum amount of water vapor that can exist in the air at a specific temperature and pressure before condensation occurs. 
# How to determine cloud potential: 
# 1. Measure the water vapor mixing ratio:
# This can be done using instruments like radiosondes, which measure atmospheric conditions at different altitudes, including humidity. 
# 2. Calculate the saturation mixing ratio:
# Use a saturation vapor pressure table based on the measured air temperature and pressure to find the saturation mixing ratio for that specific condition. 
# 3. Compare the values:
# If the measured water vapor mixing ratio is equal to or greater than the saturation mixing ratio, then the air is considered saturated and cloud formation is likely to occur. 

# mpcalc saturation mixing ratio

# Imports
import numpy as np
import pandas as pd
import os, sys
import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import SkewT
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair, ALL_TIMES)
import wrf
sys.path.append('/uufs/chpc.utah.edu/common/home/u1371671/')
from Useful_Python_Functions import serial_date_to_string
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import xarray as xr
import datetime
from datetime import timedelta
import netCDF4
import matplotlib.patheffects as path_effects

data_path = '/uufs/chpc.utah.edu/common/home/u1371671/steenburgh-group12/michael/TECPEC_Flight_level_data/'

# ######## User input arguments #############
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
print('Plotting data for run', run,)

# Which of your WRF runs are you interested in?
run_number = f'{run}'.zfill(2)

#Domain to plot data
domain = 3

# User input time to plot the skew T
if path in [12]:
    time_of_interest_list = pd.date_range(datetime.datetime(2019, 3, 22, 18), datetime.datetime(2019, 3, 22, 23), freq = '15min')        

# paths for data
base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{}/'.format(path)
WRF_path = base_path + 'wrf_runs/wrf_{}/'.format(run_number)
WPS_path = base_path + 'WPS/'

# paths for saving fig
parent_dir = '/uufs/chpc.utah.edu/common/home/u1371671/WRF' 
Fig_dir = parent_dir + '/Figures_{}/wrf_{}/Soundings/'.format(path, run_number)
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

# From daves email Feb 9 2024:

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

# Path of the data
data_path = '/uufs/chpc.utah.edu/common/home/u1371671/steenburgh-group12/michael/TECPEC_Flight_level_data/'

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
mr = pd.read_csv(data_path+'mr.csv',header = None).values.flatten()
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

# Columns relevant to soundings you will be plotting
columns_of_interest = ['datearray', 'latc', 'lonc', 'zmsl', 'tamb', 'tdp', 'uwind', 'vwind', 'pmb', 'licor_mr', 'mr']

# Simplify the flight data by only selecting some columns of interest (defined above) and times for the sounding
flight_data_simple = flight_data[columns_of_interest]

# Start and end times for the sounding are here
flight_data_simple = flight_data_simple[(flight_data_simple['datearray'] >= datetime.datetime(2019,3,22,20,7,0))
            & (flight_data_simple['datearray'] <= datetime.datetime(2019,3,22,20,14,40))]#.reset_index(drop = True)
                                                                                                                                                                               
# Define the window size for smoothing
window_size = 30

smooth = True # Do you want to smooth the data?

# Extract TECPEC Data
if smooth == True:
    # With smoothing
    p = np.convolve(flight_data_simple['pmb'].values ,np.ones(window_size)/window_size, mode='valid')* units('mbar')
    T = np.convolve(flight_data_simple.tamb.values ,np.ones(window_size)/window_size, mode='valid')* units('degC')
    tdp = np.convolve(flight_data_simple.tdp.values ,np.ones(window_size)/window_size, mode='valid')* units('degC')
    u = np.convolve(flight_data_simple.uwind.values,np.ones(window_size)/window_size, mode='valid') * units('m/s')
    v = np.convolve(flight_data_simple.vwind.values,np.ones(window_size)/window_size, mode='valid') * units('m/s')
    Qv = np.convolve(flight_data_simple.licor_mr.values,np.ones(window_size)/window_size, mode='valid') * units('g/kg')
    mr = np.convolve(flight_data_simple.mr.values,np.ones(window_size)/window_size, mode='valid') * units('g/kg')

else:
    # No smoothing
    p = flight_data_simple['pmb'].values * units('mbar')
    T = flight_data_simple.tamb.values * units('degC')
    tdp = flight_data_simple.tdp.values* units('degC')
    u = flight_data_simple.uwind.values * units('m/s')
    v = flight_data_simple.vwind.values * units('m/s')
    z = flight_data_simple.zmsl.values * units('m')
    Qv = flight_data_simple.licor_mr.values * units('g/kg')
    mr = flight_data_simple.mr.values * units('g/kg')


# Loop through all times in the time fo interest list
for time_of_interest in time_of_interest_list[:]:
    time_of_interest_str = datetime.datetime.strftime(time_of_interest, '%Y-%m-%d_%H:%M:%S') # convert to string for the WRF

    # load in  the wrf output data files
    data_file_d03 = WRF_path + 'run/' + 'wrfout_d0{}_{}'.format(domain, time_of_interest_str)

    # Load in wrf file using Dataset
    wrfin = Dataset(data_file_d03)

    print('Working on ', time_of_interest_str)

    # Get the wrf variables
    p1 = wrf.getvar(wrfin,"pressure",timeidx=0)
    T1 = wrf.getvar(wrfin,"tc",timeidx=0)
    Td1 = wrf.getvar(wrfin,"td",timeidx=0)
    u1 = wrf.getvar(wrfin,"ua",timeidx=0)
    v1 = wrf.getvar(wrfin,"va",timeidx=0)
    z1 = wrf.getvar(wrfin, 'z', timeidx=0)
    QVAPOR1 = wrf.getvar(wrfin,"QVAPOR",timeidx=0)

    # Times for the model run
    valid_time = pd.to_datetime(p1.Time.values)
    valid_time_str = datetime.datetime.strftime(valid_time, '%b %-d, %Y %H:%M UTC')
    valid_time_save = datetime.datetime.strftime(valid_time, '%Y%m%d%H%M')  # for saving figure
    valid_time_legend = valid_time.strftime('%H%M%S UTC') # For the legend

    # Initialization times
    init_time = wrfin.SIMULATION_START_DATE
    init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
    init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')

    # Initialize some arrays which will store the wrf data for each point
    init_p = []
    init_T = []
    init_Td = []
    init_u = []
    init_v = []
    init_Qv = []

    # Loop through all the latitudes for the flight data
    for i,value in enumerate(flight_data_simple['latc']):

        # Necessary flight level variables for the wrf data
        lat = value
        lon = flight_data_simple['lonc'].values[i]
        zmsl = flight_data_simple['zmsl'].values[i]

        # Interpolate necessary WRF variables for where the flight is
        p_at_level = wrf.interplevel(p1, z1, zmsl)
        T_at_level = wrf.interplevel(T1, z1, zmsl)
        Td_at_level = wrf.interplevel(Td1, z1, zmsl)
        u_at_level = wrf.interplevel(u1, z1, zmsl)
        v_at_level = wrf.interplevel(v1, z1, zmsl)
        Qv_at_level = wrf.interplevel(QVAPOR1, z1, zmsl)

        # Interpolating the lat and lon in the data to get an x and y for the gridpoints
        x_y = wrf.ll_to_xy(wrfin, lat, lon)

        # Extract the wrf value at the gridpoint corresponding with where the flight is
        # Note 10/16: I've double checked that [x_y[1],x_y[0]] is the correct way to extrat the values
        # at the correct gridpoint. It's strange, but it works
        p_WRF = p_at_level[x_y[1],x_y[0]].values.astype(float)#.item() #* units('mbar')
        T_WRF = T_at_level[x_y[1],x_y[0]].values.astype(float)#.item() #* units('degC')
        Td_WRF = Td_at_level[x_y[1],x_y[0]].values.astype(float)#.item()# * units('degC')
        u_WRF = u_at_level[x_y[1],x_y[0]].values.astype(float)#.item() #* units('m/s')
        v_WRF = v_at_level[x_y[1],x_y[0]].values.astype(float)#.item() #* units('m/s')
        Qv_WRF = Qv_at_level[x_y[1],x_y[0]].values.astype(float)

        # Add wrf values to initialization lists
        init_p.append(p_WRF)
        init_T.append(T_WRF)
        init_Td.append(Td_WRF)
        init_u.append(u_WRF)
        init_v.append(v_WRF)
        init_Qv.append(Qv_WRF)

    if smooth == True:
    
        # Smooth the data and add units
        init_p = gaussian_filter(np.array(init_p), sigma = 3) * units('mbar')
        init_T = gaussian_filter(np.array(init_T), sigma = 3)* units('degC')
        init_Td = gaussian_filter(np.array(init_Td), sigma = 3) * units('degC')
        init_u = gaussian_filter(np.array(init_u), sigma = 3)* units('m/s')
        init_v = gaussian_filter(np.array(init_v), sigma = 3)* units('m/s')
        init_Qv = (gaussian_filter(np.array(init_Qv), sigma = 3)*1000)*units('g/kg') # convert units to g/kg by multiply by 1000

    else:
        
        # Add units
        init_p = np.array(init_p) * units('mbar')
        init_T = np.array(init_T)* units('degC')
        init_Td = np.array(init_Td) * units('degC')
        init_u = np.array(init_u)* units('m/s')
        init_v = np.array(init_v)* units('m/s')
        init_Qv = np.array(init_Qv)*units('kg/kg')

    ################################################################################
    ################################################################################   
    ############################### Plotting ####################################

    # Create a figure with two subplots (1 row, 2 columns)
    fig = plt.figure(figsize=(14, 8), facecolor='white', edgecolor='k')


    ################################################################################
    ############################### First Subplot: skewT ###########################
    ################################################################################


    # Create a SkewT plot for the first subplot
    skew = SkewT(fig, rotation=45, subplot=(1, 2, 1))  # This creates the SkewT plot in the first position of 1x2 subplots

    # Plot TECPEC data
    skew.plot(p, T, 'c', linewidth=2.5, label='TECPEC Sounding\n200700 to 201440 UTC')
    skew.plot(p, tdp, 'c', linewidth=2.5)
    skew.plot_barbs(p[::15], u[::15], v[::15], color='c', xloc=0.97)

    # Plot WRF data
    skew.plot(init_p, init_T, 'm', linewidth=2.5, label=f'WRF Simulation\n{valid_time_legend}')
    skew.plot(init_p, init_Td, 'm', linewidth=2.5)
    skew.plot_barbs(init_p[::15], init_u[::15], init_v[::15], color='m', xloc=1.05)

    # Add relevant special lines
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats(alpha=0.2)
    skew.plot_mixing_lines(alpha=0.2)

    # Set limits
    skew.ax.set_xlim(-10, 20)
    skew.ax.set_ylim(900, 500)

    # Add labels
    skew.ax.set_xlabel('Temperature ($^\circ$C)', fontsize=18)
    skew.ax.set_ylabel('Pressure (hPa)', fontsize=18)

    # Set ticks
    skew.ax.set_xticks(np.arange(-10, 21, 5))
    skew.ax.tick_params(axis='x', labelsize=16)
    skew.ax.tick_params(axis='y', labelsize=16)

    # Add legend
    skew.ax.legend()

    # Add some titles
    skew.ax.set_title(f'WRF{path} Run {run_number}\nInit: {init_time_str}', loc='left')
    skew.ax.set_title(f'Domain {domain} Sounding\nValid: {valid_time_str}', loc='right')


    ################################################################################
    ############################### Second Subplot: Qvapor #########################
    ################################################################################


    # Create the second subplot for Qvapor
    # MW 10/17/24: What I'm doing here is "fudging" the axis so it aligns well with the axis for 
    # The skewT. I couldn't figure out better code to do this.
    # THis is probably bad python, but it works
    ax2 = fig.add_axes([0.60, 0.245, 0.35, 0.515])  

    # Plot TECPEC and WRF data in second subplot (ax2)
    ax2.plot(Qv, p, 'c', linewidth=2.5, label='TECPEC Licor probe\n200700 to 201440 UTC')
    ax2.plot(init_Qv, init_p, 'm', linewidth=2.5, label=f'WRF Simulation\n{valid_time_legend}')
    ax2.plot(mr, p, 'c', linestyle = '--', linewidth = 2.5, label = 'TECPEC\n200700 to 201440 UTC')

    # Set y-axis limits and scale
    ax2.set_ylim(900, 500)
    ax2.set_yscale('log')

    # Set x limits
    ax2.set_xlim(1.2,5)

    # Add labels
    ax2.set_xlabel('Water Vapor Mixing Ratio (g/kg)', fontsize=18)
    ax2.set_ylabel('Pressure (hPa)', fontsize=18)

    # Add ticks
    ax2.set_yticks([900, 800, 700, 600, 500])
    ax2.set_yticklabels([900, 800, 700, 600, 500], fontsize=14)
    ax2.tick_params(axis='x', labelsize=14)

    # Add titles
    ax2.set_title(f'WRF{path} Run {run_number}\nInit {init_time_str}', loc='left')
    ax2.set_title(f'Domain {domain}\nValid {valid_time_str}', loc='right')

    # Add a grid
    ax2.grid()

    # Add legend
    ax2.legend()

    # Adjust layout
    plt.tight_layout()

    # Save figure and show
    plt.savefig(Fig_dir + 'WRF_SkewT_Qvapor_{}.png'.format(valid_time_save), dpi = 200, bbox_inches='tight')
    #plt.show()
    plt.close()