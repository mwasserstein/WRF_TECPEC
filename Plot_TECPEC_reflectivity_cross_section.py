# Michael Wasserstein
# Plot_TECPEC_reflectivity_cross_section.py
# 10/20/2024
# Script takes in WRF outputs that had been converted
# to w-band using CRSIM and plots a cross section of simulated
# reflectivity along flight legs for domains 3 and 4 or just 3

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts_To_Use_Now/Plot_TECPEC_reflectivity_cross_section.py -r 2 -p 2
# -r represents the run number you want to plot
# -p is the wrf path (wrf1 or wrf2)

# Imports
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import from_levels_and_colors
from matplotlib.image import imread
from matplotlib.colors import LinearSegmentedColormap
from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair)
import wrf
import glob
import os, sys
import pandas as pd
import datetime
import pyart
import xarray as xr

######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")

args = parser.parse_args()

# Get user inputs
run = str(args.run)
path = int(args.path)
print('Plotting data for run', run,)

run_number = '{}'.format(run).zfill(2) # Get the run with leading zeros

plot_d04 = False # do you want to plot data for d04 (true) or d03 only (False)

# Which CRSIM variable do you want to plot?
variable_to_plot = 'Zhh' # Zhh is reflectivity at hh polarization

# TImes at which to plot the data
time_of_interest_list = pd.date_range(datetime.datetime(2019,3,22,20,15,0), datetime.datetime(2019,3,23,0,0,0), freq = '15min')
##################################### End User inputs ###############################

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

# CRSIM Data path
CRSIM_out_path = base + 'steenburgh-group12/michael/CRSIM_4/wrf{}/wrf_{}/out_files_30m/'.format(path, run_number)

# paths for saving fig
if plot_d04:
    Fig_dir = home + 'WRF/Figures_{}/wrf_{}/TECPEC_Cross_Sections_Reflectivity_QVAPOR_CRSIM_d03_d04/'.format(path,run_number)
else:
    Fig_dir = home + 'WRF/Figures_{}/wrf_{}/TECPEC_Cross_Sections_Reflectivity_QVAPOR_CRSIM_d03/'.format(path,run_number)

# make fig dir if it does not already exist
if os.path.exists(Fig_dir) == False:
    os.mkdir(Fig_dir)

# Load in TECPEC flight leg data
flight_df = pd.read_csv(base + 'steenburgh-group12/michael/TECPEC_Flight_level_data/TECPEC_leg_data.csv')

##### Plot settings (Doing here so its outside the for loop)
dbz_levels = np.arange(-30, 31., 1.0)
# Cmap for plotting
dbz_map = plt.get_cmap('jet').copy()
dbz_map.set_under('none')

for i in range(len(flight_df)):   # Loop through all flight legs
    # FLight info
    flight_leg = flight_df['Leg'][i]
    track = flight_df['Path'][i] # E or W
    print('Working on flight leg', flight_leg)

    # Extract the start and end points of the cross section based on what the track is
    if track[0] == 'E':
        # cross section start and end points
        cross_start = CoordPair(lat=flight_df['End_lat'][i], lon=flight_df['End_lon'][i])
        cross_end = CoordPair(lat=flight_df['Start_lat'][i], lon=flight_df['Start_lon'][i])
    elif track[0] == 'W':
        cross_start = CoordPair(lat=flight_df['Start_lat'][i], lon=flight_df['Start_lon'][i])
        cross_end = CoordPair(lat=flight_df['End_lat'][i], lon=flight_df['End_lon'][i])
    else:
        print('There is an issue with the track. Check your work.')
        break
    
    for time_of_interest in time_of_interest_list[:]:
        time_of_interest_str = time_of_interest.strftime('%Y-%m-%d_%H:%M:%S')

        ############ CR-SIM information #############
        # DOmain 3 CRSIM
        domain = '3'.zfill(2)
        crsim_file_name = f'crsimout_d{domain}_{time_of_interest_str}'
        crsim_file_d03 = CRSIM_out_path + crsim_file_name
        ds_CRSIM_d03 = xr.open_dataset(crsim_file_d03)
        #DOmain 3 WRF
        WRF_file = WRF_path + f'wrfout_d{domain}_{time_of_interest_str}'
        wrf_file_d03 = Dataset(WRF_file)

        if plot_d04:
            # DOmain 4 CRSIM
            domain = '4'.zfill(2)
            crsim_file_name = f'crsimout_d{domain}_{time_of_interest_str}'
            crsim_file_d04 = CRSIM_out_path + crsim_file_name
            ds_CRSIM_d04 = xr.open_dataset(crsim_file_d04)
            #DOmain 4 WRF
            WRF_file = WRF_path + f'wrfout_d{domain}_{time_of_interest_str}'
            wrf_file_d04 = Dataset(WRF_file)

            # Extract bounds for d04
            geobounds = wrf.geo_bounds(wrfin=wrf_file_d04)
            min_lat = geobounds.bottom_left.lat
            left_lon = geobounds.bottom_left.lon
            max_lat = geobounds.top_right.lat
            right_lon = geobounds.top_right.lon

        #################################################################################
        #############################   Domain 3 #######################################
        #################################################################################

        # Get the CRSIM variables for domain 3
        dbz_d03 = ds_CRSIM_d03[variable_to_plot]
        Z_d03 = 10**(dbz_d03/10.) # Use linear Z for interpolation # gets rid of the log
        ht_d03 = ds_CRSIM_d03['height']

        # Get the WRF variables for domain 3
        ter_d03 = getvar(wrf_file_d03, "ter", timeidx=-1)
        QVAPOR_d03 = getvar(wrf_file_d03, 'QVAPOR', timeidx = -1)
        ht_d03 = getvar(wrf_file_d03, "z", timeidx=-1)

        # Compute the vertical cross-section interpolation.  Also, include the
        # lat/lon points along the cross-section in the metadata by setting latlon
        # to True.
        z_cross_d03 = vertcross(Z_d03, ht_d03, wrfin=wrf_file_d03,
                            start_point=cross_start,
                            end_point=cross_end,
                            latlon=True, meta=True)

        # Convert back to dBz after interpolation
        dbz_cross_d03 = 10.0 * np.log10(z_cross_d03)

        # Add back the attributes that xarray dropped from the operations above
        dbz_cross_d03.attrs.update(z_cross_d03.attrs)
        dbz_cross_d03.attrs["description"] = "radar reflectivity cross section"
        dbz_cross_d03.attrs["units"] = "dBZ"

        # Make a copy of the z cross data. Let's use regular numpy arrays for this.
        dbz_cross_filled_d03 = np.ma.copy(to_np(dbz_cross_d03))

        try:
        # For each cross section column, find the first index with non-missing
        # values and copy these to the missing elements below.
            for i in range(dbz_cross_filled_d03.shape[-1]):
                column_vals_d03 = dbz_cross_filled_d03[:,i]
                # Let's find the lowest index that isn't filled. The nonzero function
                # finds all unmasked values greater than 0. Since 0 is a valid value
                # for dBZ, let's change that threshold to be -200 dBZ instead.
                first_idx_d03 = int(np.transpose((column_vals_d03 > -200).nonzero())[0])
                dbz_cross_filled_d03[0:first_idx_d03, i] = dbz_cross_filled[first_idx_d03, i]

        except:
            pass

        # Get the terrain heights along the cross section line
        ter_line_d03 = interpline(ter_d03, wrfin=wrf_file_d03, start_point=cross_start,
                              end_point=cross_end)

        # QVAPOR for d03
        QVAPOR_cross_d03 = vertcross(QVAPOR_d03, ht_d03, wrfin=wrf_file_d03, start_point=cross_start,
                   end_point=cross_end, latlon=True, meta=True)

        # Get the initialization times
        init_time = wrf_file_d03.SIMULATION_START_DATE
        init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
        init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')

        # Valid times
        valid_time = pd.to_datetime(ter_d03.Time.values)
        valid_time_str = datetime.datetime.strftime(valid_time, '%b %-d, %Y %H:%M UTC')

        # TImes for saving path
        YYYY = valid_time.year
        mm = valid_time.month
        DD = valid_time.day
        HH = valid_time.hour
        MM = valid_time.minute

        # Path for saving figure
        save_path = Fig_dir + 'Radar_Cross_Sect_{}_{}{:02}{:02}{:02}{:02}.png'.format(flight_leg, YYYY,mm, DD, HH, MM)

        #################################################################################
        #############################   Domain 4 #######################################
        #################################################################################
        if plot_d04:

            # Filter to include only coordinate pairs with longitude greater than left_lon
            # This line of code ensures that you wont be taking a cross section with stuff from outside of the domain
            filtered_coord_pairs = [cp for cp in w_cross_d03.xy_loc.values if cp.lon > left_lon]

            # Find the coordinate pair with the longitude closest to left_lon among the filtered pairs
            # Essentially here, you're finding the coordinate pair from the d03 cross section that is closest to the boundary of domain 4
            closest_coord_pair_min = min(filtered_coord_pairs, key=lambda cp: abs(cp.lon - left_lon))
            closest_coord_pair_max = min(w_cross_d03.xy_loc.values, key=lambda cp: abs(cp.lon - right_lon)) # Same for right lon


            # Not totally sure why steps 3 and 4 need to be done, but it was givent to me by CHATGPT. 
            # Basically, the closest one to the boundary didnt work, so you need to use one more inwards.
            # Step 3: Find the index of the closest coordinate pair in the original list
            original_list = list(w_cross_d03.xy_loc.values)  # Convert to list if not already a list
            index_of_closest = original_list.index(closest_coord_pair_min)

            # Step 4: Get the next coordinate pair in the list
            if index_of_closest < len(original_list) - 1:  # Ensure it's not the last element
                next_coord_pair = original_list[index_of_closest + 1]
            else:
                next_coord_pair = None  # Handle the case where there is no next pair

            # Extract the x values for the d04 plot that are at the start and end of the cross section
            x_min_from_d04 = closest_coord_pair_min.x
            x_max_from_d04 = closest_coord_pair_max.x


            # Coordinates for the start and end of the cross section for d04. You are doing this based on the part of d03 that is within d04
            cross_start_for_d04 = next_coord_pair
            cross_end_for_d04 = CoordPair(lon=closest_coord_pair_max.lon, lat=closest_coord_pair_max.lat)


            # Get the CRSIM variables for domain 4
            dbz_d04 = ds_CRSIM_d04[variable_to_plot]
            Z_d04 = 10**(dbz_d04/10.) # Use linear Z for interpolation # gets rid of the log
            ht_d04 = ds_CRSIM_d04['height']

            # Get the WRF variables for domain 4
            ter_d04 = getvar(wrf_file_d04, "ter", timeidx=-1)
            QVAPOR_d04 = getvar(wrf_file_d04, 'QVAPOR', timeidx = -1)
            ht_d04 = getvar(wrf_file_d04, "z", timeidx=-1)

            # Compute the vertical cross-section interpolation.  Also, include the
            # lat/lon points along the cross-section in the metadata by setting latlon
            # to True.
            z_cross_d04 = vertcross(Z_d04, ht_d04, wrfin=wrf_file_d04,
                                start_point=cross_start,
                                end_point=cross_end,
                                latlon=True, meta=True)

            # Convert back to dBz after interpolation
            dbz_cross_d04 = 10.0 * np.log10(z_cross_d04)

            # Add back the attributes that xarray dropped from the operations above
            dbz_cross_d04.attrs.update(z_cross_d04.attrs)
            dbz_cross_d04.attrs["description"] = "radar reflectivity cross section"
            dbz_cross_d04.attrs["units"] = "dBZ"

            # Make a copy of the z cross data. Let's use regular numpy arrays for this.
            dbz_cross_filled_d04 = np.ma.copy(to_np(dbz_cross_d04))

            try:
            # For each cross section column, find the first index with non-missing
            # values and copy these to the missing elements below.
                for i in range(dbz_cross_filled_d04.shape[-1]):
                    column_vals_d04 = dbz_cross_filled_d04[:,i]
                    # Let's find the lowest index that isn't filled. The nonzero function
                    # finds all unmasked values greater than 0. Since 0 is a valid value
                    # for dBZ, let's change that threshold to be -200 dBZ instead.
                    first_idx_d04 = int(np.transpose((column_vals_d04 > -200).nonzero())[0])
                    dbz_cross_filled_d04[0:first_idx_d04, i] = dbz_cross_filled[first_idx_d04, i]

            except:
                pass

            # Get the terrain heights along the cross section line
            ter_line_d04 = interpline(ter_d04, wrfin=wrf_file_d04, start_point=cross_start,
                                  end_point=cross_end)

            # QVAPOR for d04
            QVAPOR_cross_d04 = vertcross(QVAPOR_d04, ht_d04, wrfin=wrf_file_d04, start_point=cross_start,
                       end_point=cross_end, latlon=True, meta=True)

        #######################################################################################
        #################################### Plotting #########################################
        #######################################################################################
        # Create the figure
        fig, (ax ) = plt.subplots(1,1, figsize=(8,6), facecolor = 'white', edgecolor = 'k')

        # Setting the background color of the plot
        # using set_facecolor() method
        ax.set_facecolor("gainsboro")

        ####################################### d03 plot ##########################################
        # Make the cross section plot for w
        xs_d03 = np.arange(0, z_cross_d03.shape[-1], 1)
        ys_d03 = to_np(z_cross_d03.coords["vertical"])
        dbz_contours_d03 = ax.contourf(xs_d03,
                                            ys_d03,
                                            to_np(dbz_cross_filled_d03),
                                            levels=dbz_levels,
                                            cmap=dbz_map,
                                            extend="neither",zorder = 1)

        # QVAPOR Contour
        QVAPOR_contour = ax.contour(xs_d03,ys_d03,
                            to_np(QVAPOR_cross_d03)*1000,levels=np.arange(0,6,0.4),colors='black',
                            linestyles='dashed',zorder = 2)
        plt.clabel(QVAPOR_contour, inline=1, fontsize=12,zorder = 2.1)

        # Fill in the mountain area
        ht_fill = ax.fill_between(xs_d03, 0, to_np(ter_line_d03),
                                        facecolor="saddlebrown")

        # Do all the ticks and limits before you do the d04 plot
        #Set the x-ticks to use latitude and longitude labels
        coord_pairs = to_np(QVAPOR_cross_d03.coords["xy_loc"])
        x_ticks = np.arange(coord_pairs.shape[0])
        x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]

        # Set the desired number of x ticks below
        num_ticks = 5
        thin = int((len(x_ticks) / num_ticks) + .5)
        ax.set_xticks(x_ticks[::thin])
        ax.set_xticklabels(x_labels[::thin], rotation=45, fontsize=8)

        # Set the x-axis and  y-axis labels
        ax.set_xlabel("Latitude, Longitude", fontsize=12)
        ax.set_ylabel("Height (m)", fontsize=12)

        # set y limit
        ax.set_ylim(0,10000)

        ####################################### d04 plot ##########################################
        if plot_d04: # Plot data for d04 if you want it
            # Here you are extracting the x coordinates in domain 3 plot of the start and end of the cross section for domain 4
            start_x_from_d04 = wrf.ll_to_xy(wrf_file_d03, latitude=cross_start_for_d04.lat, longitude=cross_start_for_d04.lon).sel(x_y = 'x').values
            end_x_from_d04 = wrf.ll_to_xy(wrf_file_d03, latitude=cross_end_for_d04.lat, longitude=cross_end_for_d04.lon).sel(x_y = 'x').values


            # Make an array of x coordinates and y coordinates for plotting d03
            # Not sure why you need to subtract the first x from d03 cross, but it works
            if flight_leg in [1,2,3,12,13,14]:
                xs_d04 = np.linspace(start_x_from_d04-w_cross_d03.xy_loc.values[0].x, end_x_from_d04-w_cross_d03.xy_loc.values[0].x,w_cross_d04.shape[1])
            # For the legs where the end of the leg is at the same spot as d03, just use the last point from d03
            elif flight_leg in [4,5,6,7,8,9,10,11]:
                xs_d04 = np.linspace(start_x_from_d04-w_cross_d03.xy_loc.values[0].x, xs_d03[-1],w_cross_d04.shape[1])
            ys_d04 = to_np(w_cross_d04.coords["vertical"])

            # Plot w for d04
            dbz_contours_d04 = ax.contourf(xs_d04,
                                                ys_d04,
                                                to_np(dbz_cross_d04),
                                                levels=dbz_levels,
                                                cmap=my_cmap,
                                                extend="neither", zorder = 3)

            # Plot QVAPOR
            QVAPOR_contours = ax.contour(xs_d04,ys_d04,
                       to_np(QVAPOR_cross_d04)*1000,levels=np.arange(0,6,0.4),colors='black',
                       linestyles='dashed', zorder = 4)
            plt.clabel(QVAPOR_contours, inline=1, fontsize=12, zorder = 4.1)


            # Vertical Lines denoting the d04 part of the cross section
            if flight_leg in [1,2,3,12,13,14]:
                ax.vlines([start_x_from_d04-w_cross_d03.xy_loc.values[0].x, end_x_from_d04-w_cross_d03.xy_loc.values[0].x], 0,10000, colors = 'k', linestyle = '-',zorder = 50)
            if flight_leg in [4,5,6,7,8,9,10,11]: # Just use end of d04 for this one.
                ax.vlines([start_x_from_d04-w_cross_d03.xy_loc.values[0].x, xs_d04[-1]], 0,10000, colors = 'k', linestyle = '-',zorder = 50)

        # Add the color bar
        cb_dbz = fig.colorbar(dbz_contours_d03, ax=ax, label = 'w (m/s)')
        cb_dbz.ax.tick_params(labelsize=10)
        cb_dbz.ax.set_yticks(np.arange(-30,31,10))

        # Add a title
        if plot_d04:
            ax.set_title("Init {}\nFlight Leg {} d03 and d04".format(init_time_str,flight_leg), fontsize = 10, loc = 'left',)
        else:
            ax.set_title("Init {}\nFlight Leg {} d03".format(init_time_str,flight_leg), fontsize = 10, loc = 'left',)
        ax.set_title("Valid {}\n{} Track".format(valid_time_str, track), fontsize = 10, loc = 'right',)

        # Text to indicate potential temperature contours
        ax.text(0.03,0.01,'Dashed contours are QVAPOR in units g/kg',
                transform = ax.transAxes)

        # Save and show and close
        plt.savefig(save_path, dpi = 200, bbox_inches = 'tight')
        #plt.show()
        plt.close()
