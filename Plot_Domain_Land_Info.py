# Michael Wasserstein
# Plot_Domain_Land_Info.py
# 11/1/2024
# Script takes in WRF outputs and plots information about the domains being used, and the land surface use

####### Usage #########
# Conda environment - py37
# python /uufs/chpc.utah.edu/common/home/u1371671/WRF/Plotting_Scripts/Plot_Domain_Land_Info.py -r 2 -p 2
# -r represents the run number you want to plot
# -p is the path containing the wrf run

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import glob
import wrf
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim, ALL_TIMES)
from netCDF4 import Dataset
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime
import os, sys
import matplotlib as mpl
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import pyart
from matplotlib.colors import ListedColormap


######## User input arguments #############
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-r", "--run", help="WRF run of interest")
parser.add_argument("-p", "--path", help="Wrf path - where is data (1 or 2")
parser.add_argument("-H", "--hires", help="Is the terrain and land use high resolution from NCAR? (T or F)")

args = parser.parse_args()

# Get user inputs
run = str(args.run)
path = int(args.path)
hires = str(args.hires)
print('Plotting Domain for run', run)

run_number = '{}'.format(run).zfill(2) # Leading zeros included

# paths for data
if path ==1:
    base_path = '/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf/'
else:
    base_path = f'/uufs/chpc.utah.edu/common/home/steenburgh-group12/michael/wrf{path}/'
    
# Path where the WRF run is located
WRF_path = base_path + 'wrf_runs/wrf_{}/run/'.format(run_number)

# paths for saving fig
parent_dir = '/uufs/chpc.utah.edu/common/home/u1371671/'
Fig_dir = parent_dir + 'WRF/Figures_{}/wrf_{}/'.format(path,run_number)

# If you're using high resolution land use data, make that colormap
if hires == 'T':
    
    # Name list of colors
    C = np.array([
    [1,0,0],          #  1 Urban and Built-up Land
    [1,1,0],          #! 2 Dryland Cropland and Pasture
    [1,1,.2],         #  3 Irrigated Cropland and Pasture
    [1,1,.3],         #  4 Mixed Dryland/Irrigated Cropland and Pasture
    [.7,.9,.3],       #  5 Cropland/Grassland Mosaic
    [.7,.9,.3],       #  6 Cropland/Woodland Mosaic
    [0,1,0],          #  7 Grassland
    [.3,.7,0],        #  8 Shrubland
    [.82,.41,.12],    #  9 Mixed Shrubland/Grassland
    [1,.84,.0],       #  10 Savanna
    [.2,.8,.4],       #  11 Deciduous Broadleaf Forest
    [.2,.8,.2],       #  12 Deciduous Needleleaf Forest
    [0,.4,.2],        #  13 Evergreen Broadleaf Forest
    [0,.4,0],         #! 14 Evergreen Needleleaf Forest 
    [.2,.6,.2],       #  15 Mixed Forests
    [0,0,.88],        #  16 Water Bodies
    [0,1,1],          #! 17 Herbaceous Wetlands
    [.2,1,1],         #  18 Wooden Wetlands
    [.914,.914,.7],   #  19 Barren or Sparsely Vegetated
    [.86,.08,.23],    #  20 Herbaceous Tundraa
    [.86,.08,.23],    #  21 Wooded Tundra
    [.97,.5,.31],     #! 22 Mixed Tundra
    [.91,.59,.48],   #! 23 Barren Tundra
    [1,1,1],          #! 24 Snow and Ice
    [0.95686275, 0.64313725, 0.37647059], # 25 Playa
    [0.81176471, 0.0627451 , 0.1254902 ], # 26 Lava
    [1.        , 0.98039216, 0.94117647],  # 27 White Sand
    [0.8627451, 0.8627451, 0.8627451], # 28 Unassigned
    [0.8627451, 0.8627451, 0.8627451], # 29 Unassigned
    [0.8627451, 0.8627451, 0.8627451], # 30 Unassigned
    [1.        , 0.71372549, 0.75686275], # 31'Low Intensity Residential '
    [1., 0., 0.],  # 32  High Intensity Residential
    [0.50196078, 0.50196078, 0.50196078]  #33  Industrial or Commercial
    ])

    # Names of the categories
    kinds = ['Urban and Built-Up Land'
    ,'Dryland Cropland and Pasture'
    ,'Irrigated Cropland and Pasture'
    ,'Mixed Dryland/Irrigated Cropland and Pasture'
    ,'Cropland/Grassland Mosaic'
    ,'Cropland/Woodland Mosaic'
    ,'Grassland'
    ,'Shrubland'
    ,'Mixed Shrubland/Grassland'
    ,'Savanna'
    ,'Deciduous Broadleaf Forest'
    ,'Deciduous Needleleaf Forest'
    ,'Evergreen Broadleaf Forest'
    ,'Evergreen Needleleaf Forest'
    ,'Mixed Forest'
    ,'Water Bodies'
    ,'Herbaceous Wetland'
    ,'Wooded Wetland'
    ,'Barren or Sparsely Vegetated'
    ,'Herbaceous Tundra'
    ,'Wooded Tundra'
    ,'Mixed Tundra'
    ,'Bare Ground Tundra'
    ,'Snow or Ice'
    ,'Playa'
    ,'Lava'
    ,'White Sand'
    ,'Unassigned'
    ,'Unassigned'
    ,'Unassigned'
    ,'Low Intensity Residential '
    ,'High Intensity Residential'
    ,'Industrial or Commercial']

    # Array of the levels
    levs = np.arange(0.5, 34,1)
    
# If not high resolution land information
else:
    # More info here: https://ondemand.chpc.utah.edu/pun/sys/dashboard/files/fs//uufs/chpc.utah.edu/common/home/u1371671/steenburgh-group12/michael/wrf/WRF/run/LANDUSE.TBL
    C = np.array([
    [0,.4,0],      #  1 Evergreen Needleleaf Forest
    [0,.4,.2],      #! 2 Evergreen Broadleaf Forest    
    [.2,.8,.2],     #  3 Deciduous Needleleaf Forest
    [.2,.8,.4],     #  4 Deciduous Broadleaf Forest
    [.2,.6,.2],     #  5 Mixed Forests
    [.3,.7,.1],      #  6 Closed Shrublands
    [.82,.41,.12],     #  7 Open Shurblands
    [.74,.71,.41],       #  8 Woody Savannas
    [1,.84,.0],     #  9 Savannas
    [0,1,.15],        #  10 Grasslands
    [.06,1,1],        #! 11 Permanant Wetlands
    [1,.95,.05],      #  12 Croplands
    [1,.05,0],     #  13 Urban and Built-up
    [.7,.9,.3],      #! 14 Cropland/Natual Vegation Mosaic
    [1,1,1],        #! 15 Snow and Ice
    [.914,.914,.7], #  16 Barren or Sparsely Vegetated
    [.5,.7,1],        #  17 Water (like oceans)
    [.86,.08,.23],        #  18 Wooded Tundra
    [.97,.5,.31],        #! 19 Mixed Tundra
    [.91,.59,.49],     #! 20 Barren Tundra
    [0,0,.88]])      #! 21 Lake
    
    # Names of the categories
    kinds = ['Evergreen Needleleaf Forest',
                  'Evergreen Broadleaf Forest',
                  'Deciduous Needleleaf Forest',
                  'Deciduous Broadleaf Forest',
                  'Mixed Forests',
                  'Closed Shrublands',
                  'Open Shrublands',
                  'Woody Savannas',
                  'Savannas',
                  'Grasslands',
                  'Permanent Wetlands',
                  'Croplands',
                  'Urban and Built-Up',
                  'Cropland/Natural Vegetation Mosaic',
                  'Snow and Ice',
                  'Barren or Sparsely Vegetated',
                  'Water',
                  'Wooded Tundra',
                  'Mixed Tundra',
                  'Barren Tundra',
                  'Lake'] 
    



    
    # Array of the levels
    levs = np.arange(0.5, 22,1)

# Make land use colormap
cm = ListedColormap(C)


# load in all the wrf output data files
data_files_d01 = glob.glob(WRF_path + '*wrfout_d01*') # for the outermost domain
data_files_d02 = glob.glob(WRF_path + '*wrfout_d02*') # for the 2nd domain
data_files_d03 = glob.glob(WRF_path + '*wrfout_d03*') # for the innermost domain

# Get a netcdf file for each domain
ncfile_d01 = Dataset(data_files_d01[10])
ncfile_d02 = Dataset(data_files_d02[10])
ncfile_d03 = Dataset(data_files_d03[10])

# Get the grid spaing in units of meters
d_01_gridspacing = ncfile_d01.DX / 1000
d_02_gridspacing = ncfile_d02.DX / 1000
d_03_gridspacing = ncfile_d03.DX / 1000

# Get the terrain heights
ter_d01 = getvar(ncfile_d01, "ter")
ter_d02 = getvar(ncfile_d02, "ter")
ter_d03 = getvar(ncfile_d03, "ter")

# Get land info
LANDMASK_d01 = wrf.getvar(ncfile_d01, 'LANDMASK', timeidx=-1)
LANDMASK_d02 = wrf.getvar(ncfile_d02, 'LANDMASK', timeidx=-1)
LANDMASK_d03 = wrf.getvar(ncfile_d03, 'LANDMASK', timeidx=-1)

# Get land use info
LU_INDEX_d01 = wrf.getvar(ncfile_d01, 'LU_INDEX', timeidx=-1)
LU_INDEX_d02 = wrf.getvar(ncfile_d02, 'LU_INDEX', timeidx=-1)
LU_INDEX_d03 = wrf.getvar(ncfile_d03, 'LU_INDEX', timeidx=-1)

# Projection information
cart_proj = get_cartopy(wrfin=ncfile_d01)  # projection will be the same for all domains, so only do for one domain

# Limits information
xlim_d01, ylim_d01 = cartopy_xlim(wrfin=ncfile_d01), cartopy_ylim(wrfin=ncfile_d01)
xlim_d02, ylim_d02 = cartopy_xlim(wrfin=ncfile_d02), cartopy_ylim(wrfin=ncfile_d02)
xlim_d03, ylim_d03 = cartopy_xlim(wrfin=ncfile_d03), cartopy_ylim(wrfin=ncfile_d03)

# Get information about times
init_time = ncfile_d03.SIMULATION_START_DATE
init_time = datetime.datetime.strptime(init_time, '%Y-%m-%d_%H:%M:%S')
init_time_str = datetime.datetime.strftime(init_time, '%b %-d, %Y %H:%M UTC')

################## Plot all three domains in one plot #################
# Create figure and axis
fig = plt.figure(figsize = (8, 8), facecolor = 'white', edgecolor = 'k')
ax1 = plt.axes(projection=cart_proj)

# Add coasts and important borders
ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
ax1.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)

##################################
# set limits for entire plot (based on d01)
ax1.set_xlim([xlim_d01[0]-(xlim_d01[1]-xlim_d01[0])/15, xlim_d01[1]+(xlim_d01[1]-xlim_d01[0])/15])
ax1.set_ylim([ylim_d01[0]-(ylim_d01[1]-ylim_d01[0])/15, ylim_d01[1]+(ylim_d01[1]-ylim_d01[0])/15])
 
# d01 bounding box
ax1.add_patch(mpl.patches.Rectangle((xlim_d01[0], ylim_d01[0]), xlim_d01[1]-xlim_d01[0], ylim_d01[1]-ylim_d01[0],
             fill=None, lw=3, edgecolor='black', zorder=10))
ax1.text(xlim_d01[0]+(xlim_d01[1]-xlim_d01[0])*0.05, ylim_d01[0]+(ylim_d01[1]-ylim_d01[0])*1.02, '{} km'.format(d_01_gridspacing),
        size=15, color='black', zorder=10)

# Plot d01 terrain
terr = ax1.contourf(ter_d01.XLONG.values, ter_d01.XLAT.values, ter_d01.values, cmap = 'terrain', 
                    transform=ccrs.PlateCarree(), levels = np.arange(0,4001,125), extend = 'both')

##################################
# d02 bounding box
ax1.add_patch(mpl.patches.Rectangle((xlim_d02[0], ylim_d02[0]), xlim_d02[1]-xlim_d02[0], ylim_d02[1]-ylim_d02[0],
             fill=None, lw=3, edgecolor='black', zorder=10))
ax1.text(xlim_d02[0]+(xlim_d02[1]-xlim_d02[0])*0.05, ylim_d02[0]+(ylim_d02[1]-ylim_d02[0])*1.1, '{} km'.format(d_02_gridspacing),
        size=15, color='black', zorder=10)

# Plot d02 terrain
terr = ax1.contourf(ter_d02.XLONG.values, ter_d02.XLAT.values, ter_d02.values, cmap = 'terrain', 
                    transform=ccrs.PlateCarree(), levels = np.arange(0,4001,125), extend = 'both')

##################################
# d03 bounding box
ax1.add_patch(mpl.patches.Rectangle((xlim_d03[0], ylim_d03[0]), xlim_d03[1]-xlim_d03[0], ylim_d03[1]-ylim_d03[0],
             fill=None, lw=3, edgecolor='black', zorder=10))
ax1.text(xlim_d03[0], ylim_d03[0]+(ylim_d03[1]-ylim_d03[0])*1.05, '{} km'.format(d_03_gridspacing),
        size=13, color='black', zorder=10)

# Plot d03 terrain
ter = ax1.contourf(ter_d03.XLONG.values, ter_d03.XLAT.values, ter_d03.values, cmap = 'terrain', 
                    transform=ccrs.PlateCarree(), levels = np.arange(0,4001,125), extend = 'both')

# Add titles
ax1.set_title('WRF{} run {}'.format(path, run_number), loc = 'right')
ax1.set_title('Init {}'.format( init_time_str), loc = 'left')

# Add colorbar
plt.colorbar(ter,  orientation = 'horizontal', label = 'Terrain Height (m)', ax = ax1, fraction = 0.05, pad = 0.05)

plt.savefig(Fig_dir +'Domains_WRF_terrain.png', dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()

################## Plot domain 1 and lake #################
# Generate figure and axis
fig = plt.figure(figsize = (8, 8))
ax1 = plt.axes(projection=cart_proj)

# Add coasts and important borders
ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
ax1.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)

# d01 set limits
ax1.set_xlim([xlim_d01[0], xlim_d01[1]])
ax1.set_ylim([ylim_d01[0], ylim_d01[1]])

# d01 box
ax1.add_patch(mpl.patches.Rectangle((xlim_d01[0], ylim_d01[0]), xlim_d01[1]-xlim_d01[0], ylim_d01[1]-ylim_d01[0],
             fill=None, lw=3, edgecolor='black', zorder=10))

# Plot terrain height
plot = ax1.contourf(ter_d01.XLONG.values, ter_d01.XLAT.values, ter_d01.values, cmap = 'terrain', 
                    transform=ccrs.PlateCarree(), levels = np.arange(500,4001,125), extend = 'both')
# Plot lakes
plot2 = ax1.contourf(LANDMASK_d01.XLONG.values, LANDMASK_d01.XLAT.values, LANDMASK_d01.values,
                    transform=ccrs.PlateCarree(), levels = np.arange(0,0.6,0.5), colors = 'midnightblue', zorder = 100) # between 

# Add colorbar
plt.colorbar(plot, orientation = 'horizontal', label = 'Terrain Height (m)', ax = ax1, fraction = 0.05, pad = 0.05)

# Add titles
ax1.set_title('WRF{} run {} - d01'.format(path, run_number), loc = 'right')
ax1.set_title('Init {}'.format( init_time_str), loc = 'left')

# Save figure
plt.savefig(Fig_dir + 'Domain01_WRF_terrain.png', dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()

################## Plot domain 2 and lake #################
# Generate figure and axis
fig = plt.figure(figsize = (8, 8))
ax1 = plt.axes(projection=cart_proj)

# Add coasts and important borders
ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
ax1.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)

# d02 set limits
ax1.set_xlim([xlim_d02[0], xlim_d02[1]])
ax1.set_ylim([ylim_d02[0], ylim_d02[1]])

# d02 box
ax1.add_patch(mpl.patches.Rectangle((xlim_d02[0], ylim_d02[0]), xlim_d02[1]-xlim_d02[0], ylim_d02[1]-ylim_d02[0],
             fill=None, lw=3, edgecolor='black', zorder=10))

# Plot terrain height
plot = ax1.contourf(ter_d02.XLONG.values, ter_d02.XLAT.values, ter_d02.values, cmap = 'terrain', 
                    transform=ccrs.PlateCarree(), levels = np.arange(500,4001,125), extend = 'both')
# Plot lakes
plot2 = ax1.contourf(LANDMASK_d02.XLONG.values, LANDMASK_d02.XLAT.values, LANDMASK_d02.values,
                    transform=ccrs.PlateCarree(), levels = np.arange(0,0.6,0.5), colors = 'midnightblue', zorder = 100) # between 

# Add colorbar
plt.colorbar(plot, orientation = 'horizontal', label = 'Terrain Height (m)', ax = ax1, fraction = 0.05, pad = 0.05)

# Add titles
ax1.set_title('WRF{} run {} - d02'.format(path, run_number), loc = 'right')
ax1.set_title('Init {}'.format( init_time_str), loc = 'left')

# Save figure
plt.savefig(Fig_dir + 'Domain02_WRF_terrain.png', dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()

################## Plot domain 3 and lake #################
# Generate figure and axis
fig = plt.figure(figsize = (8, 8))
ax1 = plt.axes(projection=cart_proj)

# Add coasts and important borders
ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
ax1.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)

# d03 set limits
ax1.set_xlim([xlim_d03[0], xlim_d03[1]])
ax1.set_ylim([ylim_d03[0], ylim_d03[1]])

# d03 box
ax1.add_patch(mpl.patches.Rectangle((xlim_d03[0], ylim_d03[0]), xlim_d03[1]-xlim_d03[0], ylim_d03[1]-ylim_d03[0],
             fill=None, lw=3, edgecolor='black', zorder=10))

# Plot terrain height
plot = ax1.contourf(ter_d03.XLONG.values, ter_d03.XLAT.values, ter_d03.values, cmap = 'terrain', 
                    transform=ccrs.PlateCarree(), levels = np.arange(500,4001,125), extend = 'both')
# Plot lakes
plot2 = ax1.contourf(LANDMASK_d03.XLONG.values, LANDMASK_d03.XLAT.values, LANDMASK_d03.values,
                    transform=ccrs.PlateCarree(), levels = np.arange(0,0.6,0.5), colors = 'midnightblue', zorder = 100) # between 

# Add colorbar
plt.colorbar(plot, orientation = 'horizontal', label = 'Terrain Height (m)', ax = ax1, fraction = 0.05, pad = 0.05)

# Add titles
ax1.set_title('WRF{} run {} - d03'.format(path, run_number), loc = 'right')
ax1.set_title('Init {}'.format( init_time_str), loc = 'left')

# Save figure
plt.savefig(Fig_dir + 'Domain03_WRF_terrain.png', dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()

######################### Plot land use for d03 ############################
# Create figure and axis
fig = plt.figure(figsize = (8, 8), facecolor = 'white', edgecolor = 'k',)
ax1 = plt.axes(projection=cart_proj)

# Add coasts and important borders
ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5, zorder = 100)
ax1.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.25, zorder = 100)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.25, zorder = 100)


# d03 set limits
ax1.set_xlim([xlim_d03[0], xlim_d03[1]])
ax1.set_ylim([ylim_d03[0], ylim_d03[1]])
 
# d03 box
ax1.add_patch(mpl.patches.Rectangle((xlim_d03[0], ylim_d03[0]), xlim_d03[1]-xlim_d03[0], ylim_d03[1]-ylim_d03[0],
             fill=None, lw=3, edgecolor='black', zorder=10))

# Land use contour fill
plot2 = ax1.contourf(LU_INDEX_d03.XLONG.values, LU_INDEX_d03.XLAT.values, LU_INDEX_d03.values,
                    transform=ccrs.PlateCarree(), zorder = 100, boundaries = levs, levels = levs, cmap = cm) # between 

# Add colorbar
cb = plt.colorbar(plot2, orientation = 'horizontal', ax = ax1, pad = 0.05,  fraction = 0.05)
cb.ax.set_xticks((levs + 0.5)[:-1])
cb.ax.set_xticklabels(kinds,rotation = 90)

# Add titles
ax1.set_title('WRF{} run {} - d03'.format(path, run_number), loc = 'right')
ax1.set_title('Init {}'.format( init_time_str), loc = 'left')

plt.savefig(Fig_dir + 'LU_INDEX_d03.png', dpi = 300, bbox_inches = 'tight')
#plt.show()
plt.close()