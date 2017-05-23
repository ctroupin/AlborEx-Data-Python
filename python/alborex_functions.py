import os
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import six
from six.moves import xrange

import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.collections as mcollections
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png


def prepare_map(coordinates, res='i', proj='merc'):
    """Return a fig, m and ax objects
    given a set of coordinates defining a bounding box
    :param coordinates: list of coordinates (lonmin, lonmax, latmin, latmax)
    :param res: resolution in the projection ; 'i' by default (intermediate)
    :return: fig
    :type fig: Figure object
    :return m
    :type m: Basemap object
    :return ax
    :type ax: AxesSubplot object
    """
    m = Basemap(projection=proj,
                llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
                urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
                lat_ts=0.5 * (coordinates[2] + coordinates[3]), resolution=res)

    fig = plt.figure()
    ax = plt.subplot(111)
    m.ax = ax
    return fig, m, ax


def add_logo(imagepath, ax, position, zoom, zorder):
    """Add an image on the figure
    :param imagepath: path to the image to add
    :param ax: axes object
    :param position: relative position on the map
    :param zoom: zoom level
    :param zorder:
    :return:
    """
    logo2plot = read_png(imagepath)
    imagebox = OffsetImage(logo2plot, zoom=zoom)
    # coordinates to position this image

    ab = AnnotationBbox(imagebox, position,
                        xybox=(0., 0.),
                        xycoords='data',
                        pad=0.0,
                        boxcoords="offset points")
    ab.zorder = zorder
    ax.add_artist(ab)


def create_rect_patch(coordinates, m, **kwargs):
    """
    Create a rectangular patch to add on the map
    :param coordinates:
    :param m: Basemap object
    :return: patch
    """
    xr1, yr1 = m(coordinates[0], coordinates[2])
    xr2, yr2 = m(coordinates[0], coordinates[3])
    xr3, yr3 = m(coordinates[1], coordinates[3])
    xr4, yr4 = m(coordinates[1], coordinates[2])
    verts = [(xr1, yr1), (xr2, yr2), (xr3, yr3), (xr4, yr4), (xr1, yr1), ]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO,
             Path.LINETO, Path.CLOSEPOLY, ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, **kwargs)
    return patch


def extract_coastline(coordinates, filename, res='i'):
    """
    Extract the coastline in a region delimited by a bounding box
    and save the result in a text file
    :param coordinates: coordinates delimiting a bounding box (lonmin, lonmax, latmin, latmax)
    :param filename: name of the file where the coastline will be saved
    :param res: resolution for the extraction
    :return:
    """
    x, y = gshhs.get_coastline(xlim=[coordinates[0], coordinates[1]],
                               ylim=[coordinates[2], coordinates[3]],
                               res=res)
    # Save the coastline
    np.savetxt(filename, np.ma.vstack((x, y)).T)


def load_coast(coastfile, valex=-999.):
    """
    Read coastline coordinates from an existing file
    :param coastfile: name of the file containing the coastline
    :param valex: exclusion value
    :return: lon, lat
    """
    lon, lat = np.loadtxt(coastfile, usecols=(0, 1), unpack=True)
    lon[lon == valex] = np.nan
    lat[lat == valex] = np.nan
    return lon, lat


def load_coast_gshhs(coastfile, coordinates, valex=-999.):
    """
    Read coastline coordinates from an existing file in a region delimited
    by a bounding box
    :param coastfile: name of the file containing the coastline
    :param coordinates: coordinates delimiting a bounding box (lonmin, lonmax, latmin, latmax)
    :param valex: exclusion value
    :return: lon, lat
    """
    lon, lat = np.loadtxt(coastfile, usecols=(0, 1), unpack=True)
    goodcoord = np.where((lon >= coordinates[0]) & (lon <= coordinates[1]) &
                         (lat >= coordinates[2]) & (lat <= coordinates[3]))[0]
    goodcoord2 = np.where(np.logical_or((lon == valex), (lat == valex)))[0]
    goodcoord = np.union1d(goodcoord, goodcoord2)
    lon, lat = lon[goodcoord], lat[goodcoord]
    lon[lon == valex] = np.nan
    lat[lat == valex] = np.nan

    return lon, lat


def alborex_load_bathy(bathyfile, coordinates):
    """
    Load bathymetry from a netCDF file in a select region
    delimited by a list of coordinates
    :param bathyfile: name of the netCDF file
    :param coordinates: coordinates delimiting a bounding box (lonmin, lonmax, latmin, latmax)
    :return:
    """
    with netCDF4.Dataset( bathyfile, 'r') as nc:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
        depth = nc.variables['depth'][:]
        # subset
        goodlon = np.where(np.logical_and((lon >= coordinates[0]),
                                          (lon <= coordinates[1])))[0]
        goodlat = np.where(np.logical_and((lat >= coordinates[2]),
                                          (lat <= coordinates[3])))[0]
        lon = lon[goodlon]
        lat = lat[goodlat]
        depth = depth[goodlat, :]
        depth = depth[:, goodlon]
    return lon, lat, depth


def load_altimetry(altimetryfile, coordinates):
    """

    :param altimetryfile:
    :param coordinates: coordinates delimiting a bounding box (lonmin, lonmax, latmin, latmax)

    :return:
    """
    with netCDF4.Dataset(altimetryfile, 'r') as nc:
        lon = nc.variables['lon'][:] - 360.
        lat = nc.variables['lat'][:]
        u = np.squeeze(nc.variables['u'][:])
        v = np.squeeze(nc.variables['v'][:])
        # subset
        goodlon = np.where(np.logical_and((lon >= coordinates[0]),
                                          (lon <= coordinates[1])))[0]
        goodlat = np.where(np.logical_and((lat >= coordinates[2]),
                                          (lat <= coordinates[3])))[0]
        lon = lon[goodlon]
        lat = lat[goodlat]
        u = u[goodlat, :]
        u = u[:, goodlon]
        v = v[goodlat, :]
        v = v[:, goodlon]
    return lon, lat, u, v


def load_sst(sstfile, coordinates):
    """Return the coordinates, the SST, the time and the satellite name,
    given a data file and a list of coordinates delimiting a bounding box
    :param sstfile: name of the netCDF4 file containing the data
    :param coordinates: coordinates delimiting a bounding box (lonmin, lonmax, latmin, latmax)
    :return:
    """
    with netCDF4.Dataset(sstfile, 'r') as nc:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        timesst = nc.variables['time'][:]

        # subset
        goodlon = np.where(np.logical_and((lon >= coordinates[0]),
                                          (lon <= coordinates[1])))[0]
        goodlat = np.where(np.logical_and((lat >= coordinates[2]),
                                          (lat <= coordinates[3])))[0]
        lon = lon[goodlon]
        lat = lat[goodlat]
        sst = np.squeeze(nc.variables['mcsst'][:, goodlat, goodlon])
        mask = nc.variables['lsmask'][:, goodlat, goodlon].squeeze()
        sst = np.ma.masked_where(mask == 1, sst)

        timesst *= 60.
        sat = nc.satellite
        sensor = nc.sensor_name
        sensorsat = sensor.upper() + ' on ' + sat.upper()
    return lon, lat, sst, timesst, sensorsat


def load_sst_l2(filename):
    """
    Load the SST from netCDF L2 file obtained from
    https://oceancolor.gsfc.nasa.gov
    :param filename: name of the netCDF file
    :return: lon, lat, sst, sstflag, sstyear, sstday
    """
    if os.path.exists(filename):
        with netCDF4.Dataset(filename) as nc:
            # Read platform
            sat = nc.platform
            # Read time information
            # Assume all the measurements made the same day (and same year)
            year = nc.groups['scan_line_attributes'].variables['year'][0]
            day = nc.groups['scan_line_attributes'].variables['day'][0]
            # Read coordinates
            lon = nc.groups['navigation_data'].variables['longitude'][:]
            lat = nc.groups['navigation_data'].variables['latitude'][:]
            # Read geophysical variables
            try:
                sst = nc.groups['geophysical_data'].variables['sst'][:]
                sstqual = nc.groups['geophysical_data'].variables['qual_sst'][:]
            except KeyError:
                sst = nc.groups['geophysical_data'].variables['sst4'][:]
                sstqual = nc.groups['geophysical_data'].variables['qual_sst4'][:]
    else:
        lon, lat, sst, sstqual, year, day, sat = [], [], [], [], [], [], []
    return lon, lat, sst, sstqual, year, day, sat


def load_sst_l2_old(sstfile):
    """
    Load the SST from netCDF L2 file obtained from
    https://oceancolor.gsfc.nasa.gov
    :param sstfile: name of the netCDF file
    :return: lon, lat, sst, sstflag, sstyear, sstday
    """
    if 'SST4' in sstfile:
        sstname = 'Geophysical_Data_sst4'
        sstflagname = 'Geophysical_Data_qual_sst4'
    else:
        sstname = 'Geophysical_Data_sst'
        sstflagname = 'Geophysical_Data_qual_sst'

    with netCDF4.Dataset(sstfile, 'r') as nc:
        lon = nc.variables['Navigation_Data_longitude'][:]
        lat = nc.variables['Navigation_Data_latitude'][:]
        sst = nc.variables[sstname][:] * 0.005
        sstflag = nc.variables[sstflagname][:]
        sst = np.ma.masked_where(sstflag > 1, sst)
        sstyear = nc.Start_Year
        sstday = nc.Start_Day
    return lon, lat, sst, sstflag, sstyear, sstday


def load_ctd(ctdfile):
    """
    Load the coordinates (lon, lat, depth), the temperature and chlorophyll concentration
    from the selected netCDF file
    :param ctdfile: name of the netCDF file
    :return: lon, lat, depth, temp, chloro
    """
    with netCDF4.Dataset(ctdfile, 'r') as nc:
        lon = nc.variables['LON'][:]
        lat = nc.variables['LAT'][:]
        depth = nc.variables['DEPTH'][:]
        time = nc.variables['time'][:]
        temp = nc.variables['WTR_TEM_01'][:]
        chloro = nc.variables['CHLO'][:]
        chloro = np.ma.masked_where(np.isnan(chloro), chloro)
    return lon, lat, depth, time, temp, chloro


def load_glider_data(gliderfile, NN=1):
    """
    Load the coordinates and the temperature from a glider file
    :param gliderfile: name of the netCDF file
    :param NN: sub-sampling factor (keep 1 out of NN measurements)
    """
    with netCDF4.Dataset(gliderfile, 'r') as nc:
        lon = nc.variables['longitude'][::NN]
        lat = nc.variables['latitude'][::NN]
        depth = nc.variables['depth'][::NN]
        temperature = nc.variables['temperature'][::NN]
    return lon, lat, depth, temperature


def load_glider_coord(gliderfile):
    """
    Load the coordinates from a glider file
    :param gliderfile: name of the glider netCDF file
    :return: lon: longitude
    :return: lat: latitude
    :return: depth: depth
    :return: time: time
    """
    with netCDF4.Dataset(gliderfile, 'r') as nc:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
        depth = nc.variables['depth'][:]
        time = nc.variables['time'][:]
    return lon, lat, depth, time


def change_wall_prop(ax, coordinates, depths, angles):
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_xaxis.gridlines.set_linestyles(':')
    ax.w_yaxis.gridlines.set_linestyles(':')
    ax.w_zaxis.gridlines.set_linestyles(':')
    ax.view_init(angles[0], angles[1])
    ax.set_xlim(coordinates[0], coordinates[1])
    ax.set_ylim(coordinates[2], coordinates[3])
    ax.set_zlim(depths[0], depths[1])
    ax.set_zlabel('Depth (m)')

    ax.set_zticks(np.arange(depths[0], depths[1] + 10, depths[2]))
    ax.set_zticklabels(range(int(-depths[0]), -int(depths[1]) - 10, -int(depths[2])))


def read_l2_wind(windfile, coordinates):
    """
    Read the L2 wind from a netCDF file
    given a list of coordinates delimiting a bounding box
    :param windfile: netCDF file containing the data
    :param coordinates: coordinates delimiting a bounding box (lonmin, lonmax, latmin, latmax)
    :return: lon, lat, uwind, vwind, windtime
    """

    # Open NetCDF file
    with netCDF4.Dataset(windfile) as nc:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        windspeed = nc.variables['wind_speed'][:]
        winddirection = nc.variables['wind_dir'][:]
        windtime = nc.variables['time'][:]

    # Change longitudes
    lon[lon > 180] -= 360.0

    # Reduce dimensions
    lon = lon.flatten()
    lat = lat.flatten()
    windspeed = windspeed.flatten()
    winddirection = winddirection.flatten()

    # Select sub-region and check if data inside the area of interest
    goodlon = np.nonzero(np.logical_and(lon <= coordinates[1],
                                        lon >= coordinates[0]))
    goodlon = goodlon[0]

    if goodlon.size != 0:
        lat = lat[goodlon]
        lon = lon[goodlon]
        windspeed = windspeed[goodlon]
        winddirection = -winddirection[goodlon] + 90.
        goodlat = np.nonzero(np.logical_and(lat <= coordinates[3],
                                            lat >= coordinates[2]))
        goodlat = goodlat[0]
        if goodlat.size != 0:
            lat = lat[goodlat]
            lon = lon[goodlat]
            windspeed = windspeed[goodlat]
            winddirection = winddirection[goodlat]

            uwind = windspeed * np.cos(np.deg2rad(winddirection))
            vwind = windspeed * np.sin(np.deg2rad(winddirection))
            uwind = np.ma.masked_where(uwind == uwind.fill_value, uwind)
            vwind = np.ma.masked_where(vwind == vwind.fill_value, vwind)
            uwind.data[uwind.data == uwind.data.min()] = 0
            vwind.data[vwind.data == vwind.data.min()] = 0
        else:
            # print 'No value in selected region'
            # print ' '
            lon, lat, uwind, vwind = [], [], [], []
    else:
        print("No value in selected region")
        lon, lat, uwind, vwind = [], [], [], []

    return lon, lat, uwind, vwind, windtime
