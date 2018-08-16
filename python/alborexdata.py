import os
import netCDF4
import logging
import datetime
import numpy as np
import seawater
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.patches as patches
from matplotlib.path import Path
from mpl_toolkits.mplot3d import Axes3D
from geopy.distance import vincenty
import cmocean

import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

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

def configure_logging(logfile="./alborexFig2.log"):
    """Prepare the logging messages and file
    """
    logger = logging.getLogger("alborex_logger")
    logger.setLevel(logging.DEBUG)
    # Format for our loglines
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # Setup console logging
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # Setup file logging as well
    fh = logging.FileHandler(logfile)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger


def add_map_grid(m, coordinates, dlon, dlat, **kwargs):
    """Add x and y ticks (no line plotted for better visibility)
    """
    m.drawparallels(np.arange(round(coordinates[2]), coordinates[3], dlat), labels=[1, 0, 0, 0], **kwargs)
    m.drawmeridians(np.arange(round(coordinates[0]), coordinates[1], dlon), labels=[0, 0, 0, 1], **kwargs)


def load_lonloat_ctdleg(datafile):
    """Return coordinates from the file containing the information
    on the different legs
    """
    lon, lat = [], []
    with open(datafile) as f:
        line = f.readline().rsplit()
        while line:
            # print(line)
            lon.append(float(line[2]))
            lat.append(float(line[3]))
            line = f.readline().rsplit()
    return lon, lat


def read_lonlat_coast(filename, valex=999):
    """Return the coordinates of the contours
    as a list of lists (one list per contour)
    """
    with open(filename) as f:
        lonall, latall = [], []
        lon, lat = [], []
        line = f.readline().rsplit()
        while line:
            if float(line[0]) == valex:
                lonall.append(lon)
                latall.append(lat)
                lon, lat = [], []
            else:
                lon.append(float(line[0]))
                lat.append(float(line[1]))
            line = f.readline().rsplit()
    return lonall, latall


class Front(object):

    def __init__(self, lon=None, lat=None):
        self.lon = lon
        self.lat = lat

    def get_from_file(self, filename):
        """
        Read the coordinates from a text file (lon, lat)
        :param filename: file name
        :type filename: str
        """
        self.lon = []
        self.lat = []
        if os.path.exists(filename):
            with open(filename, "r") as df:
                for lines in df.readlines():
                    self.lon.append(float(lines.rstrip().split(',')[0]))
                    self.lat.append(float(lines.rstrip().split(',')[1]))

    def smooth(self, n=4, s=0.01, nest=4):
        """

        :param N: subsampling factor
        :param s: smoothness parameter
        :param nest: estimate of number of knots needed (-1 = maximal)
        :return:
        """
        npoints = len(self.lon)
        if npoints > 0:
            if npoints == len(self.lat):
                t = np.linspace(0, 1, npoints)
                t2 = np.linspace(0, 1, n * npoints)

                # find the knot points
                tckp, u = interpolate.splprep([t, self.lon, self.lat], s=s, nest=-1)

                # evaluate spline, including interpolated points
                xnew, self.lon, self.lat = interpolate.splev(t2, tckp)


class Drifter(object):

    def __init__(self, lon=None, lat=None, time=None, temperature=None,
                 qclon=None, qclat=None):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.temperature = temperature
        self.qclon = qclon
        self.qclat = qclat
        self.timeunits = None
        self.dates = None
        self.velocity = None
        self.distance2front = None

    def get_from_netcdf(self, datafile):
        """
        Read the coordinates and the temperature from existing data file
        """
        if os.path.exists(datafile):
            with netCDF4.Dataset(datafile, 'r') as nc:
                self.lon = nc.get_variables_by_attributes(standard_name='longitude')[0][:]
                self.lat = nc.get_variables_by_attributes(standard_name='latitude')[0][:]
                self.time = nc.get_variables_by_attributes(standard_name='time')[0][:]
                self.timeunits = nc.get_variables_by_attributes(standard_name='time')[0].units
                self.dates = netCDF4.num2date(self.time, self.timeunits)

                try:
                    self.qclat = nc.get_variables_by_attributes(standard_name='latitude status_flag')[0][:]
                except IndexError:
                    self.qclat = None

                try:
                    self.qclon = nc.get_variables_by_attributes(standard_name='longitude status_flag')[0][:]
                except IndexError:
                    self.qclon = None
                try:
                    self.temperature = nc.get_variables_by_attributes(standard_name='sea_water_temperature')[0][:]
                except IndexError:
                    self.temperature = None

    def apply_qc_latlon(self, QC=[1]):
        """
        Discard the measurements of which the position
        doesn't have the indicated quality flag
        """
        if (self.qclon is not None) and (self.qclat is not None):
            badlon = [qc not in QC for qc in self.qclon]
            badlat = [qc not in QC for qc in self.qclat]
            badposition = np.logical_or(np.array(badlon), np.array(badlat))
            self.lon = np.ma.masked_where(badposition, self.lon)
            self.lat = np.ma.masked_where(badposition, self.lat)

    def mask_temp(self, tmin, tmax):
        if self.temperature is not None:
            self.temperature = np.ma.masked_outside(self.temperature,
                                                    tmin,
                                                    tmax,
                                                    copy=True)

    def select_dates(self, finaldate, initialdate=None):
        """
        Mask the time outside the selected period

        finaldate and initialdate are `datetime` obects
        for example: finaldate=datatime.datetime(2017, 5, 3, 18, 30, 0)
        """
        if initialdate is not None:
            self.lon = np.ma.masked_where(np.logical_or(self.dates > finaldate,
                                                        self.dates < initialdate),
                                                        self.lon)
            self.lat = np.ma.masked_where(np.logical_or(self.dates > finaldate,
                                                        self.dates < initialdate),
                                                        self.lat)
        else:
            self.lon = np.ma.masked_where(self.dates > finaldate, self.lon)
            self.lat = np.ma.masked_where(self.dates > finaldate, self.lat)

    def scatter_plot(self, m, **kwargs):
        scat = m.scatter(self.lon, self.lat, c=self.temperature, latlon=True, **kwargs)
        return scat

    def point_plot(self, m, **kwargs):
        m.plot(self.lon.compressed(), self.lat.compressed(), latlon=True, **kwargs)

    def add_initial_position(self, m, **kwargs):
        m.plot(self.lon[0], self.lat[0], latlon=True, linewidth=0, **kwargs)

    def compute_velocity(self, velmax=5.):
        """
        Compute the velocity using the Vincenty distance

        The values above velmax are masked
        """

        distancevec = np.zeros(len(self.lon)-1)
        timevec = self.time[1:] - self.time[:-1]
        for ii in range(0, len(self.lon)-1):
            distancevec[ii] = vincenty((self.lat[ii+1], self.lon[ii+1]),
                                       (self.lat[ii], self.lon[ii])).m
        self.velocity = distancevec / timevec
        self.velocity = np.ma.masked_greater(self.velocity, velmax, copy=True)

    def get_distance_front(self, frontlon, frontlat):
        """
        For each position of the drifter, compute the distance to the front,
        specified by 2 arrays of longitudes and latitudes

        **Note:**
        Brute force approach but could also approximate the front by a parabola
        and use the formula to get the distance.
        """
        npoints = len(frontlon)
        distance2front = np.zeros(len(self.lon))
        jj = 0
        for lond, latd in zip(self.lon, self.lat):
            dd = np.zeros(npoints)
            ii = 0
            for lonf, latf in zip(frontlon, frontlat):
                dd[ii] = vincenty((lonf, latf), (lond, latf)).m
                ii += 1
            distance2front[jj] = np.min(dd)
            jj += 1
        self.distance2front = distance2front


class Thermosal(object):
    """
    Thermosalinograph (temperature and salinity measured by the
    ship near the surface)
    """
    def __init__(self, lon=None, lat=None, time=None,
                 temperature=None, salinity=None, qclon=None, qclat=None,
                 qctemp=None, qcsal=None):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.temperature = temperature
        self.salinity = salinity

    def get_from_netcdf(self, datafile):
        """
        Read the coordinates and the field values from a netCDF file
        """
        if os.path.exists(datafile):
            with netCDF4.Dataset(datafile, 'r') as nc:
                self.lon = nc.get_variables_by_attributes(standard_name='longitude')[0][:]
                self.lat = nc.get_variables_by_attributes(standard_name='latitude')[0][:]
                self.time = nc.get_variables_by_attributes(standard_name='time')[0][:]
                timeunits = nc.get_variables_by_attributes(standard_name='time')[0].units
                self.dates = netCDF4.num2date(self.time, timeunits)
                self.salinity = nc.get_variables_by_attributes(standard_name='sea_water_salinity')[0][:]
                self.temperature = nc.get_variables_by_attributes(standard_name='sea_water_temperature')[0][:]

class CTD():

    def __init__(self, lon=None, lat=None, time=None, depth=None, pressure=None,
                 temperature=None, salinity=None, qclon=None, qclat=None,
                 qctemp=None, qcsal=None, chloro=None):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.depth = depth
        self.pressure = pressure
        self.temperature = temperature
        self.salinity = salinity
        self.qclon = qclon
        self.qclat = qclat
        self.qctemp = qctemp
        self.qcsal = qcsal
        self.timeunits = None
        self.dates = None
        self.chloro = chloro

    def get_from_netcdf(self, datafile):
        """
        Read the coordinates and the temperature from existing data file
        """
        if os.path.exists(datafile):
            with netCDF4.Dataset(datafile, 'r') as nc:

                try:
                    self.pressure =  nc.get_variables_by_attributes(standard_name='sea_water_pressure')[0][:]
                except IndexError:
                    self.pressure = None

                self.lon = nc.get_variables_by_attributes(standard_name='longitude')[0][:]
                self.lat = nc.get_variables_by_attributes(standard_name='latitude')[0][:]
                self.depth = nc.get_variables_by_attributes(standard_name='depth')[0][:]
                self.time = nc.get_variables_by_attributes(standard_name='time')[0][:]
                self.timeunits = nc.get_variables_by_attributes(standard_name='time')[0].units
                self.dates = netCDF4.num2date(self.time, self.timeunits)

                try:
                    self.chloro = nc.variables["CHLO"][:]
                except KeyError:
                    self.chloro = None

                try:
                    self.qclat = nc.get_variables_by_attributes(standard_name='latitude status_flag')[0][:]
                except IndexError:
                    self.qclat = None

                try:
                    self.qclon = nc.get_variables_by_attributes(standard_name='longitude status_flag')[0][:]
                except IndexError:
                    self.qclon = None

                # Get salinity
                try:
                    salinityvar = nc.get_variables_by_attributes(standard_name='sea_water_practical_salinity')[0]
                    salinityqcvar = salinityvar.ancillary_variables
                    self.salinity = salinityvar[:]
                    self.qcsal = nc.variables[salinityqcvar][:]
                except IndexError:
                    try:
                        salinityvar = nc.get_variables_by_attributes(standard_name='sea_water_salinity')[0]
                        self.salinity = salinityvar[:]
                        salinityqcvar = salinityvar.ancillary_variables
                        self.qcsal = nc.variables[salinityqcvar][:]
                    except AttributeError:
                        self.qcsal = None


                # Get (potential) temperature and convert if necessary
                try:
                    tempvar = nc.get_variables_by_attributes(standard_name='sea_water_temperature')[0]
                    self.temperature = tempvar[:]

                except IndexError:
                    try:
                        tempvar = nc.get_variables_by_attributes(standard_name='sea_water_potential_temperature')[0]
                        potentialtemp = tempvar[:]
                        self.temperature = seawater.temp(self.salinity, potentialtemp, self.pressure)

                    except IndexError:
                        self.temperature = None
                        self.qctemp = None

                try:
                    tempqcvar = tempvar.ancillary_variables
                    self.qctemp = nc.variables[tempqcvar][:]
                except AttributeError:
                    self.qctemp = None

class Glider(CTD):

    def remove_masked_coords(self):
        """
        Remove the masked coordinates (lon, lat, time, dates)
        """
        coordmask = np.logical_not(self.lon.mask)
        self.time = self.time.compress(coordmask)
        self.dates = self.dates.compress(coordmask)
        self.lon = self.lon.compressed()
        self.lat = self.lat.compressed()

    def get_day_indices(self, ndays=1):
        """
        Get the time indices corresponding to the start of days,
        separated by "ndays"
        """
        day_indices = []
        date_list = []

        # Convert the time to datses
        datestart, dateend = self.dates[0], self.dates[-1]
        date = datetime.datetime(datestart.year, datestart.month, datestart.day,
                                 0, 0, 0)

        while date <= dateend:
            # Increment initial date
            date += datetime.timedelta(days=ndays)
            date_list.append(date)
            # Get corresponding index
            index = np.argmin(abs(self.time - netCDF4.date2num(date, self.timeunits)))
            day_indices.append(index)

        return day_indices, date_list

    def scatter_plot(self, ax, **kwargs):
        """
        Add the measurements to a 3D scatter plot
        """
        scat3D = ax.scatter(self.lon, self.lat, -self.depth, **kwargs)
        return scat3D

class Profiler(CTD):
    """
    Stores Argo profiler data
    """
    def select_dates(self, finaldate, initialdate=None):
        """
        Mask the time outside the selected period

        finaldate and initialdate are `datetime` obects
        for example: finaldate=datatime.datetime(2017, 5, 3, 18, 30, 0)
        """
        if initialdate is not None:
            dates2mask = np.logical_or(self.dates > finaldate,
                                       self.dates < initialdate)
        else:
            dates2mask = self.dates > finaldate

        ndepth = self.depth.shape[1]
        dates2mask2D = np.matlib.repmat(dates2mask, ndepth, 1).transpose()
        self.lon = np.ma.masked_where(dates2mask, self.lon)
        self.lat = np.ma.masked_where(dates2mask, self.lat)
        self.dates = np.ma.masked_where(dates2mask, self.dates)
        self.depth = np.ma.masked_where(dates2mask2D, self.depth)
        self.temperature = np.ma.masked_where(dates2mask2D, self.temperature)
        self.salinity = np.ma.masked_where(dates2mask2D, self.salinity)

class Ship(Drifter):

    def apply_qc(self, qflag=1):
        """
        Mask the coordinates with a quality flag different from the specified value
        1 = good data
        """
        badcoords = np.logical_or(self.qclon != 1, self.qclat !=1)
        self.lon = np.ma.masked_where(badcoords, self.lon)
        self.lat = np.ma.masked_where(badcoords, self.lat)

    def plot_track(self, m, **kwargs):
        m.plot(self.lon, self.lat, latlon=True, **kwargs)




class SST(object):
    """
    Sea surface temperature field
    """

    def __init__(self, lon=None, lat=None, field=None, qflag=None,
                 year=None, dayofyear=None):
        self.lon = lon
        self.lat = lat
        self.field = field
        self.qflag = qflag
        self.timeunits = year
        self.year = year
        self.dayofyear = dayofyear

    def read_from_oceancolorL2(self, filename):
        """
        Load the SST from netCDF L2 file obtained from
        https://oceancolor.gsfc.nasa.gov
        :param filename: name of the netCDF file
        :return: lon, lat, field, qflag, year, dayofyear
        """

        if os.path.exists(filename):
            with netCDF4.Dataset(filename) as nc:
                # Read platform
                sat = nc.platform
                # Read time information
                # Assume all the measurements made the same day (and same year)
                self.year = nc.groups['scan_line_attributes'].variables['year'][0]
                self.dayofyear = nc.groups['scan_line_attributes'].variables['day'][0]
                # Read coordinates
                self.lon = nc.groups['navigation_data'].variables['longitude'][:]
                self.lat = nc.groups['navigation_data'].variables['latitude'][:]
                # Read geophysical variables
                try:
                    self.field = nc.groups['geophysical_data'].variables['sst'][:]
                    self.qflag = nc.groups['geophysical_data'].variables['qual_sst'][:]
                except KeyError:
                    self.field = nc.groups['geophysical_data'].variables['sst4'][:]
                    self.qflag = nc.groups['geophysical_data'].variables['qual_sst4'][:]

    def apply_qc(self, qf=1):
        """
        Mask the sst values which don't match the mentioned quality flag
        """
        self.field = np.ma.masked_where(self.qflag != 1, self.field)

class Adcp(object):

    """
    Stores ADCP transects
    """

    def __init_(self, lon=None, lat=None, depth=None,
                u=None, v=None, qcu=None, qcv=None,
                time=None, dates=None):
        self.lon = lon
        self.lat = lat
        self.depth = depth
        self.u = u
        self.v = v
        self.qclon = qclon
        self.qclat = qclat
        self.qcu = qcu
        self.qcv = qcv
        self.time = time
        self.dates = dates

    def get_from_netcdf(self, filename):
        """
        Read the coordinates and the velocity components
        from the netCDF file
        """
        if os.path.exists(filename):
            with netCDF4.Dataset(filename) as nc:
                self.lon = nc.get_variables_by_attributes(standard_name='longitude')[0][:]
                self.lat = nc.get_variables_by_attributes(standard_name='latitude')[0][:]
                self.depth = nc.get_variables_by_attributes(standard_name='depth')[0][:]
                self.time = nc.get_variables_by_attributes(standard_name='time')[0][:]
                self.timeunits = nc.get_variables_by_attributes(standard_name='time')[0].units
                self.dates = netCDF4.num2date(self.time, self.timeunits)
                self.qclat = nc.get_variables_by_attributes(standard_name='latitude status_flag')[0][:]
                self.qclon = nc.get_variables_by_attributes(standard_name='longitude status_flag')[0][:]
                # Velocity components
                uvar = nc.get_variables_by_attributes(standard_name='eastward_sea_water_velocity')[0]
                vvar = nc.get_variables_by_attributes(standard_name='northward_sea_water_velocity')[0]
                self.u = uvar[:]
                self.v = vvar[:]
                # Quality flags for velocity
                uqcvar = uvar.ancillary_variables
                vqcvar = vvar.ancillary_variables
                self.qcu = nc.variables[uqcvar][:]
                self.qcv = nc.variables[uqcvar][:]

    def get_from_matfile(self, filename):
        """
        Read the coordinates (lon, lat, depth) and
        the velocity components from the .mat files
        """
        # Read the mat file
        dataadcp = sio.loadmat(filename)

        self.lon = dataadcp["AnFLonDeg"]
        self.lat = dataadcp["AnFLatDeg"]
        self.u = dataadcp["SerEmmpersec"]
        self.v = dataadcp["SerNmmpersec"]
        ndepth = self.u.shape[1]
        depthmin = 16.
        deltadepth = 8.
        depthmax = depthmin + (ndepth - 1) * deltadepth
        self.depth = np.linspace(depthmin, depthmax, int(nbins))

    def apply_qc(self, qf=1):
        """
        Mask the velocity values which don't match the mentioned quality flag
        """
        self.u = np.ma.masked_where(self.qcu != 1, self.u)
        self.v = np.ma.masked_where(self.qcv != 1, self.v)

    def get_norm(self):
        """
        Compute the norm of the velocity vectors
        """
        self.velnorm = np.sqrt(self.u * self.u + self.v * self.v)

    def get_time_index(self, datemin=None, datemax=None):
        """
        Return an array of indices corresponding to the dates between
        datemin and datemax
        """
        if datemin is not None:
            if datemax is not None:
                gooddates = np.where( (self.dates >= datemin) and (self.dates <= datemax))[0]
            else:
                gooddates = np.where(self.dates >= datemin)[0]
        else:
            if datemax is not None:
                gooddates = np.where(self.dates <= datemax)[0]
            else:
                gooddates = np.where(self.dates)[0]

        return gooddates

    def plot_adcp_quiver(self, m, depthindex=0, depth=None, datemin=None, datemax=None):
        """
        Plot velocity field with arrows on a map
        """

        gooddates = self.get_time_index(datemin, datemax)

        m.plot(self.lon[gooddates], self.lat[gooddates], "k--", lw=.2, latlon=True)
        llon, llat = m(self.lon[gooddates], self.lat[gooddates])
        qv = plt.quiver(llon, llat,
                   self.u[gooddates, depthindex] / self.velnorm[gooddates, depthindex],
                   self.v[gooddates, depthindex] / self.velnorm[gooddates, depthindex],
                   self.velnorm[gooddates, depthindex], headwidth=0, scale=25, cmap=cmocean.cm.speed)

        cb = plt.colorbar(qv, shrink=0.8, extend="max")
        cb.set_label("$\|v\|$\n(m/s)", rotation=0, ha="left", fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.clim(0, 1.)

        if depth:
            plt.title("Depth: {} m".format(depth), fontsize=20)

    def add_rectangle(self, N1, N2, m, dlon=0.02, dlat=0.02, label=None):
        """
        Draw a rectangle around the transect
        N1 and N2 are the indices of the extreme points
        of the considered section
        """
        lonmin = self.lon[N1:N2].min() - dlon
        lonmax = self.lon[N1:N2].max() + dlon
        latmin = self.lat[N1:N2].min() - dlat
        latmax = self.lat[N1:N2].max() + dlat
        lonrec = [lonmin, lonmax, lonmax, lonmin, lonmin]
        latrec = [latmin, latmin, latmax, latmax, latmin]
        m.plot(lonrec, latrec, "k-.", lw=1, latlon=True)
        # Add a label on top of the rectangle
        if label is not None:
            lontext = 0.5 * (lonmin + lonmax)
            lattext = latmax
            xt, yt = m(lontext, lattext)
            plt.text(xt, yt, label, fontsize=16, ha="center", va="bottom")

    @staticmethod
    def make_velocity_section(lat, depth, u, frontlat=None, title=None, xlabel=None):
        """
        Create a meridional section of zonal velocity
        Inputs:
        lat: 1-D array of latitudes
        depth: 1-D array of depths
        u: 2-D array of velocities
        """

        plt.pcolormesh(lat, depth, u, cmap=cmocean.cm.speed, vmin=0, vmax=1.)

        # Front position
        if frontlat is not None:
            plt.vlines(frontlat, 0, 400, colors='k', linestyles='--', linewidth=.5)

        if xlabel is not None:
            plt.xlabel(xlabel, fontsize=14)

        plt.ylabel("Depth\n(m)", rotation=0, ha="right", fontsize=14)
        cb = plt.colorbar(extend="max")
        cb.set_label("u\n(m/s)", rotation=0, ha="left", fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=12)

        if title is not None:
            plt.title(title, fontsize=20)
        xticks = np.arange(36.5, 37.5, 0.1)
        xticklabels = ["{}°N".format(np.round(xt,1)) for xt in xticks]
        plt.xticks(xticks, xticklabels)
        plt.xlim(lat.min(), lat.max())
        plt.ylim(0., 400.)
        plt.gca().invert_yaxis()

def prepare_3D_scat():

    fig = plt.figure(figsize=(12, 6))
    fig.patch.set_facecolor('white')

    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax1.set_aspect('equal')
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.set_aspect('equal')

    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")
    ax1.set_xticks(np.arange(-1., 0, 0.2))
    ax1.set_yticks(np.arange(36.8, 37.2, 0.2))
    ax1.set_title("Coastal glider")

    ax2.set_xlabel("Longitude")
    ax2.set_ylabel("Latitude")
    ax2.set_xticks(np.arange(-1., 0, 0.2))
    ax2.set_yticks(np.arange(36.8, 37.2, 0.2))
    ax2.set_title("Deep glider")
    fig.subplots_adjust(right=0.6)
    cbar_ax = fig.add_axes([0.65, 0.25, 0.015, 0.5])
    return fig, ax1, ax2, cbar_ax

def change_wall_prop(ax, coordinates, depths, angles):
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_xaxis.gridlines.set_linestyles(':')
    ax.w_yaxis.gridlines.set_linestyles(':')
    ax.w_zaxis.gridlines.set_linestyles(':')
    ax.view_init(angles[0], angles[1])
    ax.set_xlim(coordinates[0],coordinates[1])
    ax.set_ylim(coordinates[2],coordinates[3])
    ax.set_zlim(depths[0],depths[1])
    ax.set_zlabel('Depth (m)')

    ax.set_zticks(np.arange(depths[0],depths[1]+10,depths[2]))
    ax.set_zticklabels(range(int(-depths[0]),-int(depths[1])-10,-int(depths[2])))
