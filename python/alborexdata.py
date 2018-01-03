import os
import netCDF4
import logging
import numpy as np


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

class Drifter(object):

    def __ini__(self, lon=None, lat=None, time=None, temperature=None, qclon=None, qclat=None):
        self.lon = lon
        self.lat = lat
        self.time = time
        self.temperature = temperature
        self.qclon = qclon
        self.qclat = qclat
        self.timeunits = None

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
        don't have the indicated quality flag
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
        Mask the time outside the mission period
        """
        dates = netCDF4.num2date(self.time, self.timeunits)
        if initialdate is not None:
            self.lon = np.ma.masked_where(np.logical_or(dates > finaldate, dates < initialdate))
            self.lat = np.ma.masked_where(np.logical_or(dates > finaldate, dates < initialdate))
        else:
            self.lon = np.ma.masked_where(dates > finaldate, self.lon)
            self.lat = np.ma.masked_where(dates > finaldate, self.lat)

    def scatter_plot(self, m, **kwargs):
        scat = m.scatter(self.lon, self.lat, c=self.temperature, latlon=True, **kwargs)
        return scat

    def point_plot(self, m, **kwargs):
        m.plot(self.lon, self.lat, latlon=True, linewidth=0, **kwargs)

    def add_initial_position(self, m, **kwargs):
        m.plot(self.lon[0], self.lat[0], latlon=True, linewidth=0, **kwargs)
