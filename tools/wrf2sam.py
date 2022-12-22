import sys
import math
import numpy as np
import utils
from netCDF4 import Dataset
from wrf import getvar
from datetime import datetime


def getBoundaries(Xlon, Xlat):
    row_lat = len(Xlat) - 1
    col_lat = len(Xlat[0]) - 1

    row_long = len(Xlon) - 1
    col_long = len(Xlon[0]) - 1

    A = [Xlat[0][0], Xlon[0][0]]
    B = [Xlat[0][col_lat], Xlon[0][col_long]]
    C = [Xlat[row_lat][col_lat], Xlon[row_long][col_long]]
    D = [Xlat[row_lat][0], Xlon[row_long][0]]

    min_lat = Xlat[0][0]
    minI = 0

    ''' from A to B '''
    for i in range(col_lat,-1,-1):
        np1 = [Xlat[0][i], Xlon[0][i]]
        if np1[0] > min_lat:
            minI = i
            min_lat = np1[0]

    max_lat = Xlat[row_lat][col_lat]
    maxI = col_lat

    ''' from C to D '''
    for i in range(col_lat,-1,-1):
        np1 = [Xlat[row_lat][i], Xlon[row_long][i]]
        if np1[0] < max_lat:
            maxI = i
            max_lat = np1[0]

    min_long = Xlon[0][0]
    minJ = 0

    ''' from A to D '''
    for i in range(row_lat,-1,-1):
        np1 = [Xlat[i][0], Xlon[i][0]]
        if np1[1] > min_long:
            minJ = i
            min_long = np1[1]

    max_long = Xlon[0][col_long]
    maxJ = row_lat

    ''' from B to C '''
    for i in range(row_lat,-1,-1):
        np1 = [Xlat[i][col_lat], Xlon[i][col_lat]]
        if np1[1] < max_long:
            maxJ = i
            max_long = np1[1]
            
    minLat=np.asscalar(min_lat)
    maxLat=np.asscalar(max_lat)
    minLon=np.asscalar(min_long)
    maxLon=np.asscalar(max_long)

    return minLon,minLat,maxLon,maxLat
    

def getLatsLonsWrf5(dn, de, Xlat, Xlon):
    #Earth's radius, sphere
    R=6378137

    lon = np.average(Xlon)
    lat = np.average(Xlat)

    #Coordinate offsets in degrees
    dLat = 0.5*(dn/R)*180/math.pi
    dLon = 0.5*(de/(R*math.cos(math.pi*lat/180)))*180/math.pi

    minLon, minLat, maxLon, maxLat = getBoundaries(Xlon, Xlat)
    # Create the latitude array
    lats =  np.arange(minLat, maxLat, dLat)
    # create the longitude array
    lons =  np.arange(minLon, maxLon, dLon)
    
    return lats, lons
    

if len(sys.argv) != 4:
    print("Usage: python " + str(sys.argv[0]) + " DTM_file.nc src_his_WRF_file.nc output_directory_path")
    sys.exit(-1)

# ------------- DTM FILE -------------
dtmfile = Dataset(sys.argv[1])
grid_lat = np.linspace(np.min(dtmfile["latitude"]), np.max(dtmfile["latitude"]), len(dtmfile["latitude"]))
grid_lon = np.linspace(np.min(dtmfile["longitude"]), np.max(dtmfile["longitude"]), len(dtmfile["longitude"]))
X, Y = np.meshgrid(grid_lon, grid_lat)

# ------------- WRF5 FILE -------------
wrf5file = Dataset(sys.argv[2])
fill_value = 1.e+37

Xlon = np.array(getvar(wrf5file, "XLONG", meta=False))
Xlat = np.array(getvar(wrf5file, "XLAT", meta=False))
latswrf5, lonswrf5 = getLatsLonsWrf5(wrf5file.DY, wrf5file.DX, Xlat, Xlon)

uvmet10_wspd_wdir = getvar(wrf5file, "uvmet10_wspd_wdir", meta=False)
wspd10i = utils.interp(Xlon, Xlat, uvmet10_wspd_wdir[0], lonswrf5, latswrf5, fill_value)

slp = getvar(wrf5file, "slp", meta=False)
slp = utils.interp(Xlon, Xlat, slp, lonswrf5, latswrf5, fill_value)

latliwrf5, latuiwrf5, lonliwrf5, lonuiwrf5 = utils.getIndexMaxMinLatLon(latswrf5, lonswrf5, dtmfile)
lonswrf5 = lonswrf5[lonliwrf5:lonuiwrf5]
latswrf5 = latswrf5[latliwrf5:latuiwrf5]

wspd10i = wspd10i[latliwrf5:latuiwrf5, lonliwrf5:lonuiwrf5]
slp = slp[latliwrf5:latuiwrf5, lonliwrf5:lonuiwrf5]

X_wrf5, Y_wrf5 = np.meshgrid(lonswrf5, latswrf5)
wspd10i = utils.interp(X_wrf5, Y_wrf5, wspd10i, grid_lon, grid_lat, fill_value)
slp = utils.interp(X_wrf5, Y_wrf5, slp, grid_lon, grid_lat, fill_value)

datetimeStr = bytearray(wrf5file.variables["Times"][:][0]).decode('ascii')
dateStr = datetimeStr.split("_")[0].split("-")
timeStr = datetimeStr.split("_")[1].split(":")
datetime_current = datetime(int(dateStr[0]), int(dateStr[1]), int(dateStr[2]), int(timeStr[0]))
time = [ datetime_current ]

utils.saveNetCDF(sys.argv[3] + "/" + sys.argv[2].split("/")[-1], grid_lat, grid_lon, time, None, None, None, wspd10i, slp, fill_value)

