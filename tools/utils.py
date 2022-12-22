from netCDF4 import Dataset, date2num
from scipy.interpolate import griddata
import numpy as np

def interp(srcLons, srcLats, invar2d, dstLons, dstLats, value):
    py = srcLats.flatten()
    px = srcLons.flatten()
    z = np.array(invar2d).flatten()
    z[z == value] = 'nan'
    X, Y = np.meshgrid(dstLons, dstLats)
    outvar2d = griddata((px, py), z, (X, Y), method='nearest', fill_value=value)
    return outvar2d


def getIndexMaxMinLatLon(lats, lons, dtmFile):
    latli = np.argmin(np.abs(lats - np.min(dtmFile["latitude"])))
    latui = np.argmin(np.abs(lats - np.max(dtmFile["latitude"])))
    lonli = np.argmin(np.abs(lons - np.min(dtmFile["longitude"])))
    lonui = np.argmin(np.abs(lons - np.max(dtmFile["longitude"])))

    return latli, latui, lonli, lonui


def saveNetCDF(filename, lat, lon, time, hs, lm, dir, wspd10i, slp, value):
    ncdstfile = Dataset(filename, "w", format="NETCDF4")
    ncdstfile.createDimension("time", size=1)
    ncdstfile.createDimension("latitude", size=len(lat))
    ncdstfile.createDimension("longitude", size=len(lon))
    
    timeVar = ncdstfile.createVariable("time", "f4", "time")
    timeVar.description = "Time since initialization"
    timeVar.field = "time, scalar, series"
    timeVar.long_name = "julian day (UT)"
    timeVar.standard_name = "time"
    timeVar.calendar = "standard"
    timeVar.units = "days since 1990-01-01 00:00:00"
    timeVar.conventions = "relative julian days with decimal part (as parts of the day )"
    timeVar.axis = "T"

    lonVar = ncdstfile.createVariable("longitude", "f4", "longitude")
    lonVar.description = "Longitude"
    lonVar.long_name = "longitude"
    lonVar.units = "degrees_east"

    latVar = ncdstfile.createVariable("latitude", "f4", "latitude")
    latVar.description = "Latitude"
    latVar.long_name = "latitude"
    latVar.units = "degrees_north"

    if hs is not None and lm is not None and dir is not None:
        hsVar = ncdstfile.createVariable("hs", "f4", ("time", "latitude", "longitude"), fill_value=value)
        hsVar.description = "hs"
        hsVar.long_name = "significant height of wind and swell waves"
        hsVar.units = "m"
        
        lmVar = ncdstfile.createVariable("lm", "f4", ("time", "latitude", "longitude"), fill_value=value)
        lmVar.description = "lm"
        lmVar.long_name = "mean wave length"
        lmVar.units = "m"

        dirVar = ncdstfile.createVariable("dir", "f4", ("time", "latitude", "longitude"), fill_value=value)
        dirVar.description = "dir"
        dirVar.long_name = "wave mean direction"
        dirVar.units = "degree"

        timeVar[:] = time
        lonVar[:] = lon
        latVar[:] = lat
        hsVar[::] = hs
        lmVar[::] = lm
        dirVar[::] = dir
        
    elif wspd10i is not None and slp is not None:
        wspd10iVar = ncdstfile.createVariable("WSPD10", "f4", ("time", "latitude", "longitude"), fill_value=value)
        wspd10iVar.description = "WSPD10"
        wspd10iVar.long_name = "wind speed at 10 meters"
        wspd10iVar.units = "m s-1"
        
        slpVar = ncdstfile.createVariable("SLP", "f4", ("time", "latitude", "longitude"), fill_value=value)
        slpVar.description = "SLP"
        slpVar.long_name = "sea level pressure"
        slpVar.units = "HPa"

        time = date2num(time, units=timeVar.units)
        timeVar[:] = time
        lonVar[:] = lon
        latVar[:] = lat
        wspd10iVar[::] = wspd10i
        slpVar[::] = slp

    ncdstfile.close()