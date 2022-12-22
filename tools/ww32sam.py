import sys
import numpy as np
import utils
from netCDF4 import Dataset
from scipy.interpolate import griddata
from wrf import getvar

if len(sys.argv) != 4:
    print("Usage: python " + str(sys.argv[0]) + " DTM_file.nc src_his_WW3_file.nc output_directory_path")
    sys.exit(-1)

# ------------- DTM FILE -------------
dtmfile = Dataset(sys.argv[1])
grid_lat = np.linspace(np.min(dtmfile["latitude"]), np.max(dtmfile["latitude"]), len(dtmfile["latitude"]))
grid_lon = np.linspace(np.min(dtmfile["longitude"]), np.max(dtmfile["longitude"]), len(dtmfile["longitude"]))
X, Y = np.meshgrid(grid_lon, grid_lat)

# ------------- WW3 FILE -------------
ww3file = Dataset(sys.argv[2])
fill_value = 9.9692100e+36

latsww3, lonsww3 = ww3file["latitude"], ww3file["longitude"]
latliww3, latuiww3, lonliww3, lonuiww3 = utils.getIndexMaxMinLatLon(latsww3, lonsww3, dtmfile)
lonsww3 = lonsww3[lonliww3:lonuiww3]
latsww3 = latsww3[latliww3:latuiww3]

lm = ww3file["lm"][:, latliww3:latuiww3, lonliww3:lonuiww3]
hs = ww3file["hs"][:, latliww3:latuiww3, lonliww3:lonuiww3]
dir = ww3file["dir"][:, latliww3:latuiww3, lonliww3:lonuiww3]
t0m1 = ww3file["t0m1"][:, latliww3:latuiww3, lonliww3:lonuiww3]

X_ww3, Y_ww3 = np.meshgrid(lonsww3, latsww3)
hs = utils.interp(X_ww3, Y_ww3, hs, grid_lon, grid_lat, fill_value)
lm = utils.interp(X_ww3, Y_ww3, lm, grid_lon, grid_lat, fill_value)
dir = utils.interp(X_ww3, Y_ww3, dir, grid_lon, grid_lat, fill_value)
t0m1 = utils.interp(X_ww3, Y_ww3, t0m1, grid_lon, grid_lat, fill_value)
time = ww3file["time"][:]

utils.saveNetCDF(sys.argv[3] + "/" + sys.argv[2].split("/")[-1], grid_lat, grid_lon, time, hs, lm, dir, t0m1, None, None, fill_value)
