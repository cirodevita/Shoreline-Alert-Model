import geopandas as gpd
import numpy as np
import geopy
import sys

from sqlalchemy import create_engine
from geopy.distance import geodesic
from scipy.interpolate import griddata
from netCDF4 import Dataset

# import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap


def calculate_step(lon, lat, meters):
    x1 = geodesic(meters=meters).destination(point=geopy.Point(lat, lon), bearing=90)
    x2 = geodesic(meters=meters).destination(point=geopy.Point(lat, lon), bearing=0)
    
    step_lon = x1.longitude - lon
    step_lat = x2.latitude - lat
    
    return step_lon, step_lat


def saveasNetCDF(filename, grid_lat, grid_lon, z):
    ncdstfile = Dataset(filename + ".nc", "w", format="NETCDF4")
    ncdstfile.createDimension("latitude", size=len(grid_lat))
    ncdstfile.createDimension("longitude", size=len(grid_lon))

    lonVar = ncdstfile.createVariable("longitude", "f4", "longitude")
    lonVar.description = "Longitude"
    lonVar.long_name = "longitude"
    lonVar.units = "degrees_east"

    latVar = ncdstfile.createVariable("latitude", "f4", "latitude")
    latVar.description = "Latitude"
    latVar.long_name = "latitude"
    latVar.units = "degrees_north"

    heightVar = ncdstfile.createVariable("height", "f4", ("latitude", "longitude"))
    heightVar.description = "Height"
    heightVar.long_name = "height"
    heightVar.units = "meters"

    lonVar[:] = grid_lon
    latVar[:] = grid_lat
    heightVar[::] = z

    ncdstfile.close()
    

x_min, x_max = 14.34, 14.45
y_min, y_max = 40.69, 40.80
step = 10 #meters

if len(sys.argv)!=2:
    print("Usage: python " + str(sys.argv[0]) + " netCDFilenameOutput")
    sys.exit(-1)

engine = create_engine(f'postgresql://dtm:Dtm2022@127.0.0.1:5432/dtm?gssencmode=disable')
sql = "SELECT * FROM public.points WHERE geom && ST_MakeEnvelope({}, {}, {}, {}, 4326)".format(x_min, y_min, x_max, y_max)
points = gpd.GeoDataFrame.from_postgis(sql, engine, geom_col="geom")

lons = np.array(points["geom"].x)
lats = np.array(points["geom"].y)
depth = np.array(points["geom"].z)

min_lat, max_lat, min_lon, max_lon = np.min(lats), np.max(lats), np.min(lons), np.max(lons)
step_lon, step_lat = calculate_step(min_lon, min_lat, step)

grid_lon = np.arange(np.min(lons), np.max(lons), step_lon)
grid_lat = np.arange(np.min(lats), np.max(lats), step_lat)

xintrp, yintrp = np.meshgrid(grid_lon, grid_lat)

print("Start Interpolation...")
z = griddata((lons, lats), depth, (xintrp, yintrp), method='nearest')

print("Saving as NetCDF...")

saveasNetCDF(sys.argv[1], grid_lat, grid_lon, z)

print("Finish")

# print("Creating Plotting...")
# fig, ax = plt.subplots(figsize=(10, 10))
# m = Basemap(llcrnrlon=lons.min()-0.1,llcrnrlat=lats.min()-0.1,urcrnrlon=lons.max()+0.1,urcrnrlat=lats.max()+0.1, projection='merc', resolution='h',area_thresh=1000.,ax=ax)
# m.drawcoastlines()
# x, y = m(xintrp, yintrp)
# ln, lt = m(lons, lats)

# contour = ax.contourf(xintrp, yintrp, z, np.linspace(np.min(depth), np.max(depth), 4000, endpoint=True), cmap='jet', zorder=1)
# cb = plt.colorbar(contour, ax=ax, ticks=np.linspace(np.min(depth), np.max(depth), 30, endpoint=True))
# plt.show()
# plt.savefig(sys.argv[1] + '.png')
