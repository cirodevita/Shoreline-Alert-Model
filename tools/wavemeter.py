import geopandas as gpd
import numpy as np
import pandas as pd
import sys
import geopy.distance

from math import sin, cos, radians, degrees, atan2, pi, asin, isnan, ceil
from geojson import LineString, Feature, FeatureCollection, dump
from netCDF4 import Dataset
from shapely.geometry import Point, Polygon


# Configurations
points_distance = 25 #m
max_distance = 12.5 #km
step_size_profile = 0.0125 #Km
max_iterations = max_distance/step_size_profile
max_L = 0.1 #Km
max_H = 0.010 #Km


def calculate_bearing(prev, curr):
    lon_prev, lat_prev = prev[0], prev[1]
    lon_curr, lat_curr = curr[0], curr[1]
    y = sin(radians(lon_curr) - radians(lon_prev)) * cos(radians(lat_curr))
    x = cos(radians(lat_prev)) * sin(radians(lat_curr)) - sin(radians(lat_prev)) * cos(radians(lat_curr)) * cos(radians(lat_curr) - radians(lat_prev))
    theta = atan2(y, x)
    bearing = (theta * 180 / pi + 360) % 360
    
    return bearing


def calculate_new_position(point, distance, angle):
    R = 6378.1
    lon = point[0]
    lat = point[1]
    
    new_lat = degrees(asin(sin(radians(lat)) * cos(distance/R) + cos(radians(lat)) * sin(distance/R) * cos(radians(angle))))
    new_lon = degrees(radians(lon) + atan2(sin(radians(angle)) * sin(distance/R) * cos(radians(lat)), cos(distance/R) - sin(radians(lat)) * sin(radians(new_lat))))
    
    return Point(new_lon, new_lat)


def bilinear_interpolation_2d(netcdf, variable, point):
    values = netcdf[variable]
    if values.ndim > 2:
        slc = [0] * (values.ndim - 2)
        slc += [slice(None), slice(None)]
        values = values[slc]

    lats = np.asarray(netcdf["latitude"])
    lons = np.asarray(netcdf["longitude"])
    min_lat, max_lat, min_lon, max_lon = np.min(lats[:]), np.max(lats[:]), np.min(lons[:]), np.max(lons[:])

    lon = point.x
    lat = point.y
    i = len(lons)*(lon-min_lon)/(max_lon-min_lon)
    j = len(lats)*(lat-min_lat)/(max_lat-min_lat)

    jI, iI = int(j), int(i)
    jF, iF = j-jI, i - iI
    
    try:
        v1 = values[jI, iI] * (1.0 - iF) * (1.0-jF)
        v2=values[jI+1, iI] * (1.0-iF) * jF
        v3=values[jI+1, iI+1] * iF * jF
        v4=values[jI, iI+1] * iF * (1.0-jF)

        return v1+v2+v3+v4
    except Exception as e:
        return float('NaN')


def calculate_sea_profile(point, angle, ww3, dtm, step_size, max_iterations):
    lm_ww3 = float("nan")
    step = 0
    iter = 0
    profile = []
    positions = []

    while (isnan(lm_ww3) and iter != max_iterations):
        step = float("{:.4f}".format(step + step_size))
        new_position = calculate_new_position(point, step, angle)
        lm_ww3 = bilinear_interpolation_2d(ww3, "lm", new_position)
        depth = bilinear_interpolation_2d(dtm, "height", new_position)
        profile.append({
            "d": step*1000,
            "h": depth
        })
        positions.append(new_position)
        iter += 1

    return profile, positions


def calculate_earth_profile(point, angle, dtm, step_size, max_L, max_H, max_iterations):
    step = step_size
    depth = 0
    iter = 0
    profile = []
    positions = []

    while (depth < max_H and step < max_L and iter != max_iterations):
        new_position = calculate_new_position(point, step, angle)
        depth = bilinear_interpolation_2d(dtm, "height", new_position)
        profile.append({
            "d": -step*1000,
            "h": depth
        })
        positions.append(new_position)
        step = float("{:.4f}".format(step + step_size))
        iter += 1

    return profile, positions


"""
def PosThenNeg(L):
    L = list(filter(None,L))
    return any(a>0 and b<0 for a,b in zip(L,L[1:]))

def negThenPos(L):
    L = list(filter(None,L))
    return any(a<0 and b>0 for a,b in zip(L,L[1:]))

def getPosNegIdx(L):
    L = list(filter(None,L))
    return [[idx, idx+1] for idx, [a, b] in enumerate(zip(L, L[1:])) if a>0 and b<0]
"""

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python " + str(sys.argv[0]) + " shapefile.shp dtmFile.nc regridded_WW3_file.nc output_directory")
        sys.exit(-1)

    lines = gpd.read_file(sys.argv[1])
    dtm = Dataset(sys.argv[2])
    ww3 = Dataset(sys.argv[3])
    outputDir = sys.argv[4]

    min_lat, max_lat, min_lon, max_lon = np.min(dtm["latitude"][:]), np.max(dtm["latitude"][:]), np.min(dtm["longitude"][:]), np.max(dtm["longitude"][:])
    area = Polygon([[min_lon, min_lat], [min_lon, max_lat], [max_lon, max_lat], [max_lon, min_lat]])

    wavemeters = []
    id = 0

    for index, line in lines.iterrows():
        if line["geometry"] is not None:
            points = line["geometry"].coords[:-1]

            for i in range(1, len(points)):
                print("{}/{}".format(i, len(points)))

                i0 = points[i-1]
                i1 = points[i]

                if Point(i0).within(area) and Point(i1).within(area):
                    distance = geopy.distance.geodesic(i0, i1).m
                    points_section = [i0, i1]

                    if distance > points_distance:
                        lons = np.linspace(i0[0], i1[0], ceil(distance/points_distance), endpoint=True)
                        lats = np.linspace(i0[1], i1[1], ceil(distance/points_distance), endpoint=True)
                        points_section = list(zip(lons, lats))

                    for j in range(1, len(points_section)):
                        prev = points_section[j-1]
                        curr = points_section[j]

                        # Calculate bearing
                        bearing = calculate_bearing(prev, curr)
                        sea_direction = (bearing + 90) % 360
                        earth_direction = (bearing - 90) % 360

                        sea_profile, sea_positions = calculate_sea_profile(prev, sea_direction, ww3, dtm, step_size_profile, max_iterations)
                        h_sea = list(data['h'] for data in sea_profile)

                        earth_profile, earth_positions = calculate_earth_profile(prev, earth_direction, dtm, step_size_profile, max_L, max_H, max_iterations)
                        h_earth = list(data['h'] for data in earth_profile)
                        
                        if not np.isnan(h_sea).any() and len(h_sea) <= max_iterations/2 and not np.isnan(h_earth).any() and len(h_earth) <= max_iterations/2:
                            print("Added!")
                            sea_profile.insert(0, {
                                "d": 0,
                                "h": 0
                            })
                            earth_profile.reverse()
                            final_profile = earth_profile + sea_profile

                            for profile in final_profile:
                                angle = sea_direction
                                if profile["d"] < 0:
                                    angle = earth_direction
                                curr_point = calculate_new_position(prev, profile["d"]/1000, angle)

                                p0 = calculate_new_position(curr_point.coords[0], step_size_profile/2, sea_direction)
                                p1 = calculate_new_position(curr_point.coords[0], step_size_profile/2, earth_direction)
                                
                                d0 = bilinear_interpolation_2d(dtm, "height", p0)
                                d1 = bilinear_interpolation_2d(dtm, "height", p1)

                                h = d1 - d0
                                d = geopy.distance.geodesic(p1.coords[0], p0.coords[0]).m

                                m = atan2(h, d)
                                profile["m"] = m

                            line = LineString([list(earth_positions[0].coords[0]), list(sea_positions[-1].coords[0])])
                            properties = {
                                "id": id,
                                "angle_sea": sea_direction,
                                "angle_earth": earth_direction,
                                "Z": final_profile
                            }
                            wavemeters.append(Feature(geometry=line, properties=properties))
                            
                            id += 1

    wavemeter_collection = FeatureCollection(wavemeters)
    with open(outputDir + '/wavemeter.geojson', 'w') as f:
       dump(wavemeter_collection, f)

