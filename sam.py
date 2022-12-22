from mpi4py import MPI
from netCDF4 import Dataset
from itertools import islice
import geopandas as gpd
import numpy as np
import omp
import time


def alert(wavemeters, lons, lats, lm, hs, dir, t0m1, depth):
    cumulativeEffectClass = omp.alert(wavemeters.to_numpy(), lons, lats, lm, hs, dir, t0m1, depth)

    return cumulativeEffectClass


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
    start_time = time.time()

    #df = gpd.read_file("data/wavemeter.geojson").iloc[:1]
    df = gpd.read_file("data/wavemeter.geojson")
    totalSize = len(df)
    print(totalSize)
    wavemeters = [None] * size

    spare = totalSize % size
    data_part = [int(x) for x in np.linspace(spare, len(df), size+1)]
    for i in range(size):
        if i == 0:
            wavemeters[i] = df.iloc[data_part[i]-spare:data_part[i+1]]
        else:
            wavemeters[i] = df.iloc[data_part[i]:data_part[i+1]]

    ww3 = Dataset("data/regridded_ww33_d03_20221224Z2000.nc")
    lm = ww3["lm"][0][:]
    hs = ww3["hs"][0][:]
    dir = ww3["dir"][0][:]
    t0m1 = ww3["t0m1"][0][:]

    dtm = Dataset("data/DTM.nc")
    depth = dtm["height"][:]
    lons = dtm["longitude"][:]
    lats = dtm["latitude"][:]
else:
    totalSize = None
    wavemeters = None
    lm = None
    hs = None
    dir = None
    t0m1 = None
    depth = None
    lons = None
    lats = None

wavemeters = comm.scatter(wavemeters, root = 0)
lm = comm.bcast(lm, root=0)
hs = comm.bcast(hs, root=0)
dir = comm.bcast(dir, root=0)
t0m1 = comm.bcast(t0m1, root=0)
depth = comm.bcast(depth, root=0)
lons = comm.bcast(lons, root=0)
lats = comm.bcast(lats, root=0)
totalSize = comm.bcast(totalSize, root=0)

classVarPartition = alert(wavemeters, lons, lats, lm, hs, dir, t0m1, depth)

classVar = None
if rank == 0:
    classVar = np.zeros(len(wavemeters)*size, dtype = int)

comm.Gather(classVarPartition, classVar, root=0)

if rank == 0:
    end_time = time.time() - start_time
    print("Exectution time: ", end_time)
    print('Rank: ',rank, ', recvbuf received: ', classVar[:totalSize])
