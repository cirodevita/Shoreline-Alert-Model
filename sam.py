from mpi4py import MPI
from netCDF4 import Dataset
from itertools import islice
import geopandas as gpd
import numpy as np
import omp
import time


def alert(wavemeters, lm, hs, dir, slp, wspd10, depth):
    cumulativeEffectClass = omp.alert(wavemeters.to_numpy(), lm, hs, dir, slp, wspd10, depth)

    return cumulativeEffectClass


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
    start_time = time.time()

    df = gpd.read_file("data/virtualWavemeters.json").iloc[:1]
    totalSize = len(df)
    wavemeters = [None] * size

    spare = totalSize % size
    data_part = [int(x) for x in np.linspace(spare, len(df), size+1)]
    for i in range(size):
        if i == 0:
            wavemeters[i] = df.iloc[data_part[i]-spare:data_part[i+1]]
        else:
            wavemeters[i] = df.iloc[data_part[i]:data_part[i+1]]

    ww3 = Dataset("ww33/regridded/ww33_d03_20221201Z0300.nc.nc")
    lm = ww3["lm"][0][:]
    hs = ww3["hs"][0][:]
    dir = ww3["dir"][0][:]

    wrf5 = Dataset("wrf5/regridded/wrf5_d03_20221201Z0300.nc.nc")
    slp = wrf5['SLP'][0][:]
    wspd10 = wrf5['WSPD10'][0][:]

    dtm = Dataset("data/sam_d03_terrain.nc")
    depth = dtm["height"][:]
else:
    wavemeters = None
    lm = None
    hs = None
    dir = None
    slp = None
    wspd10 = None
    totalSize = None
    depth = None

wavemeters = comm.scatter(wavemeters, root = 0)
lm = comm.bcast(lm, root=0)
hs = comm.bcast(hs, root=0)
dir = comm.bcast(dir, root=0)
slp = comm.bcast(slp, root=0)
wspd10 = comm.bcast(wspd10, root=0)
totalSize = comm.bcast(totalSize, root=0)

classVarPartition = alert(wavemeters, lm, hs, dir, slp, wspd10, depth)

classVar = None
if rank == 0:
    classVar = np.zeros(len(wavemeters)*size, dtype = int)

comm.Gather(classVarPartition, classVar, root=0)

if rank == 0:
    end_time = time.time() - start_time
    print("Exectution time: ", end_time)
    print('Rank: ',rank, ', recvbuf received: ', classVar[:totalSize])