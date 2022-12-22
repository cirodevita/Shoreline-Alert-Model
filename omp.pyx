from cython.parallel import prange, threadid
from scipy.interpolate import griddata
from tools import utils
import numpy as np
import time

cimport numpy as np


class RangeDict(dict):
    def __getitem__(self, item):
        if not isinstance(item, range):
            for key in self:
                if item in key:
                    return self[key]
            raise KeyError(item)
        else:
            return super().__getitem__(item)


windSpeedRange = RangeDict({
    range(0, 108): 1,
    range(108, 131): 2,
    range(131, 155): 3,
    range(155, 178): 4,
    range(178, 1000): 5
})

waveHeightRange = RangeDict({
    range(0, 125): 1,
    range(125, 250): 2,
    range(250, 350): 3,
    range(350, 450): 4,
    range(450, 1000): 5
})

seaLevelPressureRange = RangeDict({
    range(1013, 2000): 1,
    range(999, 1013): 2,
    range(998, 999): 3,
    range(997, 998): 4,
    range(0, 997): 5
})


def simpleMethod(wind_speed, sea_pressure, wave_height):
    windSpeedClass = windSpeedRange[round(wind_speed) * 10]
    waveHeightClass = waveHeightRange[round(wave_height) * 100]
    seaLevelPressureClass = seaLevelPressureRange[round(sea_pressure)]

    cumulativeEffectClass = round((windSpeedClass + waveHeightClass + seaLevelPressureClass) / 3)

    return cumulativeEffectClass


def calculateIdx(wavemeter, lm, hs, dir, slp, wspd10, depth):
    alpha = wavemeter[4]

    points = wavemeter[-1].coords[:]
    point0, point1 = points[0], points[1]

    wind_speed = wspd10[wavemeter[0]][wavemeter[1]]
    sea_pressure = slp[wavemeter[0]][wavemeter[1]]
    wave_height = hs[wavemeter[2]][wavemeter[3]]
    wave_length = lm[wavemeter[2]][wavemeter[3]]
    # t0m1
    
    beta = dir[wavemeter[2]][wavemeter[3]]
    alpha = (alpha + 180) % 360
    angle = beta - alpha
    if -80 <= angle <= 80:
        time.sleep(0.01)

        return simpleMethod(wind_speed, sea_pressure, wave_height)
    else:
        return 1


cpdef np.ndarray[np.int_t,ndim=1] alert(np.ndarray wavemeters, np.ndarray lm, np.ndarray hs, np.ndarray dir, np.ndarray slp, np.ndarray wspd10, np.ndarray depth):
    cdef np.ndarray[np.int_t,ndim=1]alerts_idx
    cdef int i, num = len(wavemeters)
    alerts_idx = np.zeros(num, dtype = int)

    for i in prange(num, nogil=True):
        with gil:
            alerts_idx[i] = calculateIdx(wavemeters[i], lm, hs, dir, slp, wspd10, depth)

    return alerts_idx
