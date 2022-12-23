from cython.parallel import prange, threadid
from scipy.interpolate import griddata
from tools import utils
import numpy as np
import time
import math
import json

cimport numpy as np

G = 9.81
d = [0.667, 0.356, 0.161, 0.063, 0.021, 0.0066]

def bilinear_interpolation_2d(values, lons, lats, point):
    if values.ndim > 2:
        slc = [0] * (values.ndim - 2)
        slc += [slice(None), slice(None)]
        values = values[slc]

    min_lat, max_lat, min_lon, max_lon = np.min(lats), np.max(lats), np.min(lons), np.max(lons)

    lon = point[0]
    lat = point[1]
    
    i = len(lons)*(lon-min_lon)/(max_lon-min_lon)
    j = len(lats)*(lat-min_lat)/(max_lat-min_lat)

    jI, iI = int(j), int(i)
    jF, iF = j-jI, i - iI
    
    try:
        v1 = values[jI, iI] * (1.0 - iF) * (1.0-jF)
        v2 = values[jI+1, iI] * (1.0-iF) * jF
        v3 = values[jI+1, iI+1] * iF * jF
        v4 = values[jI, iI+1] * iF * (1.0-jF)

        return v1+v2+v3+v4
    except Exception as e:
        return float('NaN')


def calculateCOV(tmn, depth):
    sigma = (2*math.pi) / tmn
    parentesiOV = math.pow(sigma, 2) * (depth / G)
    sommaOV = 0
    for idx, j in enumerate(d):
        sommaOV += j * math.pow(parentesiOV, idx+1)
    denomOV = 1.0 + sommaOV
    frazOV = parentesiOV/denomOV
    termdexOV = frazOV + math.pow(parentesiOV, 2)
    cappaOV = math.sqrt(termdexOV) / depth
    lOV = (2*math.pi) * cappaOV
    cOV = lOV / tmn

    return cOV


def sind(d):
    return math.sin(math.radians(d))


def asind(a):
    return math.degrees(math.asin(a))


def cosd(d):
    return math.cos(math.radians(d))


def acosd(a):
    return math.degrees(math.acos(a))


def calcBreak(hs0, t0, gamma0):
    data = {}

    if gamma0 > 89 and gamma0 < -89:
        print("Error gamma0: ", gamma0)
        return data

    GammaBreak = 1000
    HsBreak = 0
    hBreak = 0
    LBreak = 0
    CBreak = 0

    C0 = G/(2*math.pi)*t0
    Cg0 = C0*0.5
    L0 = C0*t0
    Sin0 = sind(gamma0)
    Cos0 = cosd(gamma0)

    HBreak = hs0
    HBVecchio = HBreak

    NONCONV = True

    for kiteraz in range(0, 100):
        hBreak = HBreak*1.2
        CBreak = math.sqrt(G*hBreak)
        SinBreak = Sin0*CBreak/C0
        GammaBreak = asind(SinBreak)
        CosBreak = cosd(GammaBreak)
        Kr = math.sqrt(Cos0/CosBreak)
        Ks = math.sqrt(Cg0/CBreak)
        HBreak=hs0*Kr*Ks

        if abs((HBreak-HBVecchio)/HBVecchio) < 0.0001:
            NONCONV = False
            break

        HBVecchio = HBreak

    if NONCONV:
        print("Error NONCONV")
        return data

    HsBreak = HBreak
    LBreak = CBreak*t0

    data["gammaBreak"] = GammaBreak
    data["hsBreak"] = HsBreak
    data["hBreak"] = hBreak
    data["cBreak"] = CBreak

    return data


def transfBreak(data):
    hsOV = data["hs"]
    tmOV = data["tmn"]
    profOV = data["depth"]
    profTrans = -1*data["transect"][-1]["h"]
    gammaOV = data["gammaOV"]
    
    sigma = 2*math.pi/tmOV
    parentesi = math.pow(sigma, 2)*profTrans/G
    somma = 0
    for idx, j in enumerate(d):
        somma += j * math.pow(parentesi, idx+1)
    denom = 1.0 + somma
    fraz = parentesi/denom
    termdex = fraz + math.pow(parentesi, 2)
    cappaTrans = math.sqrt(termdex)/profTrans
    lTrans = 2*math.pi/cappaTrans
    cTrans = lTrans/tmOV

    parentesi = math.pow(sigma, 2)*profOV/G
    somma = 0
    for idx, j in enumerate(d):
        somma += j * math.pow(parentesi, idx+1)
    denom = 1.0 + somma
    fraz = parentesi/denom
    termdex = fraz + math.pow(parentesi, 2)
    cappaOV = math.sqrt(termdex)/profOV
    lOV = 2*math.pi/cappaOV  
    cOV=lOV/tmOV

    l0 = (G/2*math.pi)*math.pow(tmOV, 2)
    c0 = l0/tmOV
    
    hs0 = hsOV*cOV/c0
    gamma0 = gammaOV

    ksCosH = math.cosh(cappaTrans*profTrans)
    ksSinH = math.sinh(2*cappaTrans*profTrans)
    ks = math.sqrt(2*math.pow(ksCosH, 2)/(2*cappaTrans*profTrans + ksSinH))

    sGammaTrans = math.sin(gammaOV*math.pi/90)*cTrans/cOV
    gammaTrans = math.asin(sGammaTrans)*90/math.pi
    krCos1 = math.cos(gammaOV*math.pi/90)
    krCos2 = math.cos(gammaTrans*math.pi/90)
    kr = math.sqrt(krCos1/krCos2)

    hsTrans = hsOV*ks*kr
    hsTransProiett = hsTrans*math.cos(gammaTrans*math.pi/90)

    t0 = tmOV
 
    return calcBreak(hs0, t0, gamma0)


def runUp(data):
    l0 = G*math.pow(data["tmn"], 2)/(2*math.pi)
    c0 = l0/data["tmn"]
    hs0 = data["hs"]*data["cOV"]/c0
    a = 0.88
    b = 0.69

    c = data["sSlope"]/math.sqrt(hs0/l0)
    ru = a*math.pow(data["sSlope"]/math.sqrt(hs0/l0), b)*hs0
    return ru


def calculateIdx(wavemeter, lons, lats, lm, hs, dir, t0m1, depth):
    transect = json.loads(wavemeter[3])
    #m = list(filter(lambda t: t['d'] == 0, transect))[0]["m"]
    idx = next((index for (index, d) in enumerate(transect) if d["d"] == 0), None)
    m = transect[idx]["m"]
    
    earthPoint = wavemeter[-1].coords[0]
    ondameterPoint = wavemeter[-1].coords[1]
    
    hs = bilinear_interpolation_2d(hs, lons, lats, ondameterPoint)
    tmn = bilinear_interpolation_2d(t0m1, lons, lats, ondameterPoint)
    depth = bilinear_interpolation_2d(depth, lons, lats, ondameterPoint)

    dirmn = bilinear_interpolation_2d(dir, lons, lats, ondameterPoint)
    BetaTrans = wavemeter[1]
    GammaOV = dirmn - BetaTrans
    cOV = calculateCOV(tmn, abs(depth))

    kt = 0
    k = 0
    Hst = 0
    Hsl = 0
    ruback = 0

    dataItem = {
        "hs": hs,
        "tmn": tmn,
        "gammaOV": GammaOV,
        "cOV": cOV,
        "transect": transect,
        "sSlope": m
    }
    
    if -80 <= GammaOV <= 80:
        iDepth = -1.2 * hs

        try:
            for p in transect[idx+1:][::-1]:
                dataItem["depth"] = abs(p["h"])
                br = transfBreak(dataItem)
                dataItem["hs"] = br["hsBreak"]
                iDepth = -br["hBreak"]

                ru = runUp(dataItem)

            return 1
        except Exception as e:
            print("Error ", e)
            return 0
    else:
        return 1

cpdef np.ndarray[np.int_t,ndim=1] alert(np.ndarray wavemeters, np.ndarray lons, np.ndarray lats, np.ndarray lm, np.ndarray hs, np.ndarray dir, np.ndarray t0m1, np.ndarray depth):
    cdef np.ndarray[np.int_t,ndim=1]alerts_idx
    cdef int i, num = len(wavemeters)
    alerts_idx = np.zeros(num, dtype = int)

    for i in prange(num, nogil=True):
        with gil:
            alerts_idx[i] = calculateIdx(wavemeters[i], lons, lats, lm, hs, dir, t0m1, depth)

    return alerts_idx
