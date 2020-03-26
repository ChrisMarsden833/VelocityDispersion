import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import ctypes
from colossus.cosmology import cosmology
from scipy import interpolate
from scipy.interpolate import griddata
import time

cosmo = cosmology.setCosmology('planck18')
omega_m = cosmo.Om0
H = cosmo.H0
ibc = ctypes.CDLL("../cmake-build-debug/lib/libfoo.so")

def FullVelocityDispersion(ApertureSize, Beta, HalfLightRadius, SersicIndex, StellarMass, HaloRs, HaloC, OmegaM, H):
    """I am a docstring"""
    # Test that the supplied lengths are consistent
    variables = [ApertureSize, Beta, HalfLightRadius, SersicIndex, StellarMass, HaloRs, HaloC]
    variable_names = ["Aperture Size", "Beta", "Half Light Radius", "Sersic Index", "Stellar Mass", "Halo Rs", "Halo c"]
    for i, element in enumerate(variables):
        if i != 0:
            assert len(element) == previous_length, "Supplied Array lengths are inconsistent - {} has {} elements, wheras {} has {} elements".format(variable_names[i-1], previous_length, varaible_names[i], len(element)) 
        previous_length = len(element)

    # Reserve Memory for ctype arrays
    c_aperture = (ctypes.c_float * len(ApertureSize))()
    c_beta = (ctypes.c_float * len(Beta))()
    c_hlr = (ctypes.c_float * len(HalfLightRadius))()
    c_n = (ctypes.c_float * len(SersicIndex))()
    c_sm = (ctypes.c_float * len(StellarMass))()
    c_rs = (ctypes.c_float * len(HaloRs))()
    c_hc = (ctypes.c_float * len(HaloC))()

    for i in range(len(ApertureSize)):
        c_aperture[i] = ApertureSize[i]
        c_beta[i] = Beta[i]
        c_hlr[i] = HalfLightRadius[i]
        c_n[i] = SersicIndex[i]
        c_sm[i] = StellarMass[i]
        c_rs[i] = HaloRs[i]
        c_hc[i] = HaloC[i]

    c_om = float(OmegaM)
    c_H = float(H)
    c_size = int(len(StellarMass))
   
    ibc.ParallelSigma.argtypes = [ ctypes.POINTER(ctypes.c_float), 
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float), 
                               ctypes.c_float,
                               ctypes.c_float,
                               ctypes.c_int32 ]

    ibc.ParallelSigma.restype = ctypes.POINTER(ctypes.c_float)

    res = ibc.ParallelSigma(c_aperture, c_beta, c_hlr, c_n, c_sm, c_rs, c_hc, c_om, c_H, c_size)

    output = np.ones_like(StellarMass)

    for i  in range(len(output)):
        output[i] = res[i]
    
    return output 


if __name__ == "__main__":

    length = 1000

    Aperture = np.ones(length) * 0.1
    Beta = np.random.normal(0.15, 0.12, size=length)
    HLR = np.ones(length) * 10**0.9
    n = np.ones(length) * 3
    SM = np.ones(length) * 10**11.3
    DMr = HLR * 20
    DMc = np.ones(length) * 1.3
    
    om = 0.27
    H = 67.66

    t1_start = time.time()

    FullVelocityDispersion(Aperture, Beta, HLR, n, SM, DMr, DMc, om, H)

    t1_stop = time.time() 
   
    print("Time:", t1_stop-t1_start)


