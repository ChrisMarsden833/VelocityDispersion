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
ibc = ctypes.CDLL("/Users/chris/Documents/PhD/ProjectSigma/VelocityDispersion/lib/libsigma.so")

def FullVelocityDispersion(ApertureSize, Beta, HalfLightRadius, SersicIndex,
        StellarMass, z, DM = None, HaloMass = None, cpath = "../data/cM_planck18.txt" ):
    """I am a docstring"""
    # Test that the supplied lengths are consistent
    variables = [ApertureSize, Beta, HalfLightRadius, SersicIndex, StellarMass]
    variable_names = ["Aperture Size", "Beta", "Half Light Radius", "Sersic Index", "Stellar Mass"]

    for i, element in enumerate(variables):
        if i != 0:
            assert len(element) == previous_length, """Supplied Array lengths are
                inconsistent - {} has {} elements, wheras {} has {}
                elements""".format(variable_names[i-1], previous_length,
                variable_names[i], len(element))
        previous_length = len(element)
        length = len(StellarMass)

    assert np.sum(ApertureSize <= 0) == 0, "ApertudeSize has elements < 0, {}".format(ApertureSize[ApertureSize <= 0])
    assert np.sum(HalfLightRadius <= 0) == 0, "HalfLightRadius has elements < 0, {}".format(HalfLightRadius[HalfLightRadius <= 0])
    assert np.sum(SersicIndex <= 0) == 0, "Sersic Index has elements < 0 {}".format(SersicIndex[SersicIndex <=0])


    c_float_p = ctypes.POINTER(ctypes.c_float)

    ApertureSize = ApertureSize.astype(np.float32)
    c_aperture = ApertureSize.ctypes.data_as(c_float_p)

    Beta = Beta.astype(np.float32)
    c_beta = Beta.ctypes.data_as(c_float_p)

    HalfLightRadius = HalfLightRadius.astype(np.float32)
    c_hlr = HalfLightRadius.ctypes.data_as(c_float_p)

    SersicIndex = SersicIndex.astype(np.float32)
    c_n = SersicIndex.ctypes.data_as(c_float_p)

    StellarMass = StellarMass.astype(np.float32)
    c_sm = StellarMass.ctypes.data_as(c_float_p)

    if not hasattr(z, "__len__"):
        z = np.ones_like(StellarMass) * z

    if DM is None:
        assert HaloMass is None, "Halo mass should not be specified if DM is None"
        HaloMass = np.zeros_like(StellarMass)
        DM = "None"

    c_DM = DM.encode('utf-8')
    c_cpath = cpath.encode('utf-8')

    HaloMass = HaloMass.astype(np.float32)
    c_hm = HaloMass.ctypes.data_as(c_float_p)

    z = z.astype(np.float32)
    c_z = z.ctypes.data_as(c_float_p)


    #HaloRs = HaloRs.astype(np.float32)
    #c_rs = HaloRs.ctypes.data_as(c_float_p)
    #HaloRs = None

    #HaloC = HaloC.astype(np.float32)
    #c_hc = HaloC.ctypes.data_as(c_float_p)
    #HaloC = None

    #c_om = float(OmegaM)
    #c_H = float(H)
    c_size = int(length)

    ibc.ParallelSigma.argtypes = [ ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.POINTER(ctypes.c_float),
                               ctypes.c_int32,
                               ctypes.c_char_p,
                               ctypes.c_char_p]

                               #ctypes.POINTER(ctypes.c_float)]
                               #ctypes.POINTER(ctypes.c_float),
                               #ctypes.c_float,
                               #ctypes.c_float,
                               #ctypes.c_int32 ]

    ibc.ParallelSigma.restype = ctypes.POINTER(ctypes.c_float)

    res = ibc.ParallelSigma(c_aperture, c_beta, c_hlr, c_n, c_sm, c_hm, c_z,
            c_size, c_DM, c_cpath)

    b = np.ctypeslib.as_array(
    (ctypes.c_float * length).from_address(ctypes.addressof(res.contents)))

    return b


if __name__ == "__main__":

    length = 1000

    Aperture = np.ones(length) * 2
    Beta = np.random.normal(0.15, 0.12, size=length)
    HLR = np.ones(length) * 10**0.9
    n = np.ones(length) * 3
    SM = np.ones(length) * 11.3
    DMr = HLR * 20
    DMc = np.ones(length) * 1.3

    t1_start = time.time()

    FullVelocityDispersion(Aperture, Beta, HLR, n, SM, 0.0)

    t1_stop = time.time()

    print("Time:", t1_stop-t1_start)


