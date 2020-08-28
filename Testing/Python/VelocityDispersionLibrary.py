import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import ctypes
from colossus.cosmology import cosmology
from scipy import interpolate
from scipy.interpolate import griddata
import time
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True' # Potentially dangerous

cosmo = cosmology.setCosmology('planck18')
omega_m = cosmo.Om0
H = cosmo.H0
ibc = ctypes.CDLL("/Users/chris/Documents/PhD/ProjectSigma/VelocityDispersion/lib/libsigma.so")

def Sigma(Aperture, Beta, HalfLightRadius, SersicIndex, StellarMass, z,
             DM=None, HaloMass=None,
             BlackHole=False, BHMass=None,
             stars=True,
             disk_mass = 0.0,
             cpath="../data/cM_planck18.txt"):

    length = int(len(Aperture))
    c_float_p = ctypes.POINTER(ctypes.c_float)
    c_int_p = ctypes.POINTER(ctypes.c_int32)

    Aperture = Aperture.astype(np.float32).ctypes.data_as(c_float_p)
    Beta = check_make_array(Beta, length).astype(np.float32).ctypes.data_as(c_float_p)
    HalfLightRadius = check_make_array(HalfLightRadius, length).astype(np.float32).ctypes.data_as(c_float_p)
    SersicIndex = check_make_array(SersicIndex, length).astype(np.float32).ctypes.data_as(c_float_p)
    StellarMass = check_make_array(StellarMass, length).astype(np.float32).ctypes.data_as(c_float_p)
    disk_mass = check_make_array(disk_mass, length).astype(np.float32).ctypes.data_as(c_float_p)
    z = check_make_array(z, length).astype(np.float32).ctypes.data_as(c_float_p)
    if HaloMass is None:
        HaloMass = 0.0

    HaloMass = check_make_array(HaloMass, length).astype(np.float32).ctypes.data_as(c_float_p)
    if BHMass is None:
        BHMass = 0.0
    BHMass = check_make_array(BHMass, length).astype(np.float32).ctypes.data_as(c_float_p)
    cpath = cpath.encode('utf-8')

    if stars:
        stars_component = 1
    else:
        stars_component = 0
    if (DM is None) or (DM is "None"):
        dm_component = 0
    else:
        dm_component = 1
    if BlackHole:
        bh_component = 1
    else:
        bh_component = 0

    component_array = [stars_component, dm_component, bh_component]

    print(component_array)

    if DM is None:
        DM = "None"
    DM = DM.encode('utf-8')

    component_array = np.array(component_array).astype(np.int32).ctypes.data_as(c_int_p)

    ibc.ParallelSigma.argtypes = [ctypes.POINTER(ctypes.c_float), # Aperture
                                  ctypes.POINTER(ctypes.c_float), # Redshift
                                  ctypes.POINTER(ctypes.c_float), # Bulge mass
                                  ctypes.POINTER(ctypes.c_float), # Bulge radius
                                  ctypes.POINTER(ctypes.c_float), # Bulge beta
                                  ctypes.POINTER(ctypes.c_float), # Bulge sersic index
                                  ctypes.POINTER(ctypes.c_int), # Bulge component flag
                                  ctypes.POINTER(ctypes.c_float), # Disk Mass
                                  ctypes.POINTER(ctypes.c_float), # Halo Mass
                                  ctypes.c_char_p, # Profile Name
                                  ctypes.c_char_p, # c_path
                                  ctypes.POINTER(ctypes.c_float),  # Black Hole Mass
                                  ctypes.c_int32]  # Length

    ibc.ParallelSigma.restype = ctypes.POINTER(ctypes.c_float)
    res = ibc.ParallelSigma(Aperture, z, StellarMass, HalfLightRadius, Beta, SersicIndex, component_array, disk_mass, HaloMass, DM, cpath, BlackHole, length)

    b = np.ctypeslib.as_array(
    (ctypes.c_float * length).from_address(ctypes.addressof(res.contents)))

    return b

'''
def SigmaLOS(R, Beta, HalfLightRadius, SersicIndex, StellarMass, z,
             DM=None, HaloMass=None,
             BlackHole=False, BHMass=None,
             stars=True,
             cpath="../data/cM_planck18.txt"):

    length = int(len(R))
    c_float_p = ctypes.POINTER(ctypes.c_float)
    c_int_p = ctypes.POINTER(ctypes.c_int32)

    R = R.astype(np.float32).ctypes.data_as(c_float_p)
    Beta = check_make_array(Beta, length).astype(np.float32).ctypes.data_as(c_float_p)
    HalfLightRadius = check_make_array(HalfLightRadius, length).astype(np.float32).ctypes.data_as(c_float_p)
    SersicIndex = check_make_array(SersicIndex, length).astype(np.float32).ctypes.data_as(c_float_p)
    StellarMass = check_make_array(StellarMass, length).astype(np.float32).ctypes.data_as(c_float_p)
    z = check_make_array(z, length).astype(np.float32).ctypes.data_as(c_float_p)
    if HaloMass is None:
        HaloMass = 0.0

    HaloMass = check_make_array(HaloMass, length).astype(np.float32).ctypes.data_as(c_float_p)
    if BHMass is None:
        BHMass = 0.0
    BHMass = check_make_array(BHMass, length).astype(np.float32).ctypes.data_as(c_float_p)
    cpath = cpath.encode('utf-8')


    if stars:
        stars_component = 1
    else:
        stars_component = 0
    if (DM is None) or (DM is "None"):
        dm_component = 0
    else:
        dm_component = 1
    if BlackHole:
        bh_component = 1
    else:
        bh_component = 0

    component_array = [stars_component, dm_component, bh_component]

    print(component_array)

    if DM is None:
        DM = "None"
    DM = DM.encode('utf-8')

    component_array = np.array(component_array).astype(np.int32).ctypes.data_as(c_int_p)

    ibc.ParallelSigmaLos.argtypes = [ctypes.POINTER(ctypes.c_float),
                                  ctypes.POINTER(ctypes.c_float),
                                  ctypes.POINTER(ctypes.c_float),
                                  ctypes.POINTER(ctypes.c_float),
                                  ctypes.POINTER(ctypes.c_float),
                                  ctypes.POINTER(ctypes.c_float),
                                  ctypes.POINTER(ctypes.c_float),
                                  ctypes.POINTER(ctypes.c_float),
                                  ctypes.c_int32,
                                  ctypes.c_char_p,
                                  ctypes.c_char_p,
                                  ctypes.POINTER(ctypes.c_int32)]

    ibc.ParallelSigmaLos.restype = ctypes.POINTER(ctypes.c_float)

    res = ibc.ParallelSigmaLos(R, Beta, HalfLightRadius, SersicIndex, StellarMass, HaloMass, BHMass, z, length,
                               DM, cpath, component_array)

    b = np.ctypeslib.as_array((ctypes.c_float * length).from_address(ctypes.addressof(res.contents)))

    return b


def CumulativeMass(R, Beta, HalfLightRadius, SersicIndex, StellarMass, z, DM = None, HaloMass = None,
             cpath = "../data/cM_planck18.txt" , flag = 0):
    # Test that the supplied lengths are consistent
    variables = [R, Beta, HalfLightRadius, SersicIndex, StellarMass]
    variable_names = ["R", "Beta", "Half Light Radius", "Sersic Index", "Stellar Mass"]
    length = check_list(variables, variable_names)

    assert np.sum(HalfLightRadius <= 0) == 0, "HalfLightRadius has elements < 0, {}".format(
        HalfLightRadius[HalfLightRadius <= 0])
    assert np.sum(SersicIndex <= 0) == 0, "Sersic Index has elements < 0 {}".format(SersicIndex[SersicIndex <= 0])

    c_float_p = ctypes.POINTER(ctypes.c_float)

    R = R.astype(np.float32)
    c_R = R.ctypes.data_as(c_float_p)
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
    c_size = int(length)
    c_flag = int(flag)

    ibc.GetCumMass.argtypes = [ctypes.POINTER(ctypes.c_float),
                                     ctypes.POINTER(ctypes.c_float),
                                     ctypes.POINTER(ctypes.c_float),
                                     ctypes.POINTER(ctypes.c_float),
                                     ctypes.POINTER(ctypes.c_float),
                                     ctypes.POINTER(ctypes.c_float),
                                     ctypes.POINTER(ctypes.c_float),
                                     ctypes.c_int32,
                                     ctypes.c_char_p,
                                     ctypes.c_char_p,
                                     ctypes.c_int32]

    ibc.GetCumMass.restype = ctypes.POINTER(ctypes.c_float)

    res = ibc.GetCumMass(c_R, c_beta, c_hlr, c_n, c_sm, c_hm, c_z, c_size, c_DM, c_cpath, c_flag)

    b = np.ctypeslib.as_array((ctypes.c_float * length).from_address(ctypes.addressof(res.contents)))

    return b




def DM_profile(r, HaloMass, z, DM, cpath = "../data/cM_planck18.txt"):

    length = len(r)

    c_float_p = ctypes.POINTER(ctypes.c_float)

    c_DM = DM.encode('utf-8')
    c_cpath = cpath.encode('utf-8')

    r = r.astype(np.float32)
    c_r = r.ctypes.data_as(c_float_p)

    HaloMass = float(HaloMass)
    c_hm = HaloMass
    c_z = float(z)

    c_size = int(length)


    ibc.DM_Profile.argtypes = [ctypes.POINTER(ctypes.c_float),
                               ctypes.c_float,
                               ctypes.c_float,
                               ctypes.c_int32,
                               ctypes.c_char_p,
                               ctypes.c_char_p]

    ibc.DM_Profile.restype = ctypes.POINTER(ctypes.c_float)

    res = ibc.DM_Profile(c_r, c_hm, c_z, c_size, c_DM, c_cpath)

    b = np.ctypeslib.as_array( (ctypes.c_float * length).from_address(ctypes.addressof(res.contents)))

    return b


'''


def check_list(list, names):
    for i, element in enumerate(list):
        if i != 0:
            assert len(element) == previous_length, """Supplied lengths inconsistent - {} has {} elements, wheras 
                {} has {} elements""".format(names[i - 1], previous_length, names[i], len(element))
        previous_length = len(element)
        length = len(element)
    return length

def check_make_array(subject, length):
    if not hasattr(subject, "__len__"):
        subject = np.ones(length) * float(subject)
    else:
        assert len(subject) == length, "Length not consistent."
    return subject


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


