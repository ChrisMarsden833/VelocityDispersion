# Christopher Marsden, The University of Southampton

# imports
import numpy as np
import ctypes
import time
#from datetime import timedata
import os
# openmp can sometimes be already using by a jupyter notebook. This command prevents crashing.
os.environ['KMP_DUPLICATE_LIB_OK']='True' # Potentially dangerous.

def Sigma(ApertureSize,
            Bulge_mass,
              Bulge_Re,
               Bulge_n,
            Bulge_Beta=0.0,
    StellarReminantPref=1.0,
             Disk_mass=0.0,
     Disk_scale_length=0.0,
      Disk_inclination=0.0,
                     z=0,
       DarkMatter_type=None,
                HaloRs=0.0,
              HaloRhos=0.0,
           BlackHoleOn=False,
         BlackHoleMass=0.0,
               StarsOn=True,
                  mode=1,
                 debug=True,
               threads=-1,
          library_path="/home/chris/Files/ProjectSigma/VelocityDispersion/lib/libsigma.so"):

    """The Velocity Dispersion of a Galaxy (sigma).

    This is the 'master' wrapper function that allows the average python user to intuitively run Chris Marsden's
    velocity dispersion code without necessarily understanding C++. This function basically takes the desired inputs
    as numpy arrays or floats, transforms them into ctypes arrays and calls the sigma code.

    Please note that if using a jupyter notebook, the printf's from the C++ (basically just the progress of the code)
    will be piped to the notebook server terminal. A bit of work could probably fix this.

    It's okay to use single values for many arguments if they are going to be the same - this function will automatically
    convert them into tiled arrays. The first argument, Aperture, must always be the maximum length. A single value is okay.

    A galaxy is split into a bulge and disk component for the purposes of this code. For full Ellipticals or full disks,
    simply set the non-existant parameters to zero.

    Parameters:
    -----------
        ApertureSize (np array or float) : This is the aperture within which sigma will be measured, measured in kpc (not 
            arcseconds!). When the function mode is changed to mode = 2, this value is actually the radius at which sigma
            should be measured. [kpc]

        Bulge_mass (np array or float) : The Stellar Mass of the bulge [log10 Msun]
        Bulge_Re (np array or float) : The value of Re for the bulge component [kpc]
        Bulge_n (np array or float) : The Sersic Index of the bulge component [dimensionless] 
        Bulge_Beta (np array or float) : The velocity Anisotropy [dimensionless]

        Disk_mass (np array or float) : The Stellar Mass of the disk, default to 0 [log10 Msun]
        Disk_scale_length (np array or float) : The scale length R_D of the disk, default is 0 [kpc]
        Disk_inclination (np array or float) : The angle at which the disk is observed relative to aperture, 
            default is 0 [degrees]

        z (np array or float) : The redshift of the galaxy [dimensionless]
        
        DarkMatter_type (String) : String description of the dark matter structure. Default is 'None', which will turn off
            Dark Matter altogether (rendering the HaloMass and HaloC arguments redundant if assigned). The only currently 
            implemented profile is 'NFW'. This string is case-sensitive.
        HaloMass (np array or float) : The mass of the dark matter halo [log10 Msun]
        HaloC (np array or float) : The halo concentration parameter [dimensionless]
        HaloRs (np array or float): The halo scale radius [kpc]
        HaloRhos (np array or float): The halo density factor [M_sun]/kpc^3

        BlackHoleOn (bool) : If the gravitational effects of a black hole should be switched on. Defaults to False. If 
            False, the variable BlackHoleMass is redundant if assigned.
        BlackHoleMass (np array or float) : The black hole mass [log10 Msun]

        StarsOn (bool) : If the gravitational effects of the stars themselves should be switched on (used for testing). 
            defaults to True.
        mode (integer) : A backend parameter that controls how the code works. Normally, users should ignore this.
            mode = 1 (default) : Normal Behavior. Velocity dispersion is computed within an aperture.
            mode = 2 : Return LOS velocity dispersion as a function of r. In this case, ApertureSize becomes r.
        threads (integer) : The number of threads to use for multiprocessing. The default is -1, which will use the
            maximum number of threads avaiilable. Multiprocessing is performed by assigning galaxies to cores using opm.
        debug (bool) : Backend variable to activate debugging print statements. Default is False.
        library_path (string) : The (full) path to the DLL file. This defaults to:
            "/Users/chris/Documents/PhD/ProjectSigma/VelocityDispersion/lib/libsigma.so"
        
    Returns:
    --------
    (np array) : The velocity dispersions of the specified galaxies [kms^-1]

    """

    start_time = time.time() # Time the code

    ibc = ctypes.CDLL(library_path) # link to the DLL.

    # Aperture is the first variable, so we assume that it has the highest length. If it's a float, that's okay.
    if hasattr(ApertureSize, "__len__"):
        length = int(len(np.array(ApertureSize)))
    else:
        length = 1

    if debug:
        print("Length of arrays: ", length)

    # define some types for later use
    c_float_p = ctypes.POINTER(ctypes.c_float)
    c_int_p = ctypes.POINTER(ctypes.c_int32)

    # Process all the variables. This combinations of functions makes them an array if they are not,
    # and makes them valid for ctypes.
    ApertureSize = check_make_array(ApertureSize, length).astype(np.float32).ctypes.data_as(c_float_p)
    Bulge_Beta = check_make_array(Bulge_Beta, length).astype(np.float32).ctypes.data_as(c_float_p)
    Bulge_Re = check_make_array(Bulge_Re, length).astype(np.float32).ctypes.data_as(c_float_p)
    Bulge_n = check_make_array(Bulge_n, length).astype(np.float32).ctypes.data_as(c_float_p)
    Bulge_mass = check_make_array(Bulge_mass, length).astype(np.float32).ctypes.data_as(c_float_p)
    StellarReminantPref = check_make_array(StellarReminantPref, length).astype(np.float32).ctypes.data_as(c_float_p)
    Disk_mass = check_make_array(Disk_mass, length).astype(np.float32).ctypes.data_as(c_float_p)
    Disk_inclination = check_make_array(Disk_inclination, length).astype(np.float32).ctypes.data_as(c_float_p)
    Disk_scale_length = check_make_array(Disk_scale_length, length).astype(np.float32).ctypes.data_as(c_float_p)
    z = check_make_array(z, length).astype(np.float32).ctypes.data_as(c_float_p)
    HaloRs = check_make_array(HaloRs, length).astype(np.float32).ctypes.data_as(c_float_p)
    HaloRhos = check_make_array(HaloRhos, length).astype(np.float32).ctypes.data_as(c_float_p)
    BlackHoleMass = check_make_array(BlackHoleMass, length).astype(np.float32).ctypes.data_as(c_float_p)

    # Manage which components are gravitationally active, based on the inputs provided.
    if StarsOn:
        stars_component = 1
    else:
        stars_component = 0
    if (DarkMatter_type is None) or (DarkMatter_type is "None"):
        dm_component = 0
    else:
        dm_component = 1
    if BlackHoleOn:
        bh_component = 1
    else:
        bh_component = 0
    # We assemble a 'component array', which will be passed to the library to manage this succinctly.
    component_array = [stars_component, dm_component, bh_component]
    if debug:
        print("Component array [stars, dark_matter, black_hole] ", component_array)

    # Force mode and threads into integers, just in case they are not 
    # (prevents pointless errors if these are interpreted as floats)
    mode = int(mode)
    threads = int(threads)

    # Sombody might put DarkMatter_type as python None rather than the string, so we deal with that here.
    if DarkMatter_type is None:
        DarkMatter_type = "None"

    DarkMatter_type = DarkMatter_type.encode('utf-8') # Encode the dm string as a c-style string.
    component_array = np.array(component_array).astype(np.int32).ctypes.data_as(c_int_p) # Encode component array

    if debug:
        debug_var = 1
    else:
        debug_var = 0

    # Set the argument types for the function. Note the order is different here.
    ibc.ParallelSigma.argtypes = [ctypes.POINTER(ctypes.c_float),  # Aperture
                                  ctypes.POINTER(ctypes.c_float),  # Redshift
                                  ctypes.POINTER(ctypes.c_float),  # Bulge mass
                                  ctypes.POINTER(ctypes.c_float),  # Bulge radius
                                  ctypes.POINTER(ctypes.c_float),  # Bulge beta
                                  ctypes.POINTER(ctypes.c_float),  # Bulge sersic index
                                  ctypes.POINTER(ctypes.c_float),  # Stellar Remenant Prefactor
                                  ctypes.POINTER(ctypes.c_int),  # Bulge component flag
                                  ctypes.POINTER(ctypes.c_float),  # Disk Mass
                                  ctypes.POINTER(ctypes.c_float),  # Disk Inclination
                                  ctypes.POINTER(ctypes.c_float),  # Disk ScaleLength
                                  ctypes.c_char_p, # Profile Name
                                  ctypes.POINTER(ctypes.c_float), # Halo Rs
                                  ctypes.POINTER(ctypes.c_float), # Halo Rhos
                                  ctypes.POINTER(ctypes.c_float),  # Black Hole Mass
                                  ctypes.c_int32, # Threads
                                  ctypes.c_int32, # Mode
                                  ctypes.c_int32, # debug
                                  ctypes.c_int32]  # Length

    # Set the argument return type. 
    ibc.ParallelSigma.restype = ctypes.POINTER(ctypes.c_float)

    # Call the function! Note argument order is different
    res = ibc.ParallelSigma(ApertureSize, z, Bulge_mass, Bulge_Re, Bulge_Beta,
            Bulge_n, StellarReminantPref, component_array, Disk_mass, Disk_inclination, Disk_scale_length,
            DarkMatter_type, HaloRs, HaloRhos, BlackHoleMass, mode, threads, debug_var, length)

    # Convert the returned variable into readable format.
    res = np.ctypeslib.as_array((ctypes.c_float * length).from_address(ctypes.addressof(res.contents)))

    if debug:
        pass
        #print("elapsed time : ",  str(timedelta(seconds=time.time() - start_time)))

    return res

def check_list(listoflists, names):
    """ Simply function to check that a list of lists are all consistent lengths.
    Arguments
    ---------
    Returns
    -------
        length (float) : the length of (all) the lists, else assertion error.
    """
    for i, element in enumerate(listoflists):
        if i != 0:
            assert len(element) == previous_length, """Supplied lengths inconsistent - {} has {} elements, wheras
                 {} has {} elements""".format(names[i - 1], previous_length, names[i], len(element))
        previous_length = len(element)
        length = len(element)
    return length

def check_make_array(subject, length):
    """ Function to check if a variable is an array or not - if it's not, we make it a tiled (np) array of correct length.
    Arguments
    ---------
        subject (listlike, or single value) : the variable to check
        length (int) : the length of the (potentially desired) array
    Returns
    -------
        subject (numpy array) : the variable, now in array format.

    """
    if not hasattr(subject, "__len__"):
        subject = np.ones(length) * float(subject)
    else:
        subject = np.array(subject)
        assert len(subject) == length, "Length not consistent."
    return subject


