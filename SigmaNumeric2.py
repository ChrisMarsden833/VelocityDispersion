# Imports - Generic
import numpy as np
from tqdm import tqdm
import pickle

# Imports - More specialized
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.special import gamma, gammainc, gammaincc

# Imports - Astronomy
from colossus.cosmology import cosmology
from colossus.halo import profile_nfw
from colossus.halo import concentration

# Imports - Local
from SigmaLib import *

# Preprep
#  Set Cosmology
cosmo = cosmology.setCosmology("planck18")
#  Gravitational Constant
GR = 4.3009125e-6 # In units of kpc M_sun^-1 (km/s)^2

###################################################################################
# Generic Functions
def gamma_lower(a, z):
    """Lower gamma function, used in sersic profiles"""
    top = gammainc(a, z)
    cancel = gamma(a)
    return top * cancel
    
def tgamma(a, z):
    """gamma function used in sersic profiles"""
    top = gammaincc(a, z)
    cancel = gamma(a)
    return top * cancel

def sersicMassWithin(r, bulge_mass, re, n):
    """Get stellar mass within r based on a sersic profile"""
    b_n = 2.*n - (1./3.) + (.009876/n)
    p_n = 1. - .6097/n  + .00563/(n*n);
    a_s = re/(b_n**n)
    x = r/a_s
    threemp_term = (3.-p_n) * n
    res = 10.**bulge_mass 
    gamma1 = gamma_lower(threemp_term, x**(1./n))
    gamma2 = tgamma(threemp_term, 0.0)
    res *= gamma1/gamma2
    return res

def NFW_massWithin(Rs, rho, Re):
    """Mass within Re of halo with NFW profile"""
    x = Re/Rs
    res = 4 * np.pi * Rs**3 * rho * (np.log(1 + x) - x/(1+x))
    return res

def SDSS_Sizes_Fit(sm, z=0):
    """ Function to assign galaxy size according to mosleh fit, but with a redshift evolutionary term
    as shown in the velocity dispersion paper.
    params:
        sm (numpy array) : stellar mass, in log10 [m_sun]
        z (numpy array) : redshift [dimensionless]
    returns:
        numpy array, galaxy size [kpc]
    """
    args = [3.04965733e-03, 9.67029240e-02, -3.91008421e+01, 2.04401878e-01, -4.29464785e+00]
    def gammaFunc(sm, a1, a2, a3, a4, a5):
        return a1 * sm * ((sm * a2) ** a3 + (sm * a4) ** a5) ** -1

    smp = 10**(sm)
    res = (10**-0.314) * (smp**0.042) * (1 + smp/(10**10.537))**0.76
    
    isarray_sm = True if hasattr(sm, "__len__") else  False
    isarray_z = True if hasattr(z, "__len__")  else  False

    if isarray_sm and isarray_z:
        assert len(sm) == len(z), "sm and z are unequal lengths: {} and {} respectively".format(len(sm), len(z))

    gamma = gammaFunc(sm, *args)
    res = res * (1.+z)**-gamma
    return res

def SDSS_Sersic_Fit(sm, z = 0., minsm = 8, natmin = 1.6):
    res = 10**(-0.01072867 * sm**3 +\
                0.32824192 * sm**2 +\
               -3.18145729 * sm +\
               10.14652827 )
    isarray_sm = True if hasattr(sm, "__len__") else  False

    if isarray_sm:
        res[sm < minsm] = natmin
    else:
        if sm < minsm:
            res = natmin
    res *= (1 + z)**-1
    return res

def GetParameterFromMatrix(Matrix, x, y, datax, datay):
    """Function to interpolate a parameter from a 2d matrix given x and y coordinates and data.
    This is really just a wrapper function around the insane formatting that scipy's griddata requires.
    params:
        Matrix (2D numpy array) : the matrix itself.
        x, y (numpy arrays) : the coordinate data. The lengths must match the x and y dimensions of Matrix.
        datax, datay (numpy arrays) : the data to lookup.
    returns:
        (numpy array) : the values.
    """
    assert Matrix.shape[0] == len(y), "Matrix.shape[0] = {} != len(y) = {}".format(Matrix.shape[0], len(y))
    assert Matrix.shape[1] == len(x), "Matrix.shape[1] = {} != len(x) = {}".format(Matrix.shape[1], len(x))
    
    domainX = np.tile(x, len(y)).ravel()
    domainY = np.tile(np.vstack(y), len(x)).ravel()
    domain = (domainX, domainY)
    
    assert np.amin(datax) >= np.amin(x), "x data ({}) exists below domain ({})".format(np.amin(datax), np.amin(x) )
    assert np.amax(datax) <= np.amax(x), "x data ({}) exists above domain ({})".format(np.amax(datax), np.amax(x) )
    assert np.amin(datay) >= np.amin(y), "y data ({}) exists below domain ({})".format(np.amin(datay), np.amin(y) )
    assert np.amax(datay) <= np.amax(y), "y data ({}) exists above domain ({})".format(np.amax(datay), np.amax(y) )
    
    if hasattr(datax, "__len__"):
        K = np.ones_like(datax)
        for i, xi in enumerate(datax):
            sample = (xi, datay[i])
            K[i] = griddata(domain, Matrix.ravel(), sample, method='cubic')    
    else:
        sample = (datax, datay)
        K = griddata(domain, Matrix.ravel(), sample, method='cubic')   
    
    return K

def repeat_to_length(string_to_expand, length):
    """Simple function to tile a string"""
    return (string_to_expand * (int(length/len(string_to_expand))+1))[:length]

def halo_mass_to_stellar_mass(halo_mass, z, formula="Moster", scatter=0.11):
    """Function to generate stellar masses from halo masses.
    This is based on Grylls 2019, but also has the option to use the
    parameters from Moster. This is a simplified version of Pip's
    DarkMatterToStellarMass() function.
    :param halo_mass: array, of halo masses (log10)
    :param z: float, the value of redshift
    :param formula: string, the method to use. Options currently include "Grylls19" and "Moster"
    :param scatter: bool, to scatter or not
    :return array, of stellar masses (log10).
    """

    # If conditions to set the correct parameters.
    if formula == "Grylls19":
        z_parameter = np.divide(z - 0.1, z + 1)
        m_10, shm_norm_10, beta10, gamma10 = 11.95, 0.032, 1.61, 0.54
        m_11, shm_norm_11, beta11, gamma11 = 0.4, -0.02, -0.6, -0.1
    elif formula == "Moster":
        z_parameter = np.divide(z, z + 1)
        m_10, shm_norm_10, beta10, gamma10 = 11.590, 0.0351, 1.376, 0.608
        m_11, shm_norm_11, beta11, gamma11 = 1.195, -0.0247, -0.826, 0.329
    elif formula == "NoKnee":
        z_parameter = np.divide(z, z + 1)
        m_10, shm_norm_10, beta10, gamma10 = 11.590, 0.0351, 0.0, 0.608
        m_11, shm_norm_11, beta11, gamma11 = 1.195, -0.0247, -0.826, 0.329
    else:
        assert False, "Unrecognised formula"

    # Create full parameters
    m = m_10 + m_11 * z_parameter
    n = shm_norm_10 + shm_norm_11 * z_parameter
    b = beta10 + beta11 * z_parameter
    g = gamma10 + gamma11 * z_parameter
    # Full formula
    internal_stellar_mass = np.log10(np.power(10, halo_mass) *\
                                     (2 * n * np.power((np.power(np.power(10, halo_mass - m), -b)\
                                                        + np.power(np.power(10, halo_mass - m), g)), -1)))

    if formula == "Grylls19":
        internal_stellar_mass -= 0.1

    # Add scatter, if requested.
    if not scatter == False:
        internal_stellar_mass += np.random.normal(scale=scatter, size=np.shape(internal_stellar_mass))

    return internal_stellar_mass

###################################################################################
# Class that stores/runs numerical sigma code.

class NumericalSigma:
    def __init__(self, terminal_output = True):
        self.K_matrix = None
        self.L_matrix = None
        self.terminal_output = terminal_output
        
    def GenerateKmatrix(self, re_domain, n_domain, fix_ap = 1/8.):
        """Member function to generate the K matrix. Requires no prerequistes.
        Params:
            ap_domain, n_domain (numpy arrays) : the aperture/re and sersic index domain to generate the matrix values at.
            The size of these arrays will dictate the size of the generated matrix.
        returns:
            none - values are stored in the class.
        """
        z = 0
        
        if self.terminal_output:
            print("##############################################################################")
            print("Generating matrix for K of shape (x={}, y={}).".format(len(re_domain), len(n_domain)))
        
        self.fix_ap = fix_ap
        
        # Example stellar mass, example halo mass
        mstar = 10.
    
        self.re_domain_Kx = re_domain
        self.n_domain_Ky = n_domain
    
        first = True
        for n in tqdm(self.n_domain_Ky, disable = not self.terminal_output): # Iterate row by row.
                        
            sigma = Sigma(fix_ap * re_domain,
                        Bulge_mass = mstar,
                        Bulge_Re = re_domain,
                        Bulge_n = n,
                        Bulge_Beta = 0.0,
                        DarkMatter_type = 'None',
                        debug = False,
                        threads = 8,
                        library_path = "./lib/libsigma.so")
            
            LHS = GR * sersicMassWithin(re_domain, mstar, re_domain, n)/(re_domain)
            
            K = LHS/sigma**2
            
            if first:
                Kmatrix = K
                first = False
            else:
                Kmatrix = np.vstack((Kmatrix, K))
                
        self.K_matrix = Kmatrix
        
        if self.terminal_output:
            print("K matrix generated.")
        
        
    def K_LaTeX_markup(self):
        """Member function to print the latex markdown for displaying the K matrix as a table."""
        
        
        alignment = "c|c" # label, n, value
        
        midpointy = int(len(self.n_domain_Ky)/2)
        ylabel = "n"
        
        print("\\begin{array}{" + alignment + "}")
        
        first_row = "n & \mathcal{K} \\\ "
        print(first_row)
        
        for i in range(len(self.n_domain_Ky)):
            string = str(self.n_domain_Ky[i]) + " & " + str(np.round(self.K_matrix[i, 0], 2)) + " \\\ "
            print(string)
            
        print("\\end{array}")
        
    def M_LaTeX_markup(self):
        """Member function to print the latex markdown for displaying the M matrix as a table."""
        
        
        alignment = "c|c" # label, n, value
        
        midpointy = int(len(self.Beta_range)/2)
        
        print("\\begin{array}{" + alignment + "}")
        
        first_row = "\\beta & \mathcal{N} \\\ "
        print(first_row)
        
        for i in range(len(self.Beta_range)):
            string = str(self.Beta_range[i]) + " & " + str(np.round(self.M[i], 2)) + " \\\ "
            print(string)
            
        print("\\end{array}")
  

    def L_LaTeX_markup(self):
        """Member function to print the latex markdown for displaying the M matrix as a table."""
        
        print_cols = len(self.c_domain_Lx) + 2 # 1 for the column for the re labels, 1 for the label
        alignment = repeat_to_length("c", print_cols) # Centre alignment
        
        midpointx = int(len(self.c_domain_Lx)/2 + 2)
        midpointy = int(len(self.n_domain_Ly)/2)
        xlabel = "c"
        ylabel = "n"
        
        
        print("\\begin{array}{" + alignment + "}")
        
        first_row = ""
        for i in range(print_cols):
            if i+1 == midpointx:
                first_row += xlabel + " & "
            elif i == print_cols-1: # Last
                first_row += "\\\ "
                
            else:
                first_row += " & "
        
        print(first_row)
        
        second_row = ""
        for i in range(print_cols):
            if i+1 == print_cols:
                terminator = " \\\ "
            else:
                terminator = " & "
            
            if i == 0 or i == 1:
                second_row += terminator
            else:
                second_row += str(np.round(self.c_domain_Lx[i-2], 2)) + terminator
                
        print(second_row)
        
        for i in range(len(self.n_domain_Ly)):
            
            string = ""
            
            if i == midpointy:
                string += ylabel + " & " + str(np.round(self.n_domain_Ly[i], 2))  +  " & "
            else:
                string += " & " + str(np.round(self.n_domain_Ly[i], 2))  +  " & "
            
            for j in range(len(self.c_domain_Lx)):
                terminator = " & "
                
                if j == len(self.c_domain_Lx)-1:
                    terminator = "\\\ "
                
                string += str(np.round(self.L_matrix[i, j], 2)) + terminator
                
            print(string)
        
        print("\\end{array}")
  
         

        
    def GenerateLmatrix(self, c_domain, n_domain):
        
        assert self.K_matrix is not None, "K matrix must be calculated first."
        
        if self.terminal_output:
            print("##############################################################################")
            print("Generating matrix for L of shape (x={}, y={}).".format(len(c_domain), len(n_domain)))
            
        self.c_domain_Lx = c_domain
        self.n_domain_Ly = n_domain
        
        mstar = 10.0
        mhalo = 11.5 
    
        z=0

        re = SDSS_Sizes_Fit(mstar) #*(1+z)**-0.5
        
        withstars = True
        first = True
        for n in tqdm(n_domain):
                           
            nfw  = profile_nfw.NFWProfile(M = 1E12, c = 10.0, z = z, mdef = '200c') # Just for obj
            rho, rs = nfw.fundamentalParameters( (10**mhalo)*cosmo.h, c_domain, z, '200c')
            rs /= cosmo.h
            rho *= cosmo.h**2
            

            
            sigma = Sigma(self.fix_ap * re * np.ones_like(c_domain),
                            Bulge_mass = mstar,
                            Bulge_Re = re,
                            Bulge_n = n,
                            Bulge_Beta = 0.0,
                            DarkMatter_type = 'NFW',
                            HaloRs = rs,
                            HaloRhos = rho,
                            debug = False,
                            StarsOn=withstars,
                            threads = 8,
                            library_path = "./lib/libsigma.so")
            
            K = 0
            sm = 0
            if withstars:
                K = GetParameterFromMatrix(self.K_matrix, self.re_domain_Kx, self.n_domain_Ky, re, n)
                sm = sersicMassWithin(re, mstar, re, n)
                        
            dm = NFW_massWithin(rs, rho, re)
            
            tot = sm + dm
                
            LHS = GR * tot/re
            
            L = (LHS/(sigma**2) - K*sm/tot) * tot/dm
            
            if first:
                Lmatrix = L
                first = False
            else:
                Lmatrix = np.vstack((Lmatrix, L))
                    
        self.L_matrix = Lmatrix
        
        if self.terminal_output:
            print("L matrix generated.")
            
    def GenerateMVector(self, Beta_range):
        
        self.Beta_range = Beta_range
        
        mstar = 10.0
        mhalo = 11.5 
        z=0
        re = SDSS_Sizes_Fit(mstar) #*(1+z)**-0.5
        n = SDSS_Sersic_Fit(mstar)
        c = concentration.concentration((10**mhalo)*cosmo.h, '200c', z=z, model = 'ishiyama20')
        nfw  = profile_nfw.NFWProfile(M = 1E12, c = 10.0, z = z, mdef = '200c') # Just for obj
        rho, rs = nfw.fundamentalParameters( (10**mhalo)*cosmo.h, c, z, '200c')
        rs /= cosmo.h
        rho *= cosmo.h**2
                        
        sigma = Sigma(self.fix_ap * re * np.ones_like(Beta_range),
                        Bulge_mass = mstar,
                        Bulge_Re = re,
                        Bulge_n = n,
                        Bulge_Beta = Beta_range,
                        DarkMatter_type = 'NFW',
                        HaloRs = rs,
                        HaloRhos = rho,
                        debug = False,
                        threads = 8,
                        library_path = "./lib/libsigma.so")
        
        K = GetParameterFromMatrix(self.K_matrix, self.re_domain_Kx, self.n_domain_Ky, re, n)
        sm = sersicMassWithin(re, mstar, re, n)
                        
        dm = NFW_massWithin(rs, rho, re)
        L = GetParameterFromMatrix(self.L_matrix, self.c_domain_Lx, self.n_domain_Ly, c, n)
            
        tot = sm + dm
                
        LHS = GR * tot/re
        
        M = LHS/(sigma**2) - (L*dm/tot + K *sm/tot)
        
        self.M = M
        
        
        
        
        
            
    def SigmaNumeric(self, ApertureSize, Bulge_mass, Bulge_Re, Bulge_n, Bulge_Beta, HaloMass, z, dm='NFW'):
        """Function to return sigma based on the numerical matricies"""
        
        
        
            
            
        K = GetParameterFromMatrix(self.K_matrix, self.re_domain_Kx, self.n_domain_Ky, Bulge_Re, Bulge_n)
        
        if dm == 'NFW':
            c = concentration.concentration((10**HaloMass)*cosmo.h, '200c', z=z, model = 'ishiyama20')
            
            print("Haloc", c)
            
            nfw = profile_nfw.NFWProfile(M = 1E12, c = 10, z = z, mdef='200c') # Just need the obj
            rho, rs = nfw.fundamentalParameters( (10**HaloMass)*cosmo.h, c, z, '200c')
            rs /= cosmo.h
            rho *= cosmo.h**2
            
            #c2L = interp1d(self.c_rangeL, self.L)
            #L = c2L(c)
        
            L = GetParameterFromMatrix(self.L_matrix, self.c_domain_Lx, self.n_domain_Ly, c, Bulge_n)
            dm = NFW_massWithin(rs, rho, Bulge_Re)
        else:
            L = 0
            dm = 0
            
            
        beta2m = interp1d(self.Beta_range, self.M)
        M = beta2m(Bulge_Beta)
            
            
        #z2m = interp1d(self.z_range, self.M)
        #M = z2m(z)
        
        #print(self.M)
        
        #M = np.zeros_like(K)
       # M = GetParameterFromMatrix(self.M_Matrix, self.z_range_Mx, self.mhalo_range_My, z, HaloMass)
        #M[HaloMass > 0] = Mall[HaloMass > 0]
        
    
        sm = sersicMassWithin(Bulge_Re, Bulge_mass, Bulge_Re, Bulge_n)
        tot = dm + sm
        
        LHS = GR * tot/Bulge_Re
        
        RHS = (L * dm/tot + K * sm/tot + M)
        
        
        
        #RHS = (L + K)
        
        sigma = LHS/RHS
        
        return np.sqrt(sigma)
                    
        
if __name__ == "__main__":
    obj = NumericalSigma()
    
    re_domain = np.arange(1, 10, 1.0)
    n_domain = np.arange(1, 9.0, 1)
    
    obj.GenerateKmatrix(re_domain, n_domain, 0.25)
    
   
    
    c_domain = np.arange(5, 15, 1.0)
    
    obj.GenerateLmatrix(c_domain, n_domain)

    
    beta_range = np.array([-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4]) # np.linspace(-0.2, 0.49, 5)
    obj.GenerateMVector(beta_range)
    
    #mhalo_range = np.linspace(8, 15.5, 20)
    #obj.GenerateLmatrix(c_domain, mhalo_range)
    
    
    
    #mstar_range = np.linspace(7, 13, 10)
    
    #fixed_aperture = 1/8.
    #fixed_beta = 0.0
    
    #mhalo_range = np.linspace(0, 16, 20)
    #mstar_range = np.linspace(5, 13, 20)
    
    #mh_range = np.linspace(10, 15.5, 10)
    #z_range = np.linspace(0, 4, 10)
    #obj.GenerateLmatrix(mh_range, z_range)
   
    #z_domain = np.linspace(0, 5, 10)
   
    #obj.GenerateMMatrix(z_domain, mhalo_range)
    
    filehandler = open("SigmaNumeric.pkl", "wb")
    pickle.dump(obj, filehandler)
    filehandler.close()
    

    #obj.L_LaTeX_markup()
    
    #obj.K_LaTeX_markup()
    
    
    
        
        
     