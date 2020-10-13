import numpy as np
from scipy import interpolate

def GetGridVelocityDispersion(beta, size, n, stellar_mass):
    """TODO Docstring"""
    points = (beta, size, n, stellar_mass)
    n = 20
    beta_range = np.linspace(-0.19, 0.49, 10)
    sizes_range = np.linspace(0.1, 100, 20)
    sersic_index_range = np.linspace(0.5, 10, 10)
    stellar_mass_range = np.linspace(10**5, 10**15, 20)


    data = np.load("Grid.npy")

    return interpolate.interpn((beta_range, sizes_range, sersic_index_range, stellar_mass_range), data, points)


if __name__ == "__main__":
    res = GetGridVelocityDispersion(.1, 10, 4, 10**10)
    print(res)
