import numpy as np
from scipy import interpolate

def GetGridVelocityDispersion(beta, size, n, stellar_mass):
    """TODO Docstring"""
    points = (beta, size, n, stellar_mass)
    n = 20
    beta_range = np.linspace(-0.19, 0.49, n)
    sizes_range = np.linspace(2, 100, n)
    sersic_index_range = np.linspace(0.5, 10, n)
    stellar_mass_range = 10**np.linspace(5, 15, n)

    data = np.load("Grid.npy")

    return interpolate.interpn((beta_range, sizes_range, sersic_index_range, stellar_mass_range), data, points)


if __name__ == "__main__":
    res = GetGridVelocityDispersion(.1, 10, 4, 10**10)
    print(res)
