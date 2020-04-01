import numpy as np
import itertools
from VelocityDispersionLibrary import *

n = 2

aperture = 0.5

beta_range = np.linspace(-0.19, 0.49, n)
sizes_range = np.linspace(2, 100, n)
sersic_index_range = np.linspace(0.5, 10, n)
stellar_mass_range = 10**np.linspace(5, 15, n)


print(beta_range)

res = list(itertools.product(    beta_range,
                                sizes_range,
                         sersic_index_range,
                         stellar_mass_range))

npres = np.array(res)


Beta_grid = npres[:, 0]
Sizes_grid = npres[:, 1]
SersicIndex_grid = npres[:, 2]
StellarMass_grid = npres[:, 3]

ApertureSize = np.ones_like(Beta_grid) * aperture
HaloRs = 20* Sizes_grid
HaloC = np.ones_like(Beta_grid)
OmegaM = 0.29
H = 299.

VD = FullVelocityDispersion(ApertureSize, Beta_grid, Sizes_grid, SersicIndex_grid, StellarMass_grid, HaloRs, HaloC, OmegaM, H)

theGrid = np.zeros((n, n, n, n))

print(" Sorting into Array")

for i, VD_element in enumerate(VD):
    BetaCoord = list(beta_range).index(Beta_grid[i])
    HLRCoord = list(sizes_range).index(Sizes_grid[i])
    SersicCoord = list(sersic_index_range).index(SersicIndex_grid[i])
    SMCoord = list(stellar_mass_range).index(StellarMass_grid[i])
    theGrid[BetaCoord, HLRCoord, SersicCoord, SMCoord] = VD_element

np.save("Grid.npy", theGrid)

print("Completed!")