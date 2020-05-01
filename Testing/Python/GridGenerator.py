import numpy as np
import itertools
from VelocityDispersionLibrary import *

n = 1

aperture = 0.5 # Aperture size is FIXED - but it could be included in the grid quite easily.

# First set up the ranges we want to explore.
beta_range = np.linspace(-0.19, 0.49, n)
sizes_range = np.linspace(0.1, 100, 2*n)
sersic_index_range = np.linspace(0.5, 10, n)
stellar_mass_range = np.linspace(10**5, 10**15, 2*n)

# Next we combine these into a grid of all possible combinations.
res = list(itertools.product(    beta_range,
                                sizes_range,
                         sersic_index_range,
                         stellar_mass_range))
npres = np.array(res)

# Extract these into their respective arrays
Beta_grid = npres[:, 0]
Sizes_grid = npres[:, 1]
SersicIndex_grid = npres[:, 2]
StellarMass_grid = npres[:, 3]

# Booring constants
ApertureSize = np.ones_like(Beta_grid) * aperture
HaloRs = 20* Sizes_grid
HaloC = np.ones_like(Beta_grid)
OmegaM = 0.29
H = 299.

# Call the function that will actually create the array
VD = FullVelocityDispersion(ApertureSize, Beta_grid, Sizes_grid, SersicIndex_grid, StellarMass_grid, HaloRs, HaloC, OmegaM, H)

theGrid = np.zeros((len(beta_range), len(sizes_range), len(sersic_index_range), len(stellar_mass_range)))

print(" Sorting into Array")

for i, VD_element in enumerate(VD):
    BetaCoord = list(beta_range).index(Beta_grid[i])
    HLRCoord = list(sizes_range).index(Sizes_grid[i])
    SersicCoord = list(sersic_index_range).index(SersicIndex_grid[i])
    SMCoord = list(stellar_mass_range).index(StellarMass_grid[i])
    theGrid[BetaCoord, HLRCoord, SersicCoord, SMCoord] = VD_element

np.save("Grid.npy", theGrid)

print("Completed!")