# Jeans modelling of galaxy velocity dispersion
## Dr Christopher Marsden

### What is it?

Please see: https://arxiv.org/abs/2112.09720. 

Galaxies are important in modern astrophysics, and galaxy velocity dispersion (a measure of the statistical variance in the velocities of the stellar population), is especially important, especially in connection with black hole mass. See the introduction of my paper, linked above. Computing this property is quite difficult and potentially computationally expensive, so we require a more optimized solution. My code implements the solution to the Jeans equations, modelling the gravitational dynamics of the galaxy, with (uniquely) a variable (or non-variable) IMF. This code is implemented in C++, using adaptive Richardson extrapolation to compute the integrals. The code is parallelized for multiple galaxies using OpenMP, so should be quite fast. I validated this code using velocity dispersion profiles, probability distributions and large scale statistical distributions of several astronomical surveys, including SDSS and MaNGA (see paper), and the results are (in a statistical sense) very accurate, representing ‘typical’ galaxies. On individual galaxies with less than idealized profiles, your milage may vary. I used this code to perform MCMC style modelling and determine how velocity dispersion (and it’s relationship with black hole mass) varies over the lifetime of the universe.

TLDR: This code simulates a statistically derived galaxy property called velocity dispersion using C++, for outcomes such as MCMC modelling. 

### How do I use it?

#### Good news, it has a python wrapper.

Bad news, you still have to build it. 

Requirements (C++):

* OpenMP
* Boost
* cmake

Requirements (python)

* numpy
* ctypes

I also suggest Colossus for calculating your cosmology. You can use Colossus' Fundamental Parameters methods to generate the halo rho and rs from halo mass.

#### Unix

Clone the repository to the desired location, and run:

`cmake CMakeLists.txt`

Assuming this works (if it doesn’t, check your requirements are installed – cmake may also struggle to find them, so edit CMakeLists.txt as appropriate). Then run:

`make`

This should generate the .so file, so you can use the python wrapper. The python wrapper is located in SigmaLib.py, and I refer you to the docstring of the main function, Sigma().

#### Windows
Good luck. You’re on your own. 


