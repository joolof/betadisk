# betadisk

A python class to produce images of debris disks, as a function of different grain sizes. It works for scattered light images as well as for thermal emission (with some crude-ish approximation for the dust temperature).

Radiation pressure increases the eccentricity of the particles, and the increase depends on their sizes. The idea here is to consider a size distribution, and instead of creating one single image, the size distribution is binned in `ng` intervals, and `ng` images are created, to capture the different spatial scales of the debris disk. Afterwards, you can for instance do all sorts of linear combination to find the images that best reproduce your observations.

The main motivation of the approach is to get rid of as many assumptions regarding the dust properties as possible. Therefore, instead of working directly with the "true" size of the particles, the code directly uses the $\beta$ values (the ratio between radiation pressure and gravitational forces). Regardless of the luminosity or mass of the star, you should always have a distribution of $\beta$ values, ranging between 0 and 0.5 (one could argue that this might not be the case for low mass stars, but their could be stellar winds, etc). Then, you can a posteriori try to relate the $\beta$ values with the actual grain sizes. In the end, you can produce images as follows, where $\beta$ increases from the top left down to the bottom right.

![pretty](screenshots/pretty.png)

As $\beta$ increases (smaller grain sizes), the disk becomes more and more radially extended. For the last image, the orbits are so spread out because of the high eccentricity ($e \sim 0.95$) that we start to see some "noise" even in the main belt. This could be avoided by increasing the number of dust particles (though there are 10 000 000 in this example).

## Installation

Simply clone the repository, change to the directory and install using the `develop` option.

```python
python3 setup.py develop
```

The dependencies are `numpy`, `matplotlib` (though technically not necessary except for the example), and most importantly `numba` to speed things up a little. With [`numba`](https://numba.pydata.org/) you have the possibility to pre-compile parts of the code, which can significantly improve the runtime speed. In our case, since we will be launching many particles, and do very simple math operations, this is quite invaluable.

Once you cloned and installed the repository, you should therefore compile some of the code by running

```python
python3 frame.py
```

and this will create a directory called 'frame.cpython-38-x86_64-linux-gnu.so' (name will vary depending on your operating system). And that's it, you should be able to succesfully run the following

```python
python3 betadisk.py
```

## Running your first model

The code is pretty simple to use, you start by initializing the class as

```python
from betadisk import BetaDisk
disk = BetaDisk()
```

which has the following parameters (all optional)

```python
nx = 300           # number of pixels for the images
nl = 10_000_000    # number of particles to be launched
ng = 10            # number of grain size intervals
pixscale = 0.01226 # size of one pixel in arcsec
bmin = 0.001       # minimum value for beta
bmax = 0.49        # maximum value for beta
nb = 50            # number of bins for the phase function (see later)
slope = 0.5        # slope for the size distribution (might not be useful to change)
```

those should be enough (and mostly self-explanatory) to produce scattered light images. If you want to model ALMA observations, for thermal images you will need to provide the following parameters

```python
thermal = False    # False by default, need to switch it to True
lstar = None       # to compute the temperature we need a stellar luminosity
dpc = None         # we will also need a distance in pc
wave = None        # and we need to provide a wavelength in microns
dx = 0.            # possibility to have an offset in the x direction
dy = 0.            # and in the y direction
```

For thermal images, we need to provide the luminosity because the temperature of the dust grains is estimated by inverting the following equation ([Wyatt 2008](https://ui.adsabs.harvard.edu/abs/2008ARA&A..46..339W)):

```python
r = (278.3/Tdust)**2 * Lstar**0.5
```

and `r` will be in au. Since all the distances in the code are in arcseconds we also need the distance in pc for the conversion. 

Afterwards, you can call the main method

```python
disk.compute_model(a = 1.5)
```

which the following parameters

```python
a = 1.5        # the reference radius of the disk, in arcsec
dr = 0.25      # the standard deviation for the width of the main belt (normal profile)
incl = 45.     # inclination of the disk, in degrees
pa = -110.     # position angle of the disk, in degrees
opang = 0.05   # opening angle of the disk
pfunc = np.ones(nb)  # array containing the phase function
```

### A word on the phase function

The "traditional" way to deal with the phase function, either in scattered or polarized light observations, would be to use the Heynyey-Greenstein approximation, which is parametrized by a single value $g$. However, this approximation does not always work very well, and one often has to increase the complexity of the problem, for instance using a combination of two phase functions. In [Olofsson et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020A%26A...640A..12O/abstract) we proposed another approach to derive the "best" phase function directly from the observations. This is a two step process in which we run a first model wthout a phase function, accounting only for the dust density distribution, and then the phase function is derived from the brightness profile of both this first model and the observations. Afterwards, a second model can be computed using the inferred phase function.

The result of this approach is that we have two arrays, one for the scattering angle, and one for the phase function. In `compute_model` if one does not pass the parameter `pfunc`, then the array is filled with `nb` 1. In the other case, the user can provide the array `pfunc` to the method so that it is used when computing the images. In both cases, the array for the scattering angle is automatatically creating based on the length of the array for the phase function.

At the moment I do not plan on modifying the code to allow for HG approximations, but this may come in the future if I would need it. At the moment, you could still sample your HG approximation at dsicrete values and pass it to the method with `pfunc` but this is a bit brute-force because there will be some interpolation down the road.

### A minimal working example

In the `example` method you can find a minimal working example, quite similar to

```python
nx, ng = 1_000, 12
disk = BetaDisk(nx = nx, ng = ng)
disk.compute_model(a = 1.5, dr = 0.030, incl = 40.0, pa = -120.0, opang = 0.05)
print(np.shape(disk.model))
# Should return something like (12, 1000, 1000)
```

Note that I did not really notice a decrease in performance when increasing the size of the images. The position of all the particles has to be computed in any case, there is only a check towards the end to see whether there are in the field of view or not. So it should be only slightly slower if you have `nx = 200` or `nx = 1024`. This might be useful to model ALMA osbervations and want to sample your Fourier transform with high enough spatial resolution.


