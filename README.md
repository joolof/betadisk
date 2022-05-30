# betadisk

A python class to produce images of debris disks, as a function of different grain sizes. It works for scattered light images as well as for thermal emission (with some crude-ish approximation for the dust temperature).

Radiation pressure increases the eccentricity of the particles, and the increase depends on their sizes. The idea here is to consider a size distribution, and instead of creating one single image, the size distribution is binned in `ng` intervals, and `ng` images are created, to capture the different spatial scales of the debris disk. Afterwards, you can for instance do all sorts of linear combination to find the images that best reproduce your observations.

The main motivation of the approach is to get rid of as many assumptions regarding the dust properties as possible. Therefore, instead of working directly with the "true" size of the particles, the code directly uses the $beta$ values (the ratio between radiation pressure and gravitational forces).

![pretty](screenshots/pretty.png)


