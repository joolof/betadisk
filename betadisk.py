import sys
import time
import frame
import numpy as np
import matplotlib.pyplot as plt
# -----------------------------------------------------------
np.seterr(over='ignore', divide='ignore')
CC = 2.9979e+14     # in microns/s 
# -----------------------------------------------------------

class BetaDisk(object):
    """
    docstring for BetaDisk
    """
    def __init__(self, nx = 300, nl = 10_000_000, nb = 50, ng = 10, pixscale = 0.01226, bmax = 0.49, \
            bmin = 0.001, slope = 0.5, thermal = False, lstar = None, dpc = None, owave = None, \
            nwav = 100, dx = 0., dy = 0.):
        """
        description
        """
        self._nx, self._nl, self._nb, self._ng, self._nwav = nx, nl, nb, ng, nwav
        self._pixscale = pixscale
        self._slope = slope
        self._dx, self._dy = dx, dy
        self._beta = self._get_size(bmin, bmax)
        self._bgrid = np.linspace(bmin, bmax, num = ng+1)
        self.model = np.zeros((self._ng, self._nx, self._nx))
        self._a, self._dr, self._incl, self._pa, self._opang, self._pfunc = 1., 0.02, 0., 0., 0.05, None
        self._dpi, self._is_hg, self._ghg = False, False, 0.
        self._thermal, self._lstar, self._dpc, self._owave = thermal, lstar, dpc, owave
        self.wave = np.geomspace(0.1, 2000., num = self._nwav)
        self.sed = np.zeros((self._nwav, self._ng))

    def _get_size(self, smin, smax):
        """
        Draw from a power-law distribution
        with a slope in "slope"
        """
        np.random.seed(10)
        slope = self._slope + 1.
        return (((smax**slope)-(smin**slope)) * np.random.random(size=self._nl) + smin**slope)**(1./slope)

    def compute_sed(self, **kwargs):
        """
        Compute a model
        """
        self._check_parameters(kwargs)
        if self._lstar is None:
            raise ValueError('You need to provide a luminosity for the SED')
        if self._dpc is None:
            raise ValueError('You need to provide a distance in pc for the SED')
        for i in range(self._ng):
            sel = ((self._beta >= self._bgrid[i]) & (self._beta < self._bgrid[i+1]))
            btmp = self._beta[sel]
            self.sed[:,i] = frame.sed(btmp, self._a, self._dr, self._opang, self._slope,\
                                             self._dpc, self._lstar, self.wave, len(btmp))

    def compute_model(self, **kwargs):
        """
        Compute a model
        """
        self._check_parameters(kwargs)
        if self._pfunc is None:
            self._theta = np.linspace(0., np.pi, num = self._nb)
            self._pfunc = np.ones(self._nb)
        """
        Compute the model
        """
        if self._thermal:
            if self._lstar is None:
                raise ValueError('You need to provide a luminosity for thermal images')
            if self._dpc is None:
                raise ValueError('You need to provide a distance in pc for thermal images')
            if self._owave is None:
                raise ValueError('You need to provide a wavelength in microns for thermal images')
            for i in range(self._ng):
                sel = ((self._beta >= self._bgrid[i]) & (self._beta < self._bgrid[i+1]))
                btmp = self._beta[sel]
                self.model[i,:,:] = frame.alma(btmp, self._a, self._dr, self._incl, \
                                                 self._opang, self._pa, self._pixscale, self._slope,\
                                                 self._dpc, self._lstar, self._owave, \
                                                 self._dx, self._dy, len(btmp), self._nx)
        else:
            for i in range(self._ng):
                if type(self._ghg) is np.ndarray:
                    ghg = self._ghg[i]
                else:
                    ghg = self._ghg
                if self._pfunc is not None:
                    if len(np.shape(self._pfunc)) == 2:
                        pfunc = self._pfunc[:, i]
                    else:
                        pfunc = self._pfunc
                sel = ((self._beta >= self._bgrid[i]) & (self._beta < self._bgrid[i+1]))
                btmp = self._beta[sel]
                self.model[i,:,:] = frame.sphere(btmp, self._a, self._dr, self._incl, \
                                                 self._opang, self._pa, self._pixscale, self._slope,\
                                                 self._is_hg, ghg, self._theta, pfunc, self._dx, self._dy,\
                                                 len(btmp), self._nx, self._dpi)

    """
    Check the parameters that are passed as kwargs
    """
    def _check_parameters(self, kwargs):
        """
        Check parameters that are being passed
        """
        if 'a' in kwargs:
            self._a = kwargs['a']
        if 'dr' in kwargs:
            self._dr = kwargs['dr']
        if 'incl' in kwargs:
            self._incl = kwargs['incl'] * np.pi / 180.
        if 'pa' in kwargs:
            self._pa = kwargs['pa'] * np.pi / 180.
        if 'opang' in kwargs:
            self._opang = kwargs['opang']
        if 'dpi' in kwargs:
            self._dpi = kwargs['dpi']
        if 'is_hg' in kwargs:
            self._is_hg = kwargs['is_hg']
            if 'ghg' in kwargs:
                self._ghg = kwargs['ghg']
                if type(self._ghg) is np.ndarray:
                    if len(self._ghg) != self._ng:
                        print('The array for HG should have the same size as the number of grain sizes')
                        sys.exit(0)
        if 'pfunc' in kwargs:
            self._pfunc = kwargs['pfunc']
            self._theta = np.linspace(0., np.pi, num = np.shape(self._pfunc)[0])

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, a):
        self._a = a

    @property
    def opang(self):
        return self._opang

    @opang.setter
    def opang(self, opang):
        self._opang = opang

    @property
    def dr(self):
        return self._dr

    @dr.setter
    def dr(self, dr):
        self._dr = dr

    @property
    def incl(self):
        return self._incl * 180. / np.pi

    @incl.setter
    def incl(self, incl):
        self._incl = incl * np.pi/180.

    @property
    def pa(self):
        return self._pa * 180. / np.pi

    @pa.setter
    def pa(self, pa):
        self._pa = pa * np.pi/180.

def example_image():
    nx, ng = 1_000, 12
    if nx%2 ==0:
        cx = nx//2 - 0.5
    else:
        cx = nx//2
    disk = BetaDisk(nx = nx, ng = ng)
    t0 = time.perf_counter()
    disk.compute_model(a = 1.5, dr = 0.030, incl = 40.0, pa = -120.0, opang = 0.05)
    print('Whole process took: {:.2f} sec'.format(time.perf_counter()-t0))
    """
    Make a pretty plot
    """
    xlim = cx * disk._pixscale
    fig = plt.figure(figsize=(10, 10*3./4.))
    ax = plt.GridSpec(3, 4, hspace=0.0, wspace=0.0)
    ax.update(left=0.0,right=1.,top=1.,bottom=0.0,wspace=0.0,hspace=0.00)
    ct = 0
    for i in range(3):
        for j in range(4):
            ax1 = fig.add_subplot(ax[i,j])
            ax1.imshow(disk.model[ct,], origin='lower', extent = [xlim, -xlim, -xlim, xlim],
                            vmax=np.percentile(disk.model[ct,], 99.5), cmap = 'inferno')
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_axis_off()
            ct+=1
    plt.savefig('screenshots/pretty.png')
    plt.show()

def example_sed():
    ng = 12
    disk = BetaDisk(ng = ng, lstar = 2.5, dpc = 100.)
    t0 = time.perf_counter()
    disk.compute_sed(a = 0.5, dr = 0.030, opang = 0.05)
    print('Whole process took: {:.2f} sec'.format(time.perf_counter()-t0))
    sedisk = np.sum(disk.sed, axis = 1) * CC/disk.wave
    index = np.log10(sedisk[90]/sedisk[80])/np.log10(disk.wave[90]/disk.wave[80])
    print('Mm index between {:.2f} and {:.2f} is: {:.2f}'.format(disk.wave[80], disk.wave[90], index))
    """
    Make a pretty plot
    """
    fig = plt.figure(figsize=(7,6))
    ax1 = fig.add_axes([0.16, 0.14, 0.8, 0.79])
    ax1.set_xscale("log", nonpositive='clip')
    ax1.set_yscale("log", nonpositive='clip')
    for i in range(ng):
        ax1.plot(disk.wave, disk.sed[:,i] * CC/disk.wave, lw = 1)
    ax1.plot(disk.wave, sedisk, 'k', lw = 3)
    ax1.set_xlim(0.1, 2000.)
    ax1.set_xlabel('$\lambda$ [$\mu$m]')
    ax1.set_ylabel('$\\nu F_{\\nu}$ [erg/s/cm/cm]')
    plt.show()



if __name__ == "__main__":
    # example_image()
    example_sed()


