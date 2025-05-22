import numpy as np
from numba import njit
from numba.pycc import CC as compiler
# -----------------------------------------------------------
ccompile = compiler('frame')
# -----------------------------------------------------------
CC = 2.9979e+14     # in microns/s 
HH = 6.6262e-27
KK = 1.3807e-16
# -----------------------------------------------------------
@njit(nogil=True)
def get_true_anomaly(ecc, mean_anom):
    """
    Method to compute the eccentric anomaly for a given eccentricity and mean anomaly
    Input mean anomaly is in radians, and outpout is also in radians.
    """
    is_negative = False
    if ecc == 0.0:
        return mean_anom

    if ecc < .3:
        curr = np.arctan2(np.sin(mean_anom), np.cos(mean_anom) - ecc)
        err = curr - ecc * np.sin(curr) - mean_anom
        curr -= err / (1. - ecc * np.cos(curr))
        return curr

    if mean_anom < 0.:
        mean_anom = -mean_anom
        is_negative = True

    curr = mean_anom
    thresh = 1.e-8* np.abs( 1. - ecc);
    if  (((ecc > .8) & (mean_anom < np.pi / 3.)) | (ecc > 1.)):
        trial = mean_anom / np.abs(1. - ecc)
        if trial * trial > 6. * np.abs(1. - ecc):
            if mean_anom < np.pi:
                trial = np.exp(np.log(6. * mean_anom)/3.)
            else:
                trial = np.sinh( mean_anom / ecc) ### is this 'asinh' or 'sinh'???
        curr = trial

    if ecc < 1.:
        err = curr - ecc * np.sin(curr) - mean_anom
        while np.abs(err) > thresh:
            curr -= err / (1. - ecc * np.cos(curr))
            err = curr - ecc * np.sin( curr) - mean_anom
    else:
        err = ecc * np.sinh(curr) - curr - mean_anom
        while np.abs(err) > thresh:
            curr -= err / (ecc * np.cosh(curr) - 1.)
            err = ecc * np.sinh(curr) - curr - mean_anom

    if is_negative:
        return -curr
    else:
        return curr

@njit(nogil=True)
def get_launch(ec):
    """
    Method to randomly pick the launching points from the parent
    planetesimal belt
    """
    Ma = np.random.random() * 2. * np.pi - np.pi
    nu = get_true_anomaly(ec, Ma)
    return 2. * np.arctan2(np.sqrt(1.e0 + ec) * np.sin(nu/2.), np.sqrt(1.e0 - ec) * np.cos(nu/2.))
    # return nu

@njit(nogil=True)
def binarysearch(arr, h, x):
    m, l = 0, 0
    while (h-l >1):
        m = l + (h-l)//2
        if arr[m] == x:
            return m
        if arr[m] <x:
            l = m
        else:
            h = m
    return l

@njit(nogil=True)
def get_xyz(a, opang, ecc, beta):
    """
    Get the orbital parameters of the parent bodies
    """
    w = np.random.random() * 2. * np.pi - np.pi
    pa = np.random.random() * 2. * np.pi - np.pi
    deltai = np.random.normal(0., opang)
    # sma = np.random.normal(a, dr)
    nu_l = get_launch(ecc)
    """
    Get the updated orbital parameters
    """
    an = a * (1.e0 - beta) / (1.e0 - 2. * beta*(1.0 + ecc * np.cos(nu_l))/(1.0 - ecc * ecc))
    en = 1. / (1.0 - beta) * np.sqrt(ecc * ecc + 2.0 * beta * ecc * np.cos(nu_l) + beta * beta)
    wn = np.arctan2(beta * np.sin(nu_l), (beta * np.cos(nu_l) + ecc)) + w
    """
    Now populate each of those orbits
    """
    nu_d = get_launch(en)
    r = an * (1.0 - en * en) / (1.0 + en * np.cos(nu_d))
    x = r * np.cos(nu_d)
    y = r * np.sin(nu_d)

    dist = np.sqrt(x**2. + y**2.) # The real distance, not projected or anything
    xm = x*(np.cos(pa)*np.cos(wn) - np.sin(wn)*np.cos(deltai)*np.sin(pa)) -\
         y*(np.sin(wn)*np.cos(pa) +np.cos(wn)*np.cos(deltai)*np.sin(pa))
    ym = x*(np.sin(pa)*np.cos(wn) + np.sin(wn)*np.cos(deltai)*np.cos(pa)) +\
         y*(np.cos(wn)*np.cos(deltai)*np.cos(pa) - np.sin(wn)*np.sin(pa))
    zm = x*(np.sin(wn)*np.sin(deltai)) + y*(np.cos(wn)*np.sin(deltai))
    return xm, ym, zm, dist


@njit(nogil=True)
@ccompile.export('sed', 'f8[:](f8[:], f8[:], f8, f8, f8, f8, f8[:], i4)')
def sed(beta, a, opang, slope, dpc, lstar, wave, nl):
    nu = CC / wave
    bterm = 2. * HH * nu * nu * nu / (CC * CC)
    expcst = nu * HH / KK
    np.random.seed(10)
    ecc = 0.               # This only works for circular disks for now
    sedist = np.zeros(len(wave))
    """
    Here we go
    """
    for il in range(nl):
        _, _, _, dist = get_xyz(a[il], opang, ecc, beta[il])
        btemp = 278.3 / np.sqrt(dist*dpc) * lstar**0.25
        expterm = np.exp(expcst / btemp)
        bplanck = bterm / (expterm-1.0) # should be in erg/s/cm/cm/Hz
        corr_fac = (beta[il]/0.49)**(1.5-slope-2.)
        corr_fac *= ((1.-beta[il])/(1.-2.*beta[il]))**(1.5)
        sedist += (bplanck * corr_fac)
    return sedist

@njit(nogil=True)
@ccompile.export('alma', 'f8[:,:](f8[:], f8[:], f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, i4, i4)')
def alma(beta, a, incl, opang, dpa, pixscale, slope, dpc, lstar, wave, dx, dy, nl, nx):
    nu = CC / wave
    bterm = 2. * HH * nu * nu * nu / (CC * CC)
    expcst = nu * HH / KK
    np.random.seed(10)
    cx = nx / 2.
    ecc = 0.               # This only works for circular disks for now
    dpa = -(dpa+np.pi)
    image = np.zeros((nx, nx))
    cosi = np.cos(incl)
    sini = np.sin(incl)
    cospa = np.cos(dpa)
    sinpa = np.sin(dpa)
    """
    Here we go
    """
    for il in range(nl):
        xm, ym, zm, dist = get_xyz(a[il], opang, ecc, beta[il])
        btemp = 278.3 / np.sqrt(dist*dpc) * lstar**0.25
        expterm = np.exp(expcst / btemp)
        bplanck = bterm / (expterm-1.0) # the units don't matter since this is at one wavelength
        """
        Rotate along the x axis first for the inclination,
        and then rotate for the position angle
        """
        ym = ym * cosi - zm * sini
        xn = xm * cospa - ym * sinpa
        yn = xm * sinpa + ym * cospa

        xn = int(round((xn + dy)/pixscale + cx - 0.5)) # I put dy so that this would be in Declination
        yn = int(round((yn + dx)/pixscale + cx - 0.5)) # I put dx so that this would be in RA. And I'm using -dx so that it goes in the other direction
        # xn = int((xn + dy)/pixscale + cx)
        # yn = int((yn + dx)/pixscale + cx)
        if ((xn>=0) and (xn<=nx-1) and (yn>=0) and (yn<=nx-1)):
            corr_fac = (beta[il]/0.49)**(1.5-slope-2.)
            corr_fac *= ((1.-beta[il])/(1.-2.*beta[il]))**(1.5)
            image[xn,yn] += (bplanck * corr_fac)
    return image


@njit(nogil=True)
@ccompile.export('sphere', 'f8[:,:](f8[:], f8[:], f8, f8, f8, f8, f8, b1, f8, f8[:], f8[:], f8, f8, i4, i4, b1)')
def sphere(beta, a, incl, opang, dpa, pixscale, slope, is_hg, ghg, theta, pfunc, dx, dy, nl, nx, dpi):
    np.random.seed(10)
    cx = nx / 2.
    ecc = 0.               # This only works for circular disks for now
    dpa = -(dpa+np.pi)
    image = np.zeros((nx, nx))
    cosi = np.cos(incl)
    sini = np.sin(incl)
    cospa = np.cos(dpa)
    sinpa = np.sin(dpa)
    """
    Get the orbital parameters of the
    nl particules that will be launched
    """
    for il in range(nl):
        xm, ym, zm, dist = get_xyz(a[il], opang, ecc, beta[il])
        th = np.arccos((ym * sini - zm * cosi) / dist)
        if is_hg:
            hg = 1./(4.*np.pi) * (1. - ghg*ghg)/(1.+ghg*ghg - 2. * ghg * np.cos(th))**1.5
            if dpi:
                hg *= (1.-np.cos(th)**2.)/(1.+np.cos(th)**2.)
        else:
            it = binarysearch(theta, len(theta), th)
            hg = pfunc[it]
        """
        Rotate along the x axis first for the inclination,
        and then rotate for the position angle
        """
        ym = ym * cosi - zm * sini
        xn = xm * cospa - ym * sinpa
        yn = xm * sinpa + ym * cospa

        xn = int(round((xn + dy)/pixscale + cx - 0.5)) # I put dy so that this would be in Declination
        yn = int(round((yn + dx)/pixscale + cx - 0.5)) # I put dx so that this would be in RA. And I'm using -dx so that it goes in the other direction
        # xn = int((xn + dy)/pixscale + cx) # I put dy so that this would be in Declination
        # yn = int((yn - dx)/pixscale + cx) # I put dx so that this would be in RA. And this is -dx so that it goes in the other direction
        if ((xn>=0) and (xn<=nx-1) and (yn>=0) and (yn<=nx-1)):
            """
            It should be 1.5-slope and then -2 for the size**2,
            since s propto 1/beta there is a 1/beta**2
            """
            corr_fac = (beta[il])**(1.5 - slope - 2.)
            corr_fac *= ((1.-beta[il])/(1.-2.*beta[il]))**(1.5)
            # image[xn,yn] += (corr_fac * hg)
            image[xn,yn] += (corr_fac * hg/dist/dist)
    return image

if __name__ == "__main__":
    ccompile.compile()


