import math
import numpy as np
import matplotlib.pyplot as plt
import pdb


def magnitude(flux: float) -> float:
    """Converts flux to AB magnitude.

    :param flux: Source flux in erg/s/cm^2/Hz
    :type flux: float
    :return: Source AB magnitude in mag
    :rtype: float
    """
    return -2.5*math.log10(flux/3.631e-20)

def flux(magnitude: float) -> float:
    """Converts AB magnitude to flux.

    :param magnitude: Source AB magnitude in mag
    :type: float
    :return : Source flux in erg/s/cm^2/Hz
    :rtype flux: float
    """
    return (3.631e-20)*math.pow(10.0, -1.0*magnitude/2.5)

def area(major_axis: float, minor_axis: float=None):
    """Calculates area on sky.

    :param major_axis: Major axis of source in arcsec
    :type major_axis: float
    :param minor_axis: Minor axis of source in arcsec.
        Defaults to None. minor_axis == None -> spherical source
    :type minor_axis: float
    :return: Source size
    :rtype: float
    """
    if minor_axis is None:
        minor_axis = major_axis
    return math.pi*major_axis*minor_axis

def surface_brightness(magnitude: float, area: float=None):
    """Converts magnitude to surface brightness.

    :param magnitude: Source magnitude in mag
    :type magnitude: float
    :param area: Source area in arcsec^2
    :type area: float
    :return: Source AB magnitude in mag/arcsec^2
    :rtype: float
    """
    return magnitude + 2.5*math.log10(area)

def delta_surface_brightness(aperture_two: float, magnification_two: float, aperture_one: float=7.5, magnification_one:float=1.0) -> float:
    """Computes change in surface brightness moving from instrument
    with aperture_one (length unit) at magnification_one to
    instrument with aperture_two (same length unit) operating at
    magnification_two.

    :param aperture_two: Diameter of 2nd instrument
    :type aperture_two: float
    :param magnification_two: Operating magnification
        of 2nd instrument
    :type magnification_two: float
    :param aperture_one: Diameter of 1st instrument
    :type aperture_one: float
    :param magnification_one: Operating magnification
        of 1st instrument
    :type magnification_one: float
    :return: Change in surface brightness in mag/area unit^2
    :rtype: float
    """
    return -2.5*math.log10(math.pow(magnification_one*aperture_two, 2.0)/math.pow(magnification_two*aperture_one, 2.0))

def required_aperture(delta_surface_brightness: float, magnification_two: float, aperture_one: float=7.5, magnification_one:float=1.0) -> float:
    """Computes required aperture ratio for desired change in surface
    brightness given old & new magnifications.

    :param delta_surface_brightness: Desired change in surface brightness in mag/area unit^2
    :type delta_surface_brightness: float
    :param magnification_two: Operating magnification
        of 2nd instrument
    :type magnification_two: float
    :param aperture_one: Diameter of 1st instrument
    :type aperture_one: float
    :param magnification_one: Operating magnification
        of 1st instrument
    :type magnification_one: float
    :return: Aperture ratio
    :rtype: float
    """
    return aperture_one*(magnification_two/magnification_one)*math.pow(10.0,(delta_surface_brightness/5.0))

def contrast(source_surface_brightness: float, background_surface_brightness: float) -> float:
    source_flux_per_area = flux(source_surface_brightness)
    background_flux_per_area = flux(background_surface_brightness)
    contrast = source_flux_per_area/background_flux_per_area
    return contrast

def log_contrast(source_surface_brightness: float, background_surface_brightness: float) -> float:
    source_flux_per_area = flux(source_surface_brightness)
    background_flux_per_area = flux(background_surface_brightness)
    log_contrast = -2.5*math.log10(source_flux_per_area/background_flux_per_area)
    return log_contrast

if __name__=='__main__':
    plt.ion()

    delta_sbs = np.arange(0, 2, 0.05)
    ap_rats = np.array([required_aperture(delta_sb, 1.0, 1.0, 1.0) for delta_sb in delta_sbs])

    plt.figure(0)
    plt.plot(delta_sbs, ap_rats)
    plt.xlabel(r'$\Delta \sigma$ ($\mathrm{mag}/\mathrm{arcsec}^{2}$)')
    plt.ylabel(r'$d_{\mathrm{new}}/d_{\mathrm{old}}$')

    plt.figure(1)
    src_mags = np.arange(-4.0, 20.1, 0.1)
    src_area = area(178*60, 63*60) # Angular size of M31
    src_sbs = surface_brightness(src_mags, src_area)
    sky_sbs = np.arange(16.0, 22.5, 0.5)
    for sky_sb in sky_sbs:
        contrasts = np.array([contrast(src_sb, sky_sb) for src_sb in src_sbs])
        plt.plot(src_mags, contrasts, label=r'$\sigma_{\mathrm{sky}} = %f$ ($\mathrm{mag}/\mathrm{arcsec}^{2}$)'%(sky_sb))
    plt.xlabel(r'$m_{\mathrm{M31}}$ ($\mathrm{mag}$)')
    plt.ylabel(r'$c_{\mathrm{M31}}$')
    plt.legend()

    zones = [18.0, 18.5, 19.25, 20.3, 20.8, 21.3, 21.6, 21.75, 22.0]

    plt.figure(2)
    src_mag = 3.4 # Magnitude of M31
    src_area = area(178*60, 63*60) # Angular size of M31
    src_sb = surface_brightness(src_mag, src_area)
    sky_sbs = np.arange(16.0, 22.0, 0.05)
    contrasts = np.array([contrast(src_sb, sky_sb) for sky_sb in sky_sbs])
    plt.plot(sky_sbs, contrasts)
    plt.vlines(zones, ymin = 0, ymax=max(contrasts), color='black')
    plt.xlabel(r'$\sigma_{\mathrm{sky}}$ ($\mathrm{mag}/\mathrm{arcsec}^{2}$)')
    plt.ylabel(r'$c_{\mathrm{M31}}$')

    plt.show()

    pdb.set_trace()