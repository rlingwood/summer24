import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from astropy.wcs import WCS
from astropy.io import fits
from specutils import Spectrum1D
from scipy.optimize import curve_fit

path = Path("C:/Users/rubyl/Downloads/Ruby Lingwood/Data/ga20061220_46_1_reduced001.fits")

file = fits.open(path)
wcs = WCS(file[0].header)
wcs = wcs.dropaxis(2)
spec = Spectrum1D.read(file, format='wcs1d-fits')

x = 35
y = 34

ra_dec = wcs.pixel_to_world(x, y)
ra = ra_dec.ra.deg
dec = ra_dec.dec.deg

def gauss(wl, centre, width, height):
    term = (wl-centre)**2/(2*width)**2
    return height*np.e**(-1*term)

a = spec.wavelength.value
b = np.nan_to_num(spec.flux)[x, y].value

p0 = np.array([8.4041*10**6, 0.0001*10**6, 4])
popt, pcov = curve_fit(gauss, a, b, p0)

fitted = gauss(a, *popt)

linewidth = (3*10**8)*((popt[1]/2)/popt[0])
delta = linewidth/np.sqrt(8*np.log(2))
deltakms = delta/1000

fig = plt.figure(figsize=(10, 6))
plt.title(f'Flux at (RA: {ra:.2f}, Dec: {dec:.2f})')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux (erg/cm2/s/A)')
plt.plot(a, fitted, lw = '2', c = 'red')
plt.plot(a, b, lw='1', c='k')
plt.text(popt[0]+1000, max(fitted), f'Line Width: {round(popt[1], 3)} Angstrom', color='red', fontsize=12)
plt.text(popt[0]+1000, max(fitted) - 0.5, f'Velocity Dispersion: {round(deltakms, 3)} km/s', color='red', fontsize=12)
plt.xlim(popt[0] - 5000, popt[0] + 5000)
plt.show()
