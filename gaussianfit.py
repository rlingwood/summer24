import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.visualization import (SqrtStretch,
                                   ImageNormalize)

from pathlib import Path

from scipy.optimize import curve_fit

#################################################################


datadir = Path("C:/Users/rubyl/Downloads/Ruby Lingwood/Data")

outdir = Path("C:/Users/rubyl/Downloads/Ruby Lingwood/Data")

pixsize = 4.0


################################################################

# Start plotting

dset = datadir / 'afgl961_26nov_cat_mas.FIT'

# Get map

iext = datadir / 'afgl961_26nov_iext_mJysqa.fits'
mapfile = fits.open(iext)
im = np.squeeze(mapfile[0].data) # drop shallow 3rd axis from data map
var = np.squeeze(mapfile[1].data) # drop shallow 3rd axis from variance map

# Crop to central part

crop = Cutout2D(im, (round(0.5*im.shape[0]), round(0.5*im.shape[1])),
                (0.2*im.shape[1],0.2*im.shape[0]))
    
cropim = crop.data

x = np.linspace(0, cropim.shape[1], cropim.shape[1])
y = np.linspace(0, cropim.shape[0], cropim.shape[0])
x, y = np.meshgrid(x, y)

# Plot results

fig = plt.figure(figsize=(4.5,4.5))
#ax = plt.subplot(projection=cropwcs)
ax = plt.subplot()

cropim = np.nan_to_num(cropim)

norm = ImageNormalize(cropim, vmin=0., vmax=np.nanmax(cropim),
                      stretch=SqrtStretch())    
    
image = plt.imshow(cropim, origin='lower', cmap=cm.Greys, norm=norm)
        
plt.xlabel('Right Ascension (J2000)')
plt.ylabel('Declination (J2000)')

def gaussian(x, y, cx, cy, height, ratio):
    n = 0.5*np.sqrt(2)
    G = np.exp(-(x-cx)**2/(n*ratio*height**2) - (y-cy)**2/(n*height**2))
    return np.ravel(G)

cropim[cropim < 5*np.mean(cropim)] = 0
p0 = np.array([20, 30, 1, 1])
popt, pcov = curve_fit(lambda X, cx, cy, height, ratio: gaussian(X[0], X[1], cx, cy, height, ratio), (x, y), np.ravel(cropim), p0)
fitted = gaussian(x, y, *popt)
fitted = fitted.reshape(cropim.shape)
plt.contour(x, y, fitted, 0, colors='red')
print('Aspect Ratio:', round(popt[3],2))

cbar = plt.colorbar(image, fraction=0.046, pad=0.04)
cbar.set_label(r'Surface Brightness (mJy arcsec$^{-2}$)', size=12)
cbar.ax.tick_params(labelsize=11)
