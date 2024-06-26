import numpy as np
import math

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.nddata.utils import Cutout2D
from astropy.visualization import (SqrtStretch,
                                   ImageNormalize)
import copy

from pathlib import Path
from scipy.optimize import curve_fit

#################################################################


datadir = Path("C:/Users/rubyl/Downloads/Ruby Lingwood/Data")

outdir = Path("C:/Users/rubyl/Downloads/Ruby Lingwood/Data")

pixsize = 4.0


################################################################

# Start plotting

dset = datadir / 'afgl961_26nov_cat_mas.FIT'

# Get the vectors

hdulist_vec = fits.open(dset)

vec = hdulist_vec[1].data
vec_ra = vec.field('RA')
vec_dec = vec.field('DEC')
vec_x = vec.field('X')
vec_y = vec.field('Y')
vec_i = np.asarray(vec.field('I'))
vec_di = vec.field('DI')
vec_p = np.asarray(vec.field('P'))/100.
vec_dp = np.asarray(vec.field('DP'))/100.
vec_ang = vec.field('ANG') # Gets you polarization angles
vec_dang = vec.field('DANG')
vec_q = vec.field('Q')
vec_dq = vec.field('DQ')
vec_u = vec.field('U')
vec_du = vec.field('DU')


# Get map

iext = datadir / 'afgl961_26nov_iext_mJysqa.fits'
mapfile = fits.open(iext)
im = np.squeeze(mapfile[0].data) # drop shallow 3rd axis from data map
var = np.squeeze(mapfile[1].data) # drop shallow 3rd axis from variance map

# get World Coordinate System info
wcs = WCS(mapfile[0].header)
wcs = wcs.dropaxis(2) # drop shallow 3rd axis


# Crop to central part

crop = Cutout2D(im, (round(0.5*im.shape[0]), round(0.5*im.shape[1])),
                (0.2*im.shape[1],0.2*im.shape[0]), wcs=wcs)
    
cropwcs = crop.wcs
cropim = crop.data

# Get a more specific set of vectors

xco, yco = cropwcs.wcs_world2pix((360.+np.degrees(vec_ra))*u.degree,
                                 np.degrees(vec_dec)*u.degree, 0)

good = np.where(np.logical_and(np.logical_and(vec_p/vec_dp > 3.,
                                              vec_i/vec_di > 5.),
                               np.logical_and(np.logical_and(xco > 0.1*cropim.shape[0], xco < 0.9*cropim.shape[0]),
                                              np.logical_and(yco > 0.1*cropim.shape[0], yco < 0.9*cropim.shape[1]))))


vxco = copy.deepcopy(xco[good])
vyco = copy.deepcopy(yco[good])
vang = copy.deepcopy(vec_ang[good])
vdang = copy.deepcopy(vec_dang[good])
vp = copy.deepcopy(vec_p[good])
vi = copy.deepcopy(vec_i[good])

# Plot results

fig = plt.figure(figsize=(4.5,4.5))
ax = plt.subplot(projection=cropwcs)

norm = ImageNormalize(cropim, vmin=0., vmax=np.nanmax(cropim),
                      stretch=SqrtStretch())    
    
im = plt.imshow(cropim, origin='lower', cmap=cm.Greys, norm=norm)
        
plt.xlabel('Right Ascension (J2000)')
plt.ylabel('Declination (J2000)')
    
beam = plt.Circle((5.,5.), 7.05/pixsize, color='k')
    
ax.add_artist(beam)

# Plot POL-2 vectors

# note that polarization angles are measured E of N and quiver measures angles anticlockwise from image left
# so plotting polarization angles using quiver makes them appear as magnetic field angles on the map
plt.quiver(vxco, vyco,
           np.cos(np.radians(vang)),
           np.sin(np.radians(vang)),
           pivot='middle', scale=25, headwidth=1.,
           headlength=0., headaxislength=1., width=0.005, color='red', label='POL-2', zorder=12)

# Colour bar
        
cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
cbar.set_label(r'Surface Brightness (mJy arcsec$^{-2}$)', size=12)
cbar.ax.tick_params(labelsize=11)

# get vector pairs

# pairs = np.zeros(13861)

vectors = np.column_stack((vxco,vyco,vang,vdang))

def pair(array):
    length = array.shape[0]
    width = array.shape[1]
    pairs = np.zeros(((math.comb(length, 2)), 2*width))
    counter = 1
    j = 0
    k = 1
    jmax = length - counter

    while counter < (length): 
        while j < jmax:
            pairs[j,0] = array[(counter-1), 0]
            pairs[j,1] = array[(counter-1), 1]
            pairs[j,2] = array[(counter-1), 2]
            pairs[j,3] = array[(counter-1), 3]
            pairs[j,4] = array[k, 0]
            pairs[j,5] = array[k, 1]
            pairs[j,6] = array[k, 2]
            pairs[j,7] = array[k, 3]
            j = j + 1
            k = k + 1
        counter = counter + 1
        k = counter
        jmax = jmax+(length - counter)
    return pairs

vectorpairs = pair(vectors)

def distance(ax, ay, bx, by):
    dist = np.sqrt((bx-ax)**2 + (by-ay)**2)
    return dist

def dispersion(atheta, btheta):
    disp = abs(atheta-btheta)
    return disp

def uncertainty(adtheta, bdtheta):
    return np.sqrt((adtheta**2)+(bdtheta**2))

length = vectorpairs.shape[0]
data = np.zeros((length, 4))
i = 0

while i < length:
    data[i, 0] = i
    data[i, 1] = distance(vectorpairs[i, 0], vectorpairs[i, 1], vectorpairs[i, 4], vectorpairs[i, 5])
    data[i, 2] = dispersion(vectorpairs[i, 2], vectorpairs[i, 6])
    data[i, 3] = uncertainty(vectorpairs[i, 3], vectorpairs[i, 7])
    data[i, 2] = data[i, 2] - data[i, 3]
    i = i + 1

llength = 15
dispersionaverages = np.zeros((llength))
j = 1
while j < (llength+1):
    filtereddata = data[(data[:, 1] > j) & (data[:, 1] < (j+1))]
    average = np.mean(filtereddata[:, 2])
    dispersionaverages[j-1] = average
    j = j + 1
    
l = np.arange(llength)

def structure(l, b, m):
    return b**2 + (m**2)*(l**2)

popt, pcov = curve_fit(structure, l, dispersionaverages**2)
fitted = structure(l, popt[0], popt[1])

fig2 = plt.figure(figsize=(10,6))
plt.scatter(l, dispersionaverages, c='k')
plt.plot(l, np.sqrt(fitted))
plt.xlabel('l')
plt.ylabel('Dispersion')
print('b = ', round(popt[0]/np.sqrt(2),2))
