import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from specutils import Spectrum1D

path = Path("C:/Users/rubyl/Downloads/Ruby Lingwood/Data/ga20061220_46_1_reduced001.fits")
spec = Spectrum1D.read(path, format='wcs1d-fits')

x = 30
y = 30

plt.title('Flux at [' +str(x)+', '+str(y)+']')
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux (erg/cm2/s/A)')
plt.plot(spec.wavelength, np.nan_to_num(spec.flux)[x, y]) 
# plot flux against wavelength for the specified coordinate location (x, y)
plt.show()
