import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm

#########################################
# Code to create a column density and magnetic field distribution for an hourglass core
# Following Myers et al. 2020 ApJ 896 163
# Assuming an oblate spheroid with symmetry axis parallel to magnetic field axis (z)

#########################################
# Set size of array
# odd-numbered so that it has a central pixel

zco, xco = np.mgrid[:3001, :3001]

zco -= 1501
xco -= 1501

#########################################
# Set free parameters of model

x_0 = 5
z_0 = 4
inclination = np.pi/4 # inclination; max valid value = pi/2 radians; 0 for no inclination
alpha = 5 # twist scale height; zero for no twist
A = np.cos(inclination)/np.sqrt((z_0/x_0)**2-(np.sin(inclination))**2)
r_0 = 200 # characteristic size

N_u = 3*(10**7) # column density of uniform background
DeltaN_max = 3*(10**10) # central column density 


#########################################
# column density

Y = np.pi/2 # because they say so in the paper

n_u = N_u/(2.*Y) # background volume density

n_0 = DeltaN_max/(np.pi*r_0*A) # central volume density in terms of column density

nu_0 = n_0/n_u # centre-to-edge VOLUME density contrast

xi = xco/r_0 # normalised x coordinate
zeta = zco/r_0 # normalised z coordinate

twistterm = 1/np.tanh(zeta/alpha)

# change in column density
DeltaNtop = np.pi*n_0*r_0*A
extraterm = 1 + (A**2-1)*(np.sin(inclination))**2
DeltaNbottom = 1+extraterm*(xi/A)**2+zeta**2
DeltaN = DeltaNtop/np.sqrt(DeltaNbottom)

# background column density + column density of source
N = N_u + DeltaN

#########################################
# magnetic field

s = xi/zeta

omega = np.sqrt(((xi/A)*twistterm*(1/np.cos(inclination)))**2+((zeta*(1/np.cos(inclination)))**2))

t = (1 + nu_0/(1.+(omega**2.)))/(1. + 3.*nu_0*(omega**(-2.))*(1.-(1./omega)*np.arctan(omega)))

q = A*(1/np.cos(inclination))*twistterm

# this gives the magnetic field relative to the symmetry axis, i.e. relative to image North
theta_B = np.arctan((1.-t)/((1./s) + s*t*(q**(-2.)))) 

######################
# plotting

# select which vectors to plot
use = np.where(np.logical_and(xco % 100 == 0, zco % 100 == 0))

fig = plt.figure()

ax = plt.subplot(111)

plt.imshow(DeltaN, origin='lower', extent=[np.min(xco), np.max(xco), np.min(zco), np.max(zco)], cmap=cm.YlGnBu)

plt.quiver(xco[use], zco[use],
           np.cos(np.pi/2 - theta_B[use]), # pi/2 - theta_B to make it correct angles relative to image West (x axis)
           np.sin(np.pi/2 - theta_B[use]),
           pivot='middle', scale=35, headwidth=1,
           headlength=1, headaxislength=1,
           width=0.003, color='k', #color='red'
           zorder=10)

plt.show()
