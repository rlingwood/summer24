#!/usr/bin/python3

import sys

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams
from matplotlib.patches import Rectangle, Ellipse
from matplotlib.lines import Line2D

#########################################
# Code to create a column density and magnetic field distribution for an hourglass core
# Following Myers et al. 2020 ApJ 896 163
# Assuming an oblate spheroid with symmetry axis parallel to magnetic field axis (z)
# inclination angle = 0, and no toroidal twisting of the magnetic field

#########################################
# Set size of array
# odd-numbered so that it has a central pixel

zco, xco = np.mgrid[:20001, :20001]

zco -= 10001
xco -= 10001

#########################################
# Set free parameters of model

A = 2 # aspect ratio

r_0 = 200 # characteristic size

N_u = 1 # column density of uniform background
DeltaN_max = N_u*(10.**3) # central column density = N_u + DeltaN_max

#########################################
# column density

Y = np.pi/2 # because they say so in the paper

n_u = N_u/(2.*Y) # background volume density

n_0 = DeltaN_max/(np.pi*r_0*A) # central volume density in terms of column density

nu_0 = n_0/n_u # centre-to-edge VOLUME density contrast

xi = xco/r_0 # normalised x coordinate
zeta = zco/r_0 # normalised z coordinate

# change in column density
DeltaN = (np.pi*n_0*r_0*A)/np.sqrt(1+(xi/A)**2.+(zeta**2.))

# background column density + column density of source
N = N_u + DeltaN

#########################################
# magnetic field

s = xi/zeta

omega = np.sqrt((xi/A)**2.+(zeta**2.))

t = (1 + nu_0/(1.+(omega**2.)))/(1. + 3.*nu_0*(omega**(-2.))*(1.-(1./omega)*np.arctan(omega)))

# this gives the magnetic field relative to the symmetry axis, i.e. relative to image North
theta_B = np.arctan((1.-t)/((1./s) + s*t*(A**(-2.)))) 

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
           pivot='middle', scale=100, headwidth=1.,
           headlength=0., headaxislength=1.,
           width=0.002, color='k', #color='red'
           zorder=10)

plt.show()
