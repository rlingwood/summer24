import numpy as np

# planck law 

c = 3 * (10**8)
h = 6.63 * (10**(-34))
kb = 1.38 * (10**(-23))
wavelength = 853 * (10**(-6))
temp = 15
freq = c/wavelength

a = (2*h*(freq**3))/(c**2)
exp = np.e**((h*freq)/(kb*temp))
b = 1/(exp - 1)
planckflux = a*b

# mass estimate

flux = 1.409 # jansky
flux_si = flux * (10**(-26)) # W/m^2/Hz
distance = 10**6 # pc
distance_si = distance * 3.086 * (10**16) # m
beta = -1.8
kappa = 0.01*(freq/(10**12))**beta

mass = (flux_si*(distance_si**2))/(planckflux*kappa)
print('Estimated Mass = ', '{:.3e}'.format(mass), "kg")
