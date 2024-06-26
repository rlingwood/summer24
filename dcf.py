import numpy as np

velocitydispersion = 0.284 # km/s
angledispersion = 7.43 # degrees
mass = 6.65 * 10**35 # kg
c = 9.74 * (10**16) # m
a = c * 1.8 # m

volume = (4/3) * np.pi * (a**2) * c # m^3

volumecm = volume * 10**6 # cm^3

hydrogenmass = 3.33 * (10**(-27)) # kg

numberdensity = mass / (hydrogenmass*volumecm) # 1/cm^3

dcf = 9.3 * ((np.sqrt(numberdensity)*velocitydispersion)/angledispersion)
print('Magnetic Field Strength =', round(dcf, 2), 'microgauss')
