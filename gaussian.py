import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# finding the aspect ratio of a 2d gaussian

gridsize = 200
x = np.linspace(0, gridsize, gridsize)
y = np.linspace(0, gridsize, gridsize)
x, y = np.meshgrid(x, y)

def gaussian(x, y, cx, cy, height, ratio):
    n = 0.5*np.sqrt(2)
    G = np.exp(-(x-cx)**2/(n*ratio*height**2) - (y-cy)**2/(n*height**2))
    return np.ravel(G)

data = gaussian(x, y, 80, 50, 100, 2)
data = data.reshape(gridsize, gridsize)

SNR = 5
sigma = data.std() / np.sqrt(SNR)
noise = np.random.normal(0, scale=sigma, size=(gridsize, gridsize))
noisy = data + noise

popt, pcov = curve_fit(lambda X, cx, cy, height, ratio: gaussian(X[0], X[1], cx, cy, height, ratio), (x, y), np.ravel(noisy))
fitted = gaussian(x, y, *popt)
fitted = fitted.reshape(200, 200)

fig = plt.figure(figsize=(8, 8))
plt.imshow(noisy)
plt.contour(x, y, fitted, 0, colors='w')
print('Aspect Ratio:', round(popt[3],2))
