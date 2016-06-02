import numpy as np

def singlegauss(x, a, b, c, d):
    return a*np.exp(-(x - b)**2 / (2. * c**2)) + d

