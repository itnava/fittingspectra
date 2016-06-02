import numpy as np

def onegauss(x, a, b, c):
    return a*np.exp(-(x - b)**2 / (2. * c**2))

def multigauss(x, a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4, d):
    return onegauss(x, a1, b1, c1) + onegauss(x, a2, b2, c2) + onegauss(x, a3, b3, c3) + onegauss(x, a4, b4, c4) + d
