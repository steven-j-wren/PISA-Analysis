
import matplotlib.pyplot as plt
import numpy as np

from scipy.special import erf, erfinv
from scipy.stats import norm, spearmanr
from scipy.optimize import curve_fit

def pdf(x):
    return 1/np.sqrt(2*np.pi) * np.exp(-x**2/2)

def cdf(x):
    return (1 + erf(x/np.sqrt(2))) / 2

def skew(x,e=0,s=1,w=1,a=0):
    t = (x-e) * w
    return 2 * s * pdf(t) * cdf(a*t)

def crystalball(x,s,a,n,xbar,sigma):
    A = np.power(n/np.absolute(a),n)*np.exp(-np.absolute(a)**2/2)
    B = n/np.absolute(a)-np.absolute(a)
    C = n/np.absolute(a)*(1/(n-1))*np.exp(-np.absolute(a)**2/2)
    D = np.sqrt(np.pi/2)*(1+erf(np.absolute(a)/np.sqrt(2)))
    N = 1/(sigma*(C+D))
    return np.piecewise(x, [(x-xbar)/sigma > -a], [lambda x:s*N*np.exp(-np.power(x-xbar,2)/(2*sigma**2)), lambda x:s*N*A*np.power((B-(x-xbar)/sigma),-n)])

fine_binning = np.linspace(-2.5,1.0,4000)

plt.plot(fine_binning, skew(fine_binning, 0.5, 600.0, 2.5, -4.5))
plt.savefig('lol.png')

plt.close()

plt.plot(fine_binning, crystalball(fine_binning, 600.0, 1.0, 1.1, 0.1, 0.2))
plt.savefig('lol2.png')
