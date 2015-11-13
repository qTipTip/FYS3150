"""
Plots 5 Laguerrue polynomials
and 5 Lagrange polynomials
"""

import numpy as np
import matplotlib.pyplot as plt

@np.vectorize
def legendre(n, x):
    """
    Returns the legendre polynomial of degree n evaluated at x
    """
    r = 0.
    s = 1.
    t = 0.
    for i in range(n):
        t = r
        r = s
        s = (2*i + 1)*x*r - i*t
        s = s / (i + 1)
    return s

@np.vectorize
def laguerre(n, x):
    """
    Returns the laguerre polynomial of degree n evaluated at x
    """
    L_prev = 1
    L_curr = 1 - x
    L_new = 0
    for i in range(n):
        L_new = ((2*i + 1 - x)*L_curr - i*L_prev) / float(i + 1)
        L_prev = L_curr 
        L_curr = L_new
    return L_curr

x_values = np.linspace(-1, 1, 1000)

for n in range(5):
    y_values = legendre(n, x_values)
    plt.plot(x_values, y_values, label='$L_%d$' % n)
plt.grid('on')
plt.legend()
plt.savefig('../plots/legendre_polynomials.pdf')
plt.clf()

x_values = np.linspace(-5, 15, 1000)

for n in range(5):
    y_values = laguerre(n, x_values)
    plt.plot(x_values, y_values, label='$L_%d$' % n)
plt.ylim((-10, 20))
plt.grid('on')
plt.legend()
plt.savefig('../plots/laguerre_polynomials.pdf')

