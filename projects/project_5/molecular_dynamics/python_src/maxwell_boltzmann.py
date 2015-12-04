import numpy as np
import matplotlib.pyplot as plt

def maxwell_boltzmann(a, x):
    return np.sqrt(2 / np.pi) * x**2 * np.exp(-x**2 / (2.0*a**2)) / a**3

x_values = np.linspace(0, 20, 1000)
a_values = range(1, 5)

for a in a_values:
    m_values = maxwell_boltzmann(a, x_values)
    plt.plot(x_values, m_values, label='$ a = %d $' % a, alpha=0.7, linewidth=1.5)
plt.legend(loc='best')
plt.xlabel('$x$')
plt.ylabel('PDF')
plt.grid('on')
plt.savefig('maxwell_boltzmann.pdf')
