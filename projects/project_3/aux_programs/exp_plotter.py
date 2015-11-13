import numpy as np
import matplotlib.pyplot as plt

def exp_alpha(r, alpha):
    return np.exp(-alpha * r)

x_values = np.linspace(-8, 8, 1000)
y_values = exp_alpha(x_values, 2)

for i, y in enumerate(y_values):
    if y < 1.0e-6:
        print x_values[i]
