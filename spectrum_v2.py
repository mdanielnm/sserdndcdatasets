import math
import cmath
import numpy
from scipy.special import lambertw

# Spectrum of a network with delay (control region)

# parameters of Linearized system
a = 0.01
b = math.pi
# parameters of Linearized coupling
delta = 1.0
eta = 0.1


# the degree of each vertice
d = 2

# the imprimivity index number
m = 1

# define some variables
tau_min = 0.001
tau_max = 10
passo_tau = 0.001
k_min = 0.001
k_max = 3.0
passo_k = 0.002

# file to write data
f = open("spectrum_v2_oscillators_stuart_landau_non_directed_ring_wcf1.txt", "w")

# this is just auxiliary to save data, the size of r is accordingly m (index of imprimitivity)
r = [0.0]

for tau in numpy.arange(tau_min, tau_max, passo_tau):
    percentage = 100 * tau / (tau_max-passo_tau)
    print("%.4f" % percentage, "%")
    for k in numpy.arange(k_min, k_max, passo_k):
        for u in range(0, m):
            # print(u)
            s = d * cmath.exp(2 * math.pi * 1j * u / m)

            z_1 = a - k * d * delta + (b-k*d*eta) * 1j + (1 / tau) * lambertw(tau * k * s * (delta + 1j * eta) * cmath.exp(-tau * (a - k * d * delta + (b-k*d*eta) * 1j)))
            # z_2 = a - k * d - b * 1j + (1 / tau) * lambertw(tau * k * s * cmath.exp(-tau * (a - k * d - b * 1j)))
            # print(z_1)

            r[u] = z_1.real

        # print(r)
        # y = (tau, k, max(r))
        if max(r) < 0:
            f.write("%f %f %f\n" % (tau, k, max(r)))
f.close()
