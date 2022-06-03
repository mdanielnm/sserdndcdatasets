#Boundary of the control region (stability islands)

import math
import numpy as np
from matplotlib import pyplot as plt

# Parameters of the linearized isolated vector field
a = 0.0153
b = 0.2715

# Parameters of the linearized coupling function
eta = 0
delta = 1

# Parameters of the regular network
d = 1 # degree of each vertex
m = 2 # imprimitivity index

# maximal number of stability islands
island_max = math.ceil(b*m * (1-math.cos(math.pi/(2*m))) / (2*a*math.pi) - 1/2)

# minimum and maximum values for kappa (coupling parameter):
if eta == 0:
    k_min = a / (2*d*delta)
else:
    k_min = (-a*delta + a * math.sqrt(delta**2 + eta**2)) / (d * eta**2)
    
if eta != 0:
    k_max = b / (d*eta)
else: 
    k_max = 10 #manually change this number if k_max is larger

# define the bifurcations curves:
def tau_1(j, k):
    psi = math.sqrt( (k*d*eta)**2 + 2*a*d*k*delta - a**2)
    z_1 = delta*(k*d*delta-a)-eta*psi-1j*(delta*psi+eta*(k*d*delta-a))
    return (2*math.pi*j/m - np.angle(z_1)) / (b-k*d*eta-psi)
        
def tau_2(j, k):
    psi = math.sqrt( (k*d*eta)**2 + 2*a*d*k*delta - a**2)
    z_2 = delta*(k*d*delta-a)+eta*psi-1j*(-delta*psi+eta*(k*d*delta-a))
    return (2*math.pi*(j+1)/m - np.angle(z_2)) / (b-k*d*eta+psi)



#auxiliary routine for finding intersection points using Bolzano-Cauchy theorem


for j in range(0, int(island_max)):
    
    k_min_found = k_min
    k_max_found = k_max
    k_1 = np.linspace(k_min_found+0.0001, k_max_found, 10000)
    dif_min = []
    dif_max = []

    for i in range(len(k_1)):
        dif_min.append(tau_1(j, k_1[i]) - tau_2(j, k_1[i]))
        if i>0 and dif_min[i-1] * dif_min[i]<0:
            k_min_found = (k_1[i-1]+k_1[i])/2
            break
        
    k_2 = np.linspace(2*k_min_found, k_max_found, 10000)
    for i in range(len(k_2)):
        dif_max.append(tau_1(j, k_2[i]) - tau_2(j, k_2[i]))
        if i>0 and dif_max[i-1] * dif_max[i]<0:
            k_max_found = (k_2[i-1]+k_2[i])/2
            break

    print("island index = ", j, ", k_min_found = ", k_min_found, ", k_max_found = ", k_max_found)

    # plot the stability boundaries
    k = np.linspace(k_min_found, k_max_found, 1000)
    curve_1 = []
    curve_2 = []
    for i in range(len(k)):
        curve_1.append(tau_1(j, k[i]))
        curve_2.append(tau_2(j, k[i]))

    plt.plot(curve_1, k, 'blue')
    plt.plot(curve_2, k, 'red')

plt.show()
