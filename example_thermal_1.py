"""
author: Don Bruce Fox
file name: example_thermal_1.py

This example models the thermal performance of the fracture network with six
segments and two paths from the injection to production node. The fracture
network is the one depicted in Figure 4a of the following journal article:

    Fox, D. B., D. L. Koch, and J. W. Tester (2016), An analytical
    thermohydraulic model for discretely fractured geothermal reservoirs, Water
    Resources Research, 52, 6792-6817, doi: 10.1002/2016WR018666.

All values used in this example such as the thermophysical and operational
values are taken from the source above. The plotted results should match the
results in Figure 6 of the source.
"""

import matplotlib.pyplot as plt
import numpy as np

from Fluid import Fluid
from FractureNetworkThermal import FractureNetworkThermal

# fluid properties
cp_w = 4300.0
rho_w = 1000.0
mu_w = 1E-3

# reservoir properties
k_r = 2.9
cp_r = 1050.0
rho_r = 2700.0
alpha_r = k_r / (rho_r * cp_r)

# operational properties
m_inj = 50.0
P_inj = 0.0
t_end = 86400 * 365.25 * 20

# network properties
n_segs = 6
conn = [(0, 1), (1, 2), (1, 3), (2, 4), (3, 4), (4, 5)]
L = [100, 500, 500, 500, 500, 100]
H = [500, 500, 500, 500, 500, 500]
w = [1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3]
A = (2 * L[0] + 4 * L[1]) * H[0]
n_inj = 0
n_prod = 5

# create network object
fluid = Fluid(density=rho_w, viscosity=mu_w, heat_capacity=cp_w)
network = FractureNetworkThermal(conn, L, H, w, k_r, alpha_r)

# calculate flow in the network
essential_bc = {n_inj: P_inj}
point_sources = {n_prod: -m_inj}
network.calculate_flow(fluid, essential_bc, point_sources, correct=False)

# calculate temperature and plot results
segs_to_plot = (0, 1, 3, 5)
t = t_end * np.linspace(1.0 / 100, 1, num=100)
tau = k_r * rho_r * cp_r / (cp_w * m_inj / A)**2
Theta = np.zeros((len(segs_to_plot), len(t)))
f = plt.figure()

for i, seg in enumerate(segs_to_plot):
    z = np.array([L[seg]])
    Theta[i, :] = network.calculate_temperature(fluid, seg, z, t).ravel()
    plt.plot(t / tau, Theta[i, :], '--')

plt.ylim((0, 1))
plt.xlim((0, 3))
plt.ylabel('$\Theta$ (-)')
plt.xlabel('$\tau$ (-)')
f.show()
