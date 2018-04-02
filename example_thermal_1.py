import matplotlib.pyplot as plt
import numpy as np

from Fluid import Fluid
from FractureNetworkThermal import FractureNetworkThermal

# fluid properties
cp_w = 4300
rho_w = 1000.0
mu_w = 1E-3

# reservoir properties
k_r = 2.9
cp_r = 1050.0
rho_r = 2700.0
alpha_r = k_r / (rho_r * cp_r)

# operational properties
m_inj = 50.0
P_inj = 0
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

fluid = Fluid(density=rho_w, viscosity=mu_w, heat_capacity=cp_w)
network = FractureNetworkThermal(conn, L, H, w, k_r, alpha_r)
essential_bc = {n_inj: P_inj}
point_sources = {n_prod: -m_inj}

network.calculate_flow(fluid, essential_bc, point_sources, correct=True)

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
