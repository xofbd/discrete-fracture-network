from Fluid import Fluid
from FractureNetworkFlow import FractureNetworkFlow

# fluid properties
rho_w = 1000.0
mu_w = 1E-3

# operational properties
m_inj = 50.0
P_inj = 0

# network properties
n_segs = 4
conn = [(0, 1), (1, 2), (1, 2), (2, 3)]
L = [100, 250, 500, 100]
H = [500, 500, 500, 500]
w = [1E-3, 1E-3, 1E-3, 1E-3]
n_inj = 0
n_prod = 3

# create network object
fluid = Fluid(density=rho_w, viscosity=mu_w, heat_capacity=cp_w)
network = FractureNetworkThermal(conn, L, H, w, k_r, alpha_r)

# calculate flow in the network
essential_bc = {n_inj: P_inj}
point_sources = {n_prod: -m_inj}
network.calculate_flow(fluid, essential_bc, point_sources, correct=True)
print network.mass_flow
