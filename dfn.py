import numpy as np


class Fluid(object):

    def __init__(self, rho, mu):
        self.rho = rho
        self.mu = mu


class FractureNetwork(object):

    def __init__(self, conn, L, H, w):
        self.conn = np.array(conn)
        self.L = np.array(L)
        self.H = np.array(H)
        self.w = np.array(w)

        self.n_nodes = 1 + conn.max()

    def calculate_conductance(self):
        self.C = fluid.rho * self.w**3 * self.H / (12 * fluid.mu * self.L)

    def assemble_D(self):
        self.D = np.zeros((self.n_nodes, self.n_nodes))
        elemental_D = np.array([[1, -1], [-1, 1]])

        for i, seg in enumerate(self.conn):
            D_e = self.C[i] * elemental_D
            self.D[np.ix_(seg, seg)] += D_e

    def assemble_f(self):
        self.f = np.zeros(self.n_nodes)

        nodes = self.point_sources.keys()
        values = self.point_sources.values()
        self.f[nodes] = values

    def apply_EBC(self):
        nodes = self.essential_bc.keys()
        values = self.essential_bc.values()

        adj = -np.dot(self.D[:, nodes], values)
        self.f += adj
        self.f[nodes] = values

        self.D[nodes, :] = 0
        self.D[:, nodes] = 0
        self.D[nodes, nodes] = 1
        return self.D, self.f

    def solve_pressure(self):
        self.calculate_conductance()
        self.assemble_D()
        self.assemble_f()
        D, f = self.apply_EBC()

        self.P = np.linalg.solve(D, f)

    def calculate_flow(self, fluid, essential_bc, point_sources):
        """ Calculate flow throughout the fracture network.

        Parameters
        ----------
        fluid : dfn.Fluid
            Fluid object containing the fluid properties

        essential_bc : dict
            dictionary of node index to pressure at node

        point_sources : dict
            dictionary of node index to mass rate loss or gain
        """

        self.fluid = fluid
        self.essential_bc = essential_bc
        self.point_sources = point_sources

        self.solve_pressure()
        Delta_P = self.P[self.conn[:, 0]] - self.P[self.conn[:, 1]]
        return -self.C * Delta_P

if __name__ == '__main__':
    conn = [(0, 1), (1, 2), (1, 2), (2, 3)]
    L = [1, 1, 0.5, 1]
    H = [1, 1, 1, 1]
    w = [1, 1, 1, 1]

    fluid = Fluid(1000.0, 1E-3)
    network = FractureNetwork(conn, L, H, w)
    essential_bc = {0: 1}
    point_sources = {3: -10.0}
    P = network.calculate_flow(fluid, essential_bc, point_sources)
