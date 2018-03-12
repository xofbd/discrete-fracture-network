import numpy as np


class Fluid(object):

    def __init__(self, rho, mu):
        self.rho = rho
        self.mu = mu


class FractureNetworkFlow(object):

    def __init__(self, conn, length, height, width):
        self.conn = np.array(conn)
        self.length = np.array(L)
        self.height = np.array(H)
        self.width = np.array(w)

        self.n_nodes = 1 + self.conn.max()

    def __calculate_conductance(self):
        """Calculate the conduction for each segment of the network."""

        rho = self.fluid.rho
        mu = self.fluid.mu

        self.conductance = rho * self.w**3 * self.H / (12 * mu * self.L)

    def __assemble_D(self):
        """Assemble the conductance (coefficient) matrix.

        The conductance matrix results from applying mass conservation around
        each node.
        """

        self.__calculate_conductance()
        self._D = np.zeros((self.n_nodes, self.n_nodes))
        elemental_D = np.array([[1, -1], [-1, 1]])

        for i, seg in enumerate(self.conn):
            D_e = self.conductance[i] * elemental_D
            self._D[np.ix_(seg, seg)] += D_e

    def __assemble_f(self):
        """Assemble the source vector."""

        self._f = np.zeros(self.n_nodes)

        nodes = self.point_sources.keys()
        values = self.point_sources.values()
        self._f[nodes] = values

    def __apply_EBC(self):
        """Apply essential (Dirichlet) boundary conditions.

        Applying essential boundary conditions alters both the conductance
        (coefficient) matrix and the source vector.
        """

        nodes = self.essential_bc.keys()
        values = self.essential_bc.values()

        self._f -= np.dot(self._D[:, nodes], values)
        self._f[nodes] = values

        self._D[nodes, :] = 0
        self._D[:, nodes] = 0
        self._D[nodes, nodes] = 1

    def solve_pressure(self):
        """Solve for the pressure at each node of the fracture network.

        The pressure is solved by applying mass conservation around each node
        of the fracture network. The result is a system of equations in the
        form of [D]{P} = {f}, where {P} is the pressure at each node.
        """

        self.__assemble_D()
        self.__assemble_f()
        self.__apply_EBC()

        self.pressure = np.linalg.solve(self._D, self._f)

    def calculate_flow(self, fluid, essential_bc, point_sources):
        """
        Calculate the mass flow throughout the fracture network.

        Parameters
        ----------
        fluid : dfn.Fluid
            Fluid object containing the fluid properties

        essential_bc : dict
            dictionary of node index to pressure at node

        point_sources : dict
            dictionary of node index to mass rate loss or gain

        Returns
        -------
        mass flow rate : numpy.darray
            mass flow rate for each segment of the fracture network.
        """

        self.fluid = fluid
        self.essential_bc = essential_bc
        self.point_sources = point_sources

        self.solve_pressure()
        Delta_P = self.pressure[self.conn[:, 1]] - self.P[self.conn[:, 0]]

        return -self.conductance * Delta_P

if __name__ == '__main__':
    conn = [(0, 1), (1, 2), (1, 2), (2, 3)]
    L = [1, 1, 0.5, 1]
    H = [1, 1, 1, 1]
    w = [1, 1, 1, 1]

    fluid = Fluid(1000.0, 1E-3)
    network = FractureNetworkFlow(conn, L, H, w)
    essential_bc = {0: 1}
    point_sources = {3: -10.0}
    m = network.calculate_flow(fluid, essential_bc, point_sources)
