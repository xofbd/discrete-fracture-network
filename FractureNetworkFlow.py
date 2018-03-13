import numpy as np


class FractureNetworkFlow(object):

    """Discrete Fracture Network model for Flow.

    Parameters
    ----------
    conn : iterable
        Connectivity of the network, each item is an iterable of describing the
        inlet and outlet nodes of the segment. Its length is the number of
        segments in the fracture network.

    length : iterable
        The length of each segment..

    height: iterable
        The thickness of each segment, the dimension orthogonal to the fracture
        network.

    width: iterable
        The fracture width or aperture.

    n_nodes: int
        The number of nodes in the fracture network.

    Attributes
    ----------
    conductance : numpy.ndarray
        The conductance of each segment.
    pressure : numpy.ndarray
        The pressure at each node.
    """

    def __init__(self, conn, length, height, width):
        self.conn = np.array(conn)
        self.length = np.array(L)
        self.height = np.array(H)
        self.width = np.array(w)

        self.n_segments = len(conn)
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
        form of[D]{P} = {f}, where {P} is the pressure at each node.
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
        fluid:
            dfn.Fluid
            Fluid object containing the fluid properties

        essential_bc:
            dict
            dictionary of node index to pressure at node

        point_sources:
            dict
            dictionary of node index to mass rate loss or gain

        Returns
        -------
        mass flow rate:
            numpy.ndarray
            mass flow rate for each segment of the fracture network.
        """

        self.fluid = fluid
        self.essential_bc = essential_bc
        self.point_sources = point_sources

        self.solve_pressure()
        Delta_P = self.pressure[self.conn[:, 1]] - self.P[self.conn[:, 0]]

        return -self.conductance * Delta_P
