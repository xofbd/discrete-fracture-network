import numpy as np


class FractureNetworkFlow(object):

    """Discrete fracture network model for flow.

    Parameters
    ----------
    conn : iterable
        Connectivity of the network, each item is an iterable of describing the
        inlet and outlet nodes of the segment. Its length is the number of
        segments in the fracture network.

    length : iterable
        The length of each segment..

    height : iterable
        The thickness of each segment, the dimension orthogonal to the fracture
        network.

    width : iterable
        The fracture width or aperture.

    Attributes
    ----------
    n_segments : int
        The number of segments in the fracture network.

    n_nodes : int
        The number of nodes in the fracture network.

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
        """Assemble the conductance (coefficient) matrix."""

        D = np.zeros((self.n_nodes, self.n_nodes))
        elemental_D = np.array([[1, -1], [-1, 1]])
        self.__calculate_conductance()

        # assemble the global conductance matrix from elemental/segment parts
        for i, seg in enumerate(self.conn):
            D_e = self.conductance[i] * elemental_D
            D[np.ix_(seg, seg)] += D_e

        return D

    def __assemble_f(self):
        """Assemble the source vector."""

        f = np.zeros(self.n_nodes)

        nodes = self.point_sources.keys()
        values = self.point_sources.values()
        f[nodes] = values

        return f

    def __assemble_SLAE(self):
        """Assemble the system of linear algebraic equations (SLAE)."""

        D = __assemble_D()
        f = __assemble_f()

        # applying essential boundary conditions (Dirichlet)
        nodes = self.essential_bc.keys()
        values = self.essential_bc.values()

        f -= np.dot(self._D[:, nodes], values)
        f[nodes] = values

        D[nodes, :] = 0
        D[:, nodes] = 0
        D[nodes, nodes] = 1

        return D, f

    def solve_pressure(self):
        """Solve for the pressure at each node of the fracture network.

        The pressure is solved by applying mass conservation around each node
        of the fracture network. The result is a system of linear equations in
        the form of [D]{P} = {f}, where {P} is the pressure at each node.
        """

        D, f = self.__assemble_SLAE()
        self.pressure = np.linalg.solve(D, f)

    def calculate_flow(self, fluid, essential_bc, point_sources):
        """Calculate the mass flow throughout the fracture network.

        The mass flow is calculated by first finding the pressure at each node.
        The flow through each segment is equal to the pressure change times
        the conductance. The pressure change is defined as outlet minus inlet
        pressure value of the segment, where the outlet node is the second node
        in defining the segment's connectivity. Note, that negative values of
        mass flow are a result of a reverse of the specified inlet and outlet
        nodes. For example, if a segment nodes are (0, 1) and has a negative
        mass flow, then the correct ordering is (1, 0).

        Parameters
        ----------
        fluid : dfn.Fluid
            Fluid object containing the fluid's properties

        essential_bc : dict
            Dictionary of node index to pressure at node for all essential
            (Dirichlet) boundary conditions.

        point_sources : dict
            Dictionary of node index to mass rate loss or gain for all point
            sources in the system. Negative values indicate that mass is
            removed from the system; positive values indicate that mass is
            injected.

        Returns
        -------
        self : object
            Returns self
        """

        self.fluid = fluid
        self.essential_bc = essential_bc
        self.point_sources = point_sources

        self.solve_pressure()
        Delta_P = self.pressure[self.conn[:, 1]] - self.P[self.conn[:, 0]]
        self.mass_flow = -self.conductance * Delta_P

        return self
