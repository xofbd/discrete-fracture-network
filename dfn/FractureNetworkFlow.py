import numpy as np

from dfn import FractureNetwork


class FractureNetworkFlow(FractureNetwork):

    """Discrete fracture network model for flow.

    Parameters
    ----------
    connectivity : array-like
        The inlet and outlet node pairing for each segment, describing the
        network connectivity. The first and second elements of the pair are
        inlet and outlet node, respectively. Its length is the number of
        segments in the fracture network.

    length : array-like
        The length of each segment.

    thickness : array-like
        The thickness of each segment, the dimension orthogonal to the fracture
        network.

    width : array-like
        The fracture width or aperture.

    Attributes
    ----------
    n_segments : int
        The number of segments in the fracture network.

    n_nodes : int
        The number of nodes in the fracture network.

    fluid : dfn.Fluid
        The fluid flowing through the network.

    conductance : numpy.ndarray
        The conductance of each segment.

    pressure : numpy.ndarray
        The pressure at each node.

    mass_flow : numpy.ndarray
        The mass flow rate for each segment.

    corrected_network : Boolean
        Whether the designation of inlet and outlet nodes of the segments have
        been checked and corrected.
    """

    def __init__(self, connectivity, length, thickness, width):
        super(FractureNetworkFlow, self).__init__(connectivity, length,
                                                  thickness, width)
        self.fluid = None
        self.conductance = None
        self.pressure = None
        self.mass_flow = None
        self.corrected_network = False

    def calculate_flow(self, fluid, essential_bc, point_sources, correct=False):
        """Calculate the mass flow throughout the fracture network.

        The mass flow is calculated by first finding the pressure at each node.
        The flow through each segment is equal to the pressure change times
        the conductance. The pressure change is defined as outlet minus inlet
        pressure value of the segment, where the outlet node is the second node
        in defining the segment's connectivity. Note, that negative values of
        mass flow are a result of a reverse of the specified inlet and outlet
        nodes. For example, if a segment's nodes are (0, 1) and has a negative
        mass flow, then the correct ordering is (1, 0).

        Parameters
        ----------
        fluid : dfn.Fluid
            Fluid object containing the fluid's properties.

        essential_bc : dict
            Dictionary of node index to pressure at node for all essential
            (Dirichlet) boundary conditions.

        point_sources : dict
            Dictionary of node index to mass rate loss or gain for all point
            sources in the system. Negative values indicate that mass is
            removed from the system; positive values indicate that mass is
            injected.

        correct : Boolean
            Whether to correct the inlet and outlet node designation for the
            segments.

        Returns
        -------
        self : object
            Returns self.
        """

        # check inputs
        ebc_node_set = set(essential_bc.keys())
        ps_node_set = set(point_sources.keys())
        problem_nodes = ebc_node_set.intersection(ps_node_set)

        if bool(problem_nodes):
            error_str = ("A node cannot have an essential boundary condition "
                         "imposed and be a point source/sink. The problem "
                         "nodes are {}.").format(list(problem_nodes))
            raise ValueError(error_str)

        # solve for pressure
        self.fluid = fluid
        self._solve_pressure(essential_bc, point_sources)
        outlets = self.connectivity[:, 1]
        inlets = self.connectivity[:, 0]

        # calculate mass flow
        Delta_P = self.pressure[outlets] - self.pressure[inlets]
        self.mass_flow = -self.conductance * Delta_P

        if correct:
            return self.correct_direction()

        return self

    def _solve_pressure(self, essential_bc, point_sources):
        """Solve for the pressure at each node of the fracture network.

        The pressure is solved by applying mass conservation around each node
        of the fracture network. The result is a system of linear equations in
        the form of [D]{P} = {f}, where {P} is the pressure at each node.
        """

        D = self._assemble_D()
        f = self._assemble_f(point_sources)

        # apply essential boundary conditions (Dirichlet)
        nodes = list(essential_bc.keys())
        values = list(essential_bc.values())

        f -= np.dot(D[:, nodes], values)
        f[nodes] = values

        D[nodes, :] = 0
        D[:, nodes] = 0
        D[nodes, nodes] = 1

        self.pressure = np.linalg.solve(D, f)

        return self

    def _assemble_D(self):
        """Assemble the conductance (coefficient) matrix."""

        D = np.zeros((self.n_nodes, self.n_nodes))
        elemental_D = np.array([[1, -1], [-1, 1]])
        self._calculate_conductance()

        # assemble the global conductance matrix
        for i, seg in enumerate(self.connectivity):
            D_e = self.conductance[i] * elemental_D
            D[np.ix_(seg, seg)] += D_e

        return D

    def _assemble_f(self, point_sources):
        """Assemble the source vector."""

        f = np.zeros(self.n_nodes)

        nodes = list(point_sources.keys())
        values = list(point_sources.values())
        f[nodes] = values

        return f

    def _calculate_conductance(self):
        """Calculate the conduction for each segment of the network."""

        num = self.fluid.rho * self.width**3 * self.thickness
        denom = 12 * self.fluid.mu * self.length
        self.conductance = num / denom

        return self

    def correct_direction(self):
        """Correct the order of the inlet and outlet nodes (direction).

        The first entry in a segment's connectivity is the inlet node and the
        second is the segment's outlet. However, the connectivity array is
        usually defined before one knows the flow structure in the network. If
        the calculated flow in the segment is negative, then the designation of
        the inlet and outlet nodes is reversed and flow value is made positive.

        Returns
        -------
        self : object
            Returns self.
        """

        # raise error if mass flow needs to be calculated
        if self.mass_flow is None:
            error_str = ("Network has not had the mass flow calculated, call "
                         "'calculate_flow' before calling this method.")
            raise TypeError(error_str)

        # flip nodes if the direction is wrong (negative flow)
        for i, seg in enumerate(self.connectivity):
            if self.mass_flow[i] < 0:
                self.connectivity[i] = seg[::-1]
                self.mass_flow[i] = -self.mass_flow[i]

        self.corrected_network = True

        return self
