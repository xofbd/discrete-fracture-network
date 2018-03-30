from itertools import product

import networkx as nx
import numpy as np
from scipy.special import erf

from FractureNetworkFlow import FractureNetworkFlow


class FractureNetworkThermal(FractureNetworkFlow):

    """Discrete fracture network model for flow and thermal performance.

    Parameters
    ----------
    connectivity : array-like
        Connectivity of the network, each item is an array-like of describing the
        inlet and outlet nodes of the segment. Its length is the number of
        segments in the fracture network.

    length : array-like
        The length of each segment..

    thickness : array-like
        The thickness of each segment, the dimension orthogonal to the fracture
        network.

    width : array-like
        The fracture width or aperture.

    thermal_cond : float
        Thermal conductivity of the reservoir

    thermal_diff : float
        Thermal diffusivity of the reservoir

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

    graph : networkx.MultiDiGraph
        Directed graph with possible multiple edges representation of the network.
    """

    def __init__(self, connectivity, length, thickness, width, thermal_cond,
                 thermal_diff):
        super(FractureNetworkThermal, self).__init__(
            connectivity, length, thickness, width)
        self.thermal_cond = thermal_cond
        self.thermal_diff = thermal_diff
        self.graph = None

    def calculate_temperature(self, fluid, segment, distance, time):
        """Calculate the (dimensionless) temperature for a segment.

        Parameters
        ----------
        fluid : dfn.Fluid
            Fluid injected into network.

        segment : int
            Index of segment of interest.

        length : array-like
            Distances to calculate temperature for the segment, measured from
            the segment's inlet node.

        time : array-like
            Time values to calculate the temperature for the segment, starting
            from the start of operation.

        Returns
        -------
        Theta : numpy.array, shape = (len(distance), len(time))
            Dimensionless temperature for each distance and time pairing.
        """

        # raise error from mass flow type or value errors
        if self.mass_flow is None:
            raise TypeError("Network has not had the mass flow calculated, "
                            "call 'calculate_flow' before calling this method.")

        if (self.mass_flow < 0).sum() > 0:
            raise ValueError("Network has negative mass flow values, need to correct"
                             " node designation using 'correct_direction' method or"
                             " set 'correct' to True in 'calculate_flow' method.")

        if self.graph is None:
            self._construct_graph()

        self.fluid = fluid
        inj_nodes = self._find_injection_nodes()

        z, t = np.meshgrid(distance, time)
        inlet = self.connectivity[segment, 0]

        # create local variables for ease
        chi = self._mass_contribution()
        k_r = self.thermal_cond
        alpha_r = self.thermal_diff
        m = self.mass_flow
        H = self.thickness
        L = self.length
        cp_f = self.fluid.c_f

        # useful dimensionless parameters
        beta = 2 * k_r * H / (m * cp_f)
        xi = np.einsum('i,jk -> ijk', beta * L, 1 / (2 * np.sqrt(alpha_r * t)))

        # loop through each path to the segment
        paths = self.find_paths(inj_nodes, inlet)
        Theta = 0

        for S_k in paths:
            S_k = list(S_k)

            chi_prod = chi[S_k].prod()
            xi_eff = xi[S_k, :].sum(axis=0) + beta[
                segment] * z / (2 * np.sqrt(alpha_r * t))
            Theta += chi_prod * erf(xi_eff)

        return Theta

    def _construct_graph(self):
        """Construct a NetworkX graph object that represents the network.

        A graph representation of the fracture network facilitates operations
        such as finding the flow paths from the injection and a specified
        segment. The fracture network with flow is a directed graph with
        possible multiedges, more than one edge/segment from one node to
        another.

        Returns
        -------
        self : object
            Returns self.
        """

        self.graph = nx.MultiDiGraph()
        edge_data = [(seg[0], seg[1], {'index': i})
                     for i, seg in enumerate(self.connectivity)]

        self.graph.add_edges_from(edge_data)

        return self

    def _find_injection_nodes(self):
        """Find injection node of the graph."""

        preds = self.graph
        inj_nodes = [node for node, p in preds.items() if bool(p) is False]

        return inj_nodes

    def _mass_contribution(self):
        """Calculate the segments' relative mass flow contribution.

        The mass flow contribution is the relative mass flow that a segment
        provides to its immediate upstream segment.
        """

        chi = np.zeros(self.n_segments)

        for i, seg in enumerate(self.connectivity):
            indices = self.connectivity[:, 1] == seg[1]
            chi[i] = self.mass_flow[i] / self.mass_flow[indices].sum()

        return chi

    def find_paths(self, n_inj, inlet):
        """Find all the paths from injection node to a segment's inlet node.

        This method relies on NetworkX's all_simple_paths function to find the
        paths. However, the function returns the nodal description for each of
        the paths. What is needed is segments that are part of each path.

        Parameters
        ----------
        n_inj : int
            The injection node.

        inlet : int
            The inlet node of the segment of interest.

        Returns
        -------
        path_segments : list
            A list of all paths to the segment's inlet node, where each entry
            is a tuple containing the segments that constitute the path.
        """

        # nx.all_simple_paths returns a generator and replicate paths as a
        # result of multiedges, need unique paths and as a list of tuples.
        path_nodes = set()
        for n_inj in inj_nodes:
            temp = [p for p in nx.all_simple_paths(self.graph, n_inj, inlet)]
            temp = map(tuple, temp)
            path_nodes = path_nodes.union(set(temp))

        # get segments in each path described by the nodes
        path_segments = []

        for nodes in path_nodes:
            # get all segments described by each inlet-outlet node pairing
            segment_choices = []

            for i in xrange(len(nodes) - 1):
                seg_inlet, seg_outlet = nodes[i], nodes[i + 1]

                # dictionary of segments with the inlet-outlet node paring
                d = self.graph.get_edge_data(seg_inlet, seg_outlet)

                segments = [val['index'] for val in d.values()]
                segment_choices.append(segments)

            # transform [[0], [1, 2]] to [(0, 1), (0, 2)]
            paths = list(product(*segment_choices))
            path_segments.extend(paths)

        return path_segments
