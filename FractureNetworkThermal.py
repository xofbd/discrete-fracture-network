from FractureNetworkFlow import FractureNetworkFlow
import networkx as nx
import numpy as np
from scipy.special import erf
from itertools import product


class FractureNetworkThermal(FractureNetworkFlow):

    def __init__(self, connectivity, length, thickness, width, thermal_cond,
                 thermal_diff):
        self.connectivity = np.array(connectivity)
        self.length = np.array(length)
        self.thickness = np.array(thickness)
        self.width = np.array(width)
        self.thermal_cond = thermal_cond
        self.thermal_diff = thermal_diff

        self.n_segments = len(connectivity)
        self.n_nodes = 1 + self.connectivity.max()

    def correct_direction(self):
        """Correct the order of the inlet and outlet nodes (direction).

        The first entry in a segment's connectivity is the inlet node and the
        second is the segment's outlet. However, the connectivity array is
        usually defined before one knows the flow structure in the network. If
        the calculated flow in the segment is negative, then the designation of
        the inlet and outlet nodes is reversed.

        Returns
        -------
        self : object
            Returns self
        """

        # raise error if mass flow needs to be calculated
        try:
            self.mass_flow
        except AttributeError:
            print """Network has not had the mass flow calculated, call
                     'calculate_flow' before calling this method."""

        # flip nodes if the direction is wrong
        for i, seg in enumerate(self.connectivity):
            if self.mass_flow[i] < 0:
                self.connectivity[i] = seg[::-1]

        self.corrected_network = True

        return self

    def construct_graph(self):
        """Construct a NetworkX graph object that represents the network.

        A graph representation of the fracture network facilitates operations
        such as finding the flow paths from the injection and a specified
        segment. The fracture network with flow is a directed graph with
        possible multiedges, more than one edge/segment from one node to
        another.

        Returns
        -------
        self : object
            Returns self
        """

        self.graph = nx.MultiDiGraph()
        conn = [(seg[0], seg[1], {'index': i})
                for i, seg in enumerate(self.connectivity)]

        self.graph.add_edges_from(conn)

        return self

    def __mass_contribution(self):
        chi = np.zeros(self.n_segments)

        for i, seg in enumerate(self.connectivity):
            indices = self.connectivity[:, 1] == seg[1]
            chi[i] = self.mass_flow[i] / self.mass_flow[indices].sum()

        return chi

    def find_paths(self, n_inj, inlet):
        path_nodes = [tuple(p)
                      for p in nx.all_simple_paths(self.graph, n_inj, inlet)]
        path_nodes = set(path_nodes)

        # get segments in each path described by nodes
        path_segments = []
        for nodes in path_nodes:
            segment_choices = []
            for i in xrange(len(nodes) - 1):
                seg_inlet, seg_outlet = nodes[i], nodes[i + 1]
                d = self.graph.get_edge_data(seg_inlet, seg_outlet)
                segments = [val['index'] for val in d.values()]
                segment_choices.append(segments)

            paths = list(product(*segment_choices))
            path_segments.extend(paths)

        return path_segments

    def calculate_temperature(self, fluid, segment, length, time):
        z, t = np.meshgrid(length, time)

        inlet = self.connectivity[segment, 0]

        chi = self.__mass_contribution()
        k_r = self.thermal_cond
        alpha_r = self.thermal_diff
        m = self.mass_flow
        H = self.thickness
        L = self.length
        cp_f = fluid.c_f

        beta = 2 * k_r * H / (m * cp_f)
        xi = np.einsum('i,jk -> ijk', beta * L, 1 / (2 * np.sqrt(alpha_r * t)))
        paths = self.find_paths(0, inlet)
        Theta = 0

        for S_k in paths:
            # convert to list
            S_k = list(S_k)

            chi_prod = chi[S_k].prod()
            xi_eff = xi[S_k, :].sum(axis=0) + beta[
                segment] * z / (2 * np.sqrt(alpha_r * t))
            Theta += chi_prod * erf(xi_eff)

        return Theta
