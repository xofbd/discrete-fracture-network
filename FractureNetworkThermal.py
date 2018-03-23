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

        # raise error if mass flow needs to be calculated
        try:
            self.mass_flow
        except AttributeError:
            print "Network has not had mass flow calculated, call \
            'calculate_flow' before calling this method."

        # flip nodes if direction is wrong
        for i, seg in enumerate(self.connectivity):
            if self.mass_flow[i] < 0:
                self.connectivity[i] = seg[::-1]

        self.corrected_network = True

        return self

    def construct_graph(self):
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

    def __find_paths(self, inlet):
        n_inj = self.injection_node
        path_nodes = [tuple(p)
                      for p in nx.all_simple_paths(self.graph, n_inj, inlet)]
        path_nodes = set(path_nodes)

        # get segments in each path described by nodes
        path_segments = []
        for nodes in path_nodes:
            segment_choices = []
            for i in xrange(len(nodes) - 1):
                seg_inlet, seg_outlet = nodes[i], nodes[i + 1]
                d = G.get_edge_data(seg_inlet, seg_outlet)
                segment = [val['index'] for val in d.values()]
                segment_choices.append(segment)

                path = list(itertools.product(*segment_choices))
                path_segments.extend(path)

        return path_segments

    def calculate_temperature(self, segment, length, time):
        z, t = np.meshgrid(length, time)

        inlet = self.connectivity[segment, 0]
        outlet = self.connectivity[segment, 1]

        chi = __mass_contribution()
        k_r = self.thermal_cond
        alpha_r = self.thermal_diff
        m = self.mass_flow
        H = self.thickness
        L = self.length
        cp_f = self.c_f

        beta = 2 * k_r * H / (m * cp_f)
        xi = beta * L / (2 * np.sqrt(alpha_r * t))

        for path in nx.all_simple_paths(self.graph, inlet):
            S_k = self.__path_segments(path)
            chi_prod = chi[S_k].prod()
            xi_eff = xi[S_k].sum() + beta[
                segment] * z / (2 * np.sqrt(alpha_r * t))
            Theta += chi_prod * erf(xi_eff)
