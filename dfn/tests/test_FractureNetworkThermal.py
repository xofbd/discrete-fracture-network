import copy
import unittest

import networkx as nx
import numpy as np
from scipy.special import erf

from dfn import Fluid, FractureNetworkThermal


class TestFractureNetworkThermal(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestFractureNetworkThermal, self).__init__(*args, **kwargs)

        # fluid properties
        cp_w = 4300.0
        rho_w = 1000.0
        mu_w = 1E-3
        self.fluid = Fluid(density=rho_w, viscosity=mu_w, heat_capacity=cp_w)

        # reservoir properties
        k_r = 2.9
        cp_r = 1050.0
        rho_r = 2700.0
        alpha_r = k_r / (rho_r * cp_r)

        # first network
        conn_1 = [(0, 1), (1, 2), (1, 3), (2, 4), (3, 4), (4, 5)]
        L_1 = [100, 500, 500, 500, 500, 100]
        H_1 = [500, 500, 500, 500, 500, 500]
        w_1 = [1E-3, 1E-3, 1E-3, 1E-3, 1E-3, 1E-3]
        self.network_1 = FractureNetworkThermal(conn_1, L_1, H_1, w_1, k_r,
                                                alpha_r)

        # second network
        conn_2 = [(0, 1), (1, 2), (2, 3), (1, 4), (2, 5), (3, 6), (4, 5),
                  (5, 6), (4, 7), (5, 8), (6, 9), (7, 8), (8, 9), (9, 10)]
        L_2 = 250 * np.ones(len(conn_2))
        L_2[0] = 100
        L_2[-1] = 100
        H_2 = 500 * np.ones(len(conn_2))
        w_2 = 1E-3 * np.ones(len(conn_2))
        self.network_2 = FractureNetworkThermal(conn_2, L_2, H_2, w_2, k_r,
                                                alpha_r)

    def copy_networks(self):
        """Return a copy of the fracture networks."""
        return copy.copy(self.network_1), copy.copy(self.network_2)

    def networks_with_flow(self):
        """Return networks with the mass flow calculated."""

        network_1, network_2 = self.copy_networks()

        P_0 = 0.0
        m_inj = 50.0

        network_1.calculate_flow(self.fluid, {0: P_0}, {5: -m_inj})
        network_2.calculate_flow(self.fluid, {0: P_0}, {10: -m_inj})

        return network_1, network_2

    def reverse_nodes(self, network, segments):
        """Reverse the node order for given segments."""

        conn = network.connectivity

        for seg in segments:
            inlet, outlet = conn[seg]
            conn[seg, :] = outlet, inlet

        network.connectivity = conn
        return network

    def test_no_mass_flow(self):
        """Test if TypeError is raised for networks without flow calculated."""

        with self.assertRaises(TypeError):
            self.network_1._check_if_calculated()

        with self.assertRaises(TypeError):
            self.network_2._check_if_calculated()

    def test_neg_mass_flow(self):
        """Test if valueError is raised for networks with negative flow."""

        network_1, network_2 = self.networks_with_flow()

        network_1 = self.reverse_nodes(network_1, [1])
        network_2 = self.reverse_nodes(network_2, [1])

        network_1.calculate_flow(self.fluid, {0: 0}, {5: -1.0})
        network_2.calculate_flow(self.fluid, {0: 0}, {10: -1.0})

        with self.assertRaises(ValueError):
            network_1.calculate_temperature(self.fluid, 0, [0], [1])

        with self.assertRaises(ValueError):
            network_2.calculate_temperature(self.fluid, 0, [0], [1])

    def test_construct_graph(self):
        """Test _construct_graph method."""

        network_1, network_2 = self.networks_with_flow()

        network_1._construct_graph()
        network_2._construct_graph()

        # construct graph for network 1
        G_1 = nx.MultiDiGraph()
        edge_data_1 = [(0, 1, {'index': 0}), (1, 2, {'index': 1}),
                       (1, 3, {'index': 2}), (2, 4, {'index': 3}),
                       (3, 4, {'index': 4}), (4, 5, {'index': 5})]
        G_1.add_edges_from(edge_data_1)

        # construct graph for network 2
        G_2 = nx.MultiDiGraph()
        edge_data_2 = [(0, 1, {'index': 0}), (1, 2, {'index': 1}),
                       (2, 3, {'index': 2}), (1, 4, {'index': 3}),
                       (2, 5, {'index': 4}), (3, 6, {'index': 5}),
                       (4, 5, {'index': 6}), (5, 6, {'index': 7}),
                       (4, 7, {'index': 8}), (5, 8, {'index': 9}),
                       (6, 9, {'index': 10}), (7, 8, {'index': 11}),
                       (8, 9, {'index': 12}), (9, 10, {'index': 13})]
        G_2.add_edges_from(edge_data_2)

        # return True if graphs are the same
        is_isomorphic_1 = nx.is_isomorphic(network_1.graph, G_1)
        is_isomorphic_2 = nx.is_isomorphic(network_2.graph, G_2)

        self.assertTrue(is_isomorphic_1)
        self.assertTrue(is_isomorphic_2)

    def test_find_injection_nodes(self):
        """Test _find_injection_nodes method."""

        network_1, network_2 = self.networks_with_flow()
        network_1._construct_graph()
        network_2._construct_graph()

        self.assertEqual(network_1._find_injection_nodes(), [0])
        self.assertEqual(network_2._find_injection_nodes(), [0])

    def test_mass_contribution(self):
        """Test _mass_contribution method."""

        network_1, network_2 = self.networks_with_flow()

        chi_1 = network_1._mass_contribution()
        chi_2 = network_2._mass_contribution()

        # first network
        for i in (0, 1, 2, 5):
            self.assertAlmostEqual(chi_1[i], 1.0, 12)

        self.assertAlmostEqual(chi_1[3] + chi_1[4], 1.0, 12)

        # second network
        for i in (0, 1, 2, 3, 8, 13):
            self.assertAlmostEqual(chi_2[i], 1.0, 12)

        for i, j in [(4, 6), (5, 7), (9, 11), (10, 12)]:
            self.assertAlmostEqual(chi_2[i] + chi_2[j], 1.0, 12)

    def test_find_paths(self):
        """Test find_paths method."""

        network_1, network_2 = self.networks_with_flow()

        path_1 = {(0, 1, 3), (0, 2, 4)}
        path_2 = {(0, 1, 2, 5, 10), (0, 1, 4, 7, 10), (0, 3, 6, 7, 10),
                  (0, 3, 6, 9, 12), (0, 3, 8, 11, 12), (0, 1, 4, 9, 12)}

        self.assertEqual(path_1, set(network_1.find_paths(0, 4)))
        self.assertEqual(path_2, set(network_2.find_paths(0, 9)))

    def test_calculate_temperature_inlet_segment(self):
        """Test calculate_temperature ability to handle the inlet segment."""

        # operational parameters for temperature
        t_end = 86400 * 365.25 * 20
        time = t_end * np.linspace(1.0 / 100, 1.0, 100)

        distance = np.linspace(0.0, 100.0, 100)
        z, t = np.meshgrid(distance, time)

        network_1, network_2 = self.networks_with_flow()

        # create parameters for temperature manually
        m_1 = network_1.mass_flow[0]
        m_2 = network_2.mass_flow[0]
        beta_1 = 2 * network_1.thermal_cond * network_1.thickness[0] / \
            (m_1 * network_1.fluid.c_f)
        beta_2 = 2 * network_2.thermal_cond * network_2.thickness[0] / \
            (m_2 * network_2.fluid.c_f)

        xi_1 = beta_1 * z / (2 * np.sqrt(network_1.thermal_diff * t))
        xi_2 = beta_2 * z / (2 * np.sqrt(network_2.thermal_diff * t))

        Theta_1 = erf(xi_1)
        Theta_2 = erf(xi_2)

        # difference between manual and automatic construction
        diff_1 = Theta_1 - network_1.calculate_temperature(self.fluid, 0,
                                                           distance, time)
        diff_2 = Theta_2 - network_2.calculate_temperature(self.fluid, 0,
                                                           distance, time)

        self.assertAlmostEqual((diff_1**2).sum() / (Theta_1**2).sum(), 0, 12)
        self.assertAlmostEqual((diff_2**2).sum() / (Theta_2**2).sum(), 0, 12)

    def test_calculate_temperature(self):
        """Test calculate_temperature by constructing manual the equations."""

        # operational parameters for temperature
        t_end = 86400 * 365.25 * 20
        time = t_end * np.linspace(1.0 / 100, 1.0, 100)

        distance = np.linspace(0.0, 100.0, 100)
        z, t = np.meshgrid(distance, time)

        network_1, network_2 = self.networks_with_flow()

        # create parameters for temperature manually
        chi_1 = np.array([1.0, 1.0, 1.0, 0.5, 0.5, 1.0])
        chi_2 = np.ones(network_2.n_segments)
        chi_2[4:8] = 0.5
        chi_2[9:13] = 0.5

        m_1 = network_1.mass_flow
        m_2 = network_2.mass_flow
        beta_1 = 2 * network_1.thermal_cond * network_1.thickness / \
            (m_1 * network_1.fluid.c_f)
        beta_2 = 2 * network_2.thermal_cond * network_2.thickness / \
            (m_2 * network_2.fluid.c_f)

        xi_1 = np.einsum('i,jk->ijk', beta_1 * network_1.length,
                         1 / (2 * np.sqrt(network_1.thermal_diff * t)))
        xi_2 = np.einsum('i,jk->ijk', beta_2 * network_2.length,
                         1 / (2 * np.sqrt(network_2.thermal_diff * t)))
        a = xi_1[[0, 2, 4], :, :].sum(axis=0)
        b = xi_1[[0, 1, 3], :, :].sum(axis=0)
        xi_seg = beta_1[-1] * z / (2 * np.sqrt(network_1.thermal_diff * t))

        Theta_1 = chi_1[0] * chi_1[2] * chi_1[4] * erf(a + xi_seg) + \
            chi_1[0] * chi_1[1] * chi_1[3] * erf(b + xi_seg)

        a = xi_2[[0, 1, 2, 5, 10], :, :].sum(axis=0)
        b = xi_2[[0, 1, 4, 7, 10], :, :].sum(axis=0)
        c = xi_2[[0, 3, 6, 7, 10], :, :].sum(axis=0)
        d = xi_2[[0, 3, 6, 9, 12], :, :].sum(axis=0)
        e = xi_2[[0, 3, 8, 11, 12], :, :].sum(axis=0)
        f = xi_2[[0, 1, 4, 9, 12], :, :].sum(axis=0)

        C_1 = chi_2[0] * chi_2[1] * chi_2[2] * chi_2[5] * chi_2[10]
        C_2 = chi_2[0] * chi_2[1] * chi_2[4] * chi_2[7] * chi_2[10]
        C_3 = chi_2[0] * chi_2[3] * chi_2[6] * chi_2[7] * chi_2[10]
        C_4 = chi_2[0] * chi_2[3] * chi_2[6] * chi_2[9] * chi_2[12]
        C_5 = chi_2[0] * chi_2[3] * chi_2[8] * chi_2[11] * chi_2[12]
        C_6 = chi_2[0] * chi_2[1] * chi_2[4] * chi_2[9] * chi_2[12]
        xi_seg = beta_2[-1] * z / (2 * np.sqrt(network_2.thermal_diff * t))

        Theta_2 =  C_1 * erf(a + xi_seg) + C_2 * erf(b + xi_seg) + \
            C_3 * erf(c + xi_seg) + C_4 * erf(d + xi_seg) + \
            C_5 * erf(e + xi_seg) + C_6 * erf(f + xi_seg)

        # difference between manual and automatic construction
        diff_1 = Theta_1 - network_1.calculate_temperature(self.fluid, 5,
                                                           distance, time)
        diff_2 = Theta_2 - network_2.calculate_temperature(self.fluid, 13,
                                                           distance, time)

        self.assertAlmostEqual((diff_1**2).sum() / (Theta_1**2).sum(), 0, 12)
        self.assertAlmostEqual((diff_2**2).sum() / (Theta_2**2).sum(), 0, 12)

if __name__ == '__main__':
    unittest.main()
