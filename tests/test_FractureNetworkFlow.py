import copy
import unittest

import numpy as np

from dfn import Fluid, FractureNetworkFlow


class TestFractureNetworkFlow(unittest.TestCase):

    """Test FractureNetworkFlow class."""

    def __init__(self, *args, **kwargs):
        super(TestFractureNetworkFlow, self).__init__(*args, **kwargs)

        conn = [(0, 1), (1, 2), (1, 2), (2, 3)]
        L = [1.0, 1.0, 1.0, 1.0]
        H = [1.0, 1.0, 1.0, 1.0]
        w = [1.0, 1.0, 1.0, 1.0]

        self.fluid = Fluid(density=1.0, viscosity=2.0, heat_capacity=3.0)
        self.network = FractureNetworkFlow(conn, L, H, w)
        self.network.fluid = self.fluid

        P_0 = 0.0
        m_inj = 1.0
        self.essential_bc = {0: P_0}
        self.point_sources = {3: -m_inj}

    def test_ebc_and_point_sources_inputs(self):
        """Test raising error for a point source and EBC on the same node."""

        m_inj = -list(self.point_sources.values())[0]
        point_sources = {0: -m_inj}

        with self.assertRaises(ValueError):
            self.network.calculate_flow(self.fluid, self.essential_bc,
                                        point_sources)

    def test_pressure_ebc(self):
        """Test applying pressure EBC on a node."""

        P_0 = list(self.essential_bc.values())[0]
        self.network._solve_pressure(self.essential_bc, self.point_sources)
        self.assertEqual(P_0, self.network.pressure[0])

    def test_solve_pressure(self):
        """Test behavior of solving for the pressure."""

        m_inj = -list(self.point_sources.values())[0]
        point_sources_double = {3: -2 * m_inj}

        # pressure should double when doubling the flow
        self.network._solve_pressure(self.essential_bc, self.point_sources)
        Delta_P = self.network.pressure[-1] - self.network.pressure[0]

        self.network._solve_pressure(self.essential_bc, point_sources_double)
        Delta_P_new = self.network.pressure[-1] - self.network.pressure[0]

        ratio = Delta_P_new / Delta_P

        self.assertAlmostEqual(ratio, 2, places=12)

    def test_mass_flow(self):
        """Test calculating mass flow by changing network parameters."""

        # base case
        self.network.calculate_flow(self.fluid, self.essential_bc,
                                    self.point_sources)
        m = self.network.mass_flow
        self.assertEqual(m[1], m[2])
        self.assertAlmostEqual(m[0], m[1] + m[2], places=12)

        # double length, double thickness
        network = copy.copy(self.network)
        network.length = np.array([1.0, 2.0, 1.0, 1.0])
        network.thickness = np.array([1.0, 2.0, 1.0, 1.0])
        network.calculate_flow(self.fluid, self.essential_bc,
                               self.point_sources)
        m = network.mass_flow
        self.assertAlmostEqual(m[1], m[2])
        self.assertAlmostEqual(m[0], m[1] + m[2], places=12)

        # double width, octuple length
        network = copy.copy(self.network)
        network.length = np.array([1.0, 8.0, 1.0, 1.0])
        network.width = np.array([1.0, 2.0, 1.0, 1.0])
        network.calculate_flow(self.fluid, self.essential_bc,
                               self.point_sources)
        m = network.mass_flow
        self.assertAlmostEqual(m[1], m[2], places=12)
        self.assertAlmostEqual(m[0], m[1] + m[2], places=12)

    def test_correct_flow(self):
        """Test that the flow and connectivity is corrected."""

        # check TypeError when mass flow has not been calculated
        with self.assertRaises(TypeError):
            self.network.correct_direction()

        # rearrange inlet-outlet designation to get negative flow values
        network = copy.copy(self.network)
        inlet, outlet = network.connectivity[1]
        network.connectivity[1] = (outlet, inlet)
        inlet, outlet = network.connectivity[-1]
        network.connectivity[-1] = (outlet, inlet)

        network.calculate_flow(self.fluid, self.essential_bc,
                               self.point_sources, correct=True)

        self.assertTrue((network.mass_flow >= 0).any())


if __name__ == '__main__':
    unittest.main()
