import unittest
from dfn import Fluid


class TestFluid(unittest.TestCase):

    """Tests for dfn.Fluid."""

    def __init__(self, *args, **kwargs):
        super(TestFluid, self).__init__(*args, **kwargs)
        self.rho = 1.0
        self.mu = 2.0
        self.c_f = 3.0

        self.fluid_1 = Fluid(density=self.rho, viscosity=self.mu)
        self.fluid_2 = Fluid(density=self.rho, viscosity=self.mu,
                             heat_capacity=self.c_f)

    def test_density_attr(self):
        self.assertEqual(self.rho, self.fluid_1.rho)

    def test_viscosity_attr(self):
        self.assertEqual(self.mu, self.fluid_1.mu)

    def test_default_value(self):
        self.assertIsNone(self.fluid_1.c_f)

    def test_heat_capacity_attr(self):
        self.assertEqual(self.c_f, self.fluid_2.c_f)


if __name__ == '__main__':
    unittest.main()
