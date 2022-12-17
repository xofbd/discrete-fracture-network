import unittest
from dfn import Fluid


class TestBase(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestBase, self).__init__(*args, **kwargs)
        self.fluid_1 = Fluid(1, 2, 3)
        self.fluid_2 = Fluid(1, 2, 3)
        self.fluid_3 = Fluid(1, 2, 2)

    def test_equality(self):
        """Test equality of objects when all attrs are equal."""
        self.assertEqual(self.fluid_1, self.fluid_2)

    def test_inequality(self):
        """Test inequality of objects when attrs differ."""
        self.assertNotEqual(self.fluid_1, self.fluid_3)


if __name__ == '__main__':
    unittest.main()
