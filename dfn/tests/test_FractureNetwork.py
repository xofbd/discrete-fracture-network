import copy
import unittest

import numpy as np

from dfn import FractureNetwork


class TestFractureNetwork(unittest.TestCase):

    """Test dfn.FractureNetwork."""

    def __init__(self, *args, **kwargs):
        super(TestFractureNetwork, self).__init__(*args, **kwargs)

        conn = [(0, 1), (1, 2), (1, 2), (2, 3)]
        L = range(1, 5)
        H = range(1, 5)
        w = range(1, 5)
        self.network = FractureNetwork(conn, L, H, w)

    def test_set_number_of_nodes(self):
        """Test determining the number of nodes in the fracture network."""
        self.assertEqual(self.network.n_nodes, 4)

    def test_check_parameters(self):
        """Test checking valid inputs, same size and greater than zero."""

        for attr in ('length', 'thickness', 'width'):
            # create FractureNetwork objects with bad attributes (ValueError)
            net_wrong_n_nodes = copy.copy(self.network)
            net_neg_value = copy.copy(self.network)
            net_zero_value = copy.copy(self.network)

            setattr(net_wrong_n_nodes, attr, np.array(range(1, 4)))
            setattr(net_neg_value, attr, np.array([-1, 2, 3, 4]))
            setattr(net_zero_value, attr, np.array([0, 2, 3, 4]))

            for network in (net_wrong_n_nodes, net_neg_value, net_zero_value):
                with self.assertRaises(ValueError):
                    network._check_parameters()

if __name__ == '__main__':
    unittest.main()
