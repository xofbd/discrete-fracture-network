import numpy as np

from base import Base


class FractureNetwork(Base):

    """Discrete fracture network model.

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
    """

    def __init__(self, connectivity, length, thickness, width):
        self.connectivity = np.array(connectivity)
        self.length = np.array(length)
        self.thickness = np.array(thickness)
        self.width = np.array(width)
        self.n_segments = len(connectivity)

        self._set_number_of_nodes()
        self._check_parameters()

    def _set_number_of_nodes(self):
        """Get number of unique nodes in the fracture network."""

        node_set = set()
        for a, b in self.connectivity:
            node_set = node_set.union({a, b})
        self.n_nodes = len(node_set)

        return self

    def _check_parameters(self):
        """Check that the object's parameters are valid."""

        params = ('length', 'thickness', 'width')
        error_str_1 = "The parameter(s) {} needs to be greater than zero."
        error_str_2 = ("The size of parameter(s) {} needs to be equal to the "
                       "number of segments as indicated by the connectivity.")

        # parameters must be positive
        bad_params = [p for p in params if (getattr(self, p) <= 0).any()]

        if len(bad_params) > 0:
            raise ValueError(error_str_1.format(", ".join(bad_params)))

        # parameter size needs to be equal to the number of segments
        bad_params = [p for p in params
                      if getattr(self, p).size != self.n_segments]

        if len(bad_params) > 0:
            raise ValueError(error_str_2.format(", ".join(bad_params)))

        return self

    def __eq__(self, other):
        attributes = ['connectivity', 'length', 'thickness', 'width']
        return _check_attr_equality(self, other, attributes)
