import numpy as np


class Base(object):

    def _check_attr_equality(self, other, attributes):
        """"Return True if all attributes have equal value."""

        for attr in attributes:
            if isinstance(getattr(self, attr), np.ndarray):
                if not np.all(getattr(self, attr) == getattr(other, attr)):
                    return False
            elif getattr(self, attr) != getattr(other, attr):
                return False

        return True
