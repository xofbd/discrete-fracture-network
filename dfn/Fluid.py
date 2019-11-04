from dfn import Base


class Fluid(Base):

    """Fluid to flow in the fracture network.

    Attributes
    ----------
    rho : float
        Density of the fluid.

    mu : float
        Dynamic viscosity of the fluid.

    c_f : float
        Specific heat capacity of fluid. The heat capacity can either be
        constant pressure or constant value, whatever makes sense for the fluid
        and the process.
    """

    def __init__(self, density=None, viscosity=None, heat_capacity=None):
        self.rho = density
        self.mu = viscosity
        self.c_f = heat_capacity

    def __eq__(self, other):
        attributes = ['rho', 'mu', 'c_f']
        return self._check_attr_equality(other, attributes)
