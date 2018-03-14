class Fluid(object):

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

    def __init__(self, rho, mu, c_f):
        self.rho = rho
        self.mu = mu
        self.c_f = c_f
