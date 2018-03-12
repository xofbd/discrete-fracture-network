class Fluid(object):

    """Fluid to flow in the fracture network.

    Attributes
    ----------
    rho : float
        Density of the fluid.

    mu : float
        Dynamic viscosity of the fluid.
    """

    def __init__(self, rho, mu):
        self.rho = rho
        self.mu = mu
