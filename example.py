if __name__ == '__main__':
    from dfn import FluidB
    from dfn import FractureNetworkFlow

    conn = [(0, 1), (1, 2), (1, 2), (2, 3)]
    L = [1, 1, 0.5, 1]
    H = [1, 1, 1, 1]
    w = [1, 1, 1, 1]

    fluid = Fluid(1000.0, 1E-3)
    network = FractureNetworkFlow(conn, L, H, w)
    essential_bc = {0: 1}
    point_sources = {3: -10.0}
    m = network.calculate_flow(fluid, essential_bc, point_sources)
