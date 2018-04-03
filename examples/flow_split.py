"""
author: Don Fox
file name: flow_fraction.py

This script calculates the flow split fraction in a fracture network with two
paths. The network has four segments, were the flow first splits into two
segments and then combines again. The script generates an ensemble of networks,
each with a different set of aperture values, sampled from a log-normal
distribution.
"""

from math import log, sqrt

import numpy as np

from dfn import Fluid, FractureNetworkFlow


def calculate_flow_fraction(network, essential_bc, point_sources):
    # simulation parameters
    n_sims = 10000
    w_0 = 1E-3
    sigma_over_w_0 = 0.25
    flow_fraction = []

    # mean and variance of underlying normal distribution
    mu = log(w_0 / sqrt(1 + sigma_over_w_0**2))
    var = log(1 + sigma_over_w_0**2)
    w = np.random.lognormal(mu, sqrt(var), (4, n_sims))

    for i in xrange(n_sims):
        network.width = w[:, i]
        network.calculate_flow(fluid, essential_bc, point_sources)
        m = network.mass_flow
        flow_fraction.append(m[1] / (m[1] + m[2]))

    return flow_fraction

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import time

    start_time = time.time()

    # fluid properties
    rho_w = 1000.0
    mu_w = 1E-3

    # operational properties
    m_inj = 50.0
    P_inj = 0.0

    # network properties
    n_segs = 4
    conn = [(0, 1), (1, 2), (1, 2), (2, 3)]
    L = [100, 500, 500, 100]
    H = [500, 500, 500, 500]
    w = [1E-3, 1E-3, 1E-3, 1E-3]
    n_inj = 0
    n_prod = 3

    # create network object
    fluid = Fluid(density=rho_w, viscosity=mu_w, heat_capacity=None)
    network = FractureNetworkFlow(conn, L, H, w)

    # calculate flow in the network
    essential_bc = {n_inj: P_inj}
    point_sources = {n_prod: -m_inj}
    network.calculate_flow(fluid, essential_bc, point_sources, correct=True)

    # generate and plot results
    flow_fraction = calculate_flow_fraction(network, essential_bc,
                                            point_sources)

    f = plt.figure()
    plt.hist(flow_fraction, bins=50, normed=True)
    plt.xlabel('flow fraction (-)')
    plt.ylabel('probability density (-)')
    f.show()

    print "Elapsed time: {:.5} seconds".format(time.time() - start_time)
