import time
from math import log, sqrt

import matplotlib.pyplot as plt
import numpy as np

from Fluid import Fluid
from FractureNetworkThermal import FractureNetworkThermal


def run_simulations(network, t, z, w_0, lognormal_sigma, segment):
    results = []
    n_segs = network.n_segments

    # mean and variance of underlying normal distribution
    mu = log(w_0 / sqrt(1 + lognormal_sigma**2))
    var = log(1 + lognormal_sigma**2)

    for i in xrange(n_sims):
        # randomly sample apertures
        w = w_0 * np.random.lognormal(mean=mu, sigma=sqrt(var), size=n_segs)
        network.width = w

        # recalculate flow and calculate temperature
        network.calculate_flow(
            fluid, essential_bc, point_sources, correct=True)
        Theta = network.calculate_temperature(fluid, segment, [z], [t])
        results.append(Theta)

    return results


def plot_results(x, Y):
    # get percentile of results
    Y_50 = np.percentile(Y, 50, axis=1)
    Y_25 = np.percentile(Y, 25, axis=1)
    Y_75 = np.percentile(Y, 75, axis=1)

    # create lower and upper bar range
    yerr = np.zeros((2, len(Y_50)))
    yerr[0, :] = Y_50 - Y_25
    yerr[1, :] = Y_75 - Y_50

    # plot results
    f = plt.figure()
    plt.errorbar(x, Y_50, yerr=yerr, marker='o', linewidth=2.0)
    plt.ylim((0.68, 0.81))
    plt.xlim((0, 0.425))
    plt.ylabel('$\Theta_p$ (-)')
    plt.xlabel('$\sigma_w / w_0$ (-)')
    f.show()

if __name__ == '__main__':
    start_time = time.time()
    np.random.seed(0)

    # fluid properties
    cp_w = 4300.0
    rho_w = 1000.0
    mu_w = 1E-3

    # reservoir properties
    k_r = 2.9
    cp_r = 1050.0
    rho_r = 2700.0
    alpha_r = k_r / (rho_r * cp_r)

    # operational properties
    m_inj = 50.0
    P_inj = 0.0
    t_end = 86400 * 365.25 * 17.33

    # network properties
    conn = [(0, 1), (1, 2), (2, 3), (1, 4), (2, 5), (3, 6), (4, 5), (5, 6),
            (4, 7), (5, 8), (6, 9), (7, 8), (8, 9), (9, 10)]
    n_segs = len(conn)
    L = 250 * np.ones(n_segs)
    L[0] = 100
    L[-1] = 100
    H = 500 * np.ones(n_segs)
    w_0 = 1E-3
    w = w_0 * np.ones(n_segs)
    n_inj = 0
    n_prod = 10

    essential_bc = {n_inj: P_inj}
    point_sources = {n_prod: -m_inj}

    # simulation parameters
    n_sims = 10000
    sigma_over_w_0 = np.linspace(0.01, 0.4, 40)
    Theta_results = np.zeros((len(sigma_over_w_0), n_sims))

    # create network object
    fluid = Fluid(density=rho_w, viscosity=mu_w, heat_capacity=cp_w)
    network = FractureNetworkThermal(conn, L, H, w, k_r, alpha_r)

    # run simulations and plot results
    for i, sigma in enumerate(sigma_over_w_0):
        Theta_results[i, :] = run_simulations(network, t_end, L[-1], w_0,
                                              sigma, 13)
    plot_results(sigma_over_w_0, Theta_results)

    print "Elapsed time: {0:.4f} seconds".format(time.time() - start_time)
