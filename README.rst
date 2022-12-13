|Python| |Release| |License| |CI| |Codecov|

.. |Python| image:: https://shields.io/badge/Python-3.8%20%7C%203.9%20%7C%203.10-blue
.. |Release| image:: https://img.shields.io/github/v/release/xofbd/discrete-fracture-network.svg
   :target: https://github.com/xofbd/discrete-fracture-network.svg/releases
.. |License| image:: https://img.shields.io/github/license/xofbd/discrete-fracture-network.svg
   :target: https://opensource.org/licenses/MIT
.. |CI| image:: https://github.com/xofbd/discrete-fracture-network/actions/workflows/ci.yaml/badge.svg?branch=master
.. |Codecov| image:: https://codecov.io/gh/xofbd/discrete-fracture-network/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/xofbd/discrete-fracture-network

Discrete Fracture Network
=========================
An analytical thermohydraulic model for discretely fractured geothermal
reservoirs. This package is an implementation of the model described in
Fox, D. B., D. L. Koch, and J. W. Tester (2016), An analytical thermohydraulic
model for discretely fractured geothermal reservoirs, Water Resources Research,
52, 6792-6817, doi:10.1002/2016WR018666

Installation
----------
To install, run::

    pip install git+https://github.com/xofbd/discrete-fracture-network.git

or::

    git clone https://github.com/xofbd/discrete-fracture-network.git
    pip install discrete-fracture-network/

or::

    git clone https://github.com/xofbd/discrete-fracture-network.git
    cd discrete-fracture-network/
    make install

Note: you can use the ``--user`` flag for a local installation as opposed to system-wide.

Usage
-----
Here is a short example of creating a simple network.
  >>> import dfn
  >>> conn = [(0, 1), (1, 2), (1, 2), (2, 3)]
  >>> L = [100., 500., 500., 100.]
  >>> H = [500., 500., 500., 500.]
  >>> w = [1E-3, 1E-3, 1E-3, 1E-3]
  >>> network = dfn.FractureNetwork(conn, L, H, w)

For more extensive examples of usage, look at the scripts in the ``example`` directory.

Tests
-----
Unit tests are in ``tests`` and can be run with ``make test-unit``. To lint, run ``make test-lint``. Both unit testing and linting can be run with ``make tests``.

License
-------
This project is distributed under the MIT license. Please see ``LICENSE`` for more information
