import sys

from dfn.base import Base
from dfn.Fluid import Fluid
from dfn.FractureNetwork import FractureNetwork
from dfn.FractureNetworkFlow import FractureNetworkFlow
from dfn.FractureNetworkThermal import FractureNetworkThermal
from dfn.read_json import (
        read_json,
        read_flow_json,
        read_fluid_json,
        read_network_json,
        read_thermal_json
    )


if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    import importlib_metadata as metadata

__all__ = [
        "Base",
        "Fluid",
        "FractureNetwork",
        "FractureNetworkFlow",
        "FractureNetworkThermal",
        "read_json",
        "read_flow_json",
        "read_fluid_json",
        "read_network_json",
        "read_thermal_json",
    ]

__version__ = metadata.version("dfn")
