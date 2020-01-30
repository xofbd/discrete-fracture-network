import json

from dfn import (Fluid, FractureNetwork, FractureNetworkFlow,
                 FractureNetworkThermal)


def read_json(func):
    """Decorator for passing allowing reading data from JSON."""
    def wrapper(path_to_file):
        with open(path_to_file, 'r') as f:
            data = json.load(f)

        return func(data)
    return wrapper


@read_json
def read_fluid_json(data):
    return Fluid(**data)


@read_json
def read_network_json(data):
    return FractureNetwork(**data)


@read_json
def read_flow_json(data):
    return FractureNetworkFlow(**data)


@read_json
def read_thermal_json(data):
    return FractureNetworkThermal(**data)
