import json
import os
import unittest

from dfn import (Fluid, FractureNetwork, FractureNetworkFlow,
                 FractureNetworkThermal, read_network_json, read_flow_json,
                 read_fluid_json, read_thermal_json)

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class TestReadJSON(unittest.TestCase):
    """Test reading from JSON."""

    @staticmethod
    def get_file_path(filename):
        return os.path.join(DATA_DIR, filename)

    def test_read_fluid_json(self):
        file_path = self.get_file_path('fluid_data.json')
        fluid_1 = read_fluid_json(file_path)
        with open(file_path, 'r') as f:
            data = json.load(f)

        fluid_2 = Fluid(data['density'], data['viscosity'],
                        data['heat_capacity'])
        self.assertEqual(fluid_1, fluid_2)

    def test_read_network_json(self):
        file_path = self.get_file_path('network_data.json')
        network_1 = read_network_json(file_path)
        with open(file_path, 'r') as f:
            data = json.load(f)

        network_2 = FractureNetwork(data['connectivity'], data['length'],
                                    data['thickness'], data['width'])
        self.assertEqual(network_1, network_2)

    def test_read_flow_json(self):
        file_path = self.get_file_path('network_data.json')
        network_1 = read_flow_json(file_path)
        with open(file_path, 'r') as f:
            data = json.load(f)

        network_2 = FractureNetworkFlow(data['connectivity'], data['length'],
                                        data['thickness'], data['width'])
        self.assertEqual(network_1, network_2)

    def test_read_thermal_json(self):
        file_path = self.get_file_path('network_thermal_data.json')
        network_1 = read_thermal_json(file_path)
        with open(file_path, 'r') as f:
            data = json.load(f)

        network_2 = FractureNetworkThermal(data['connectivity'],
                                           data['length'],
                                           data['thickness'],
                                           data['width'],
                                           data['thermal_cond'],
                                           data['thermal_diff'])
        self.assertEqual(network_1, network_2)


if __name__ == '__main__':
    unittest.main()
