from .dfn import FractureNetworkFlow
import networkx as nx

# create the multi-direction graph from the flow
# the usage should be:
# network = FractureNetworkThermal(conn, L, H, w)
# network.calculate_flow(fluid, essential_bc, point_sources)
# network.calculate_temperature()
# construct_graph


class FractureNetworkThermal(FractureNetworkFlow):

    def construct_graph(self):
        self.graph = nx.MultiDiGraph()

        for seg in conn:
            # invert node order if flow is negative
            if self.mass_flow < 0:
                self.graph.add_edge(seg[::-1])
            else:
                self.graph.add_edge(seg)

        self.graph.add_edges_from(conn)

    def __flow_split(self):
        pass

    def calculate_temperature(self):
        pass
