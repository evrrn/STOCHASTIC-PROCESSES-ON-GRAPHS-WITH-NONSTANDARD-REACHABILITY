import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import math


class NSRGraph:
    def __init__(self):
        self.G = nx.read_edgelist("input/graph.txt", create_using=nx.DiGraph)

        self.edges = [(int(edge[0]), int(edge[1])) for edge in self.G.edges]
        self.vertexes = set(np.array(self.edges).flat)

        self.n = len(self.vertexes)
        self.m = len(self.edges)

        self.stochastic_matrix = self.get_stochastic_matrix()
        self.original_stochastic_matrix = self.stochastic_matrix.copy()
        self.sinks = self.get_sinks()
        self.apply_sinks()

        self.k = int(open("input/k.txt").readline())
        self.distribution = self.get_distribution()

        self.aG = self.get_auxiliary_graph()

    def get_stochastic_matrix(self):
        probabilities = open("input/st_matrix.txt").readlines()
        st_matrix = {self.edges[i]: float(probabilities[i]) for i in range(self.m)}
        return st_matrix

    def get_sinks(self):
        pairs = [line.split() for line in open("input/sinks.txt").readlines()]
        sinks = {int(pair[0]): float(pair[1]) for pair in pairs}
        return sinks

    def apply_sinks(self):
        for sink in self.sinks:
            for j in self.vertexes:
                p = self.stochastic_matrix.get((sink, j), 0.0)
                if p:
                    self.stochastic_matrix[(sink, j)] *= (1 - self.sinks[sink])
                    self.stochastic_matrix[(sink, -1)] = self.sinks[sink]

        self.stochastic_matrix[(-1, -1)] = 1.0

    def get_distribution(self):
        pairs = [line.split() for line in open("input/dist.txt").readlines()]
        dist = {int(pair[0]): int(pair[1]) for pair in pairs}
        return dist

    def get_auxiliary_graph(self):
        level = int(math.pow(10, math.ceil(math.log(max(self.vertexes), 10))))
        edges = []

        for e in self.edges:
            for i in range(self.k):
                edges.append((i * level + e[0], (i + 1) % self.k * level + e[1]))

        aG = nx.DiGraph(edges)
        self.plot(aG)
        return aG

    def plot(self, graph):
        nx.draw_networkx(graph, node_color='aqua')
        plt.show()


gr = NSRGraph()
