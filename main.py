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

        self.level = int(math.pow(10, math.ceil(math.log(max(self.vertexes), 10))))

        self.a_G = self.get_auxiliary_graph()
        self.a_stochastic_matrix = self.get_auxiliary_stochastic_matrix()

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
        edges = []

        for e in self.edges:
            for i in range(self.k):
                edges.append((i * self.level + e[0], (i + 1) % self.k * self.level + e[1]))

        a_G = nx.DiGraph(edges)
        #self.plot(a_G)
        return a_G

    def get_edge(self, l, i, j):
        return l * self.level + i, (l + 1) % self.k * self.level + j

    def get_auxiliary_stochastic_matrix(self):
        st_matrix = {}

        for l in range(1, self.k):
            for i in self.vertexes:
                for j in self.vertexes:
                    st_matrix[(self.get_edge(l, i, j))] = self.original_stochastic_matrix.get((i, j), 0.0)

        for i in self.vertexes:
            for j in self.vertexes:
                st_matrix[(self.get_edge(0, i, j))] = self.stochastic_matrix.get((i, j), 0.0)

        for sink in self.sinks:
            st_matrix[(sink, -1)] = self.stochastic_matrix[(sink, -1)]

        st_matrix[(-1, -1)] = 1.0
        return st_matrix

    def count_final_probabilities(self, time=10):
        st_matrix = self.a_stochastic_matrix.copy()
        current_matrix = {}

        vertexes = set(np.array(list(self.a_stochastic_matrix.keys())).flat)
        list(vertexes).sort()

        for _ in range(time):
            for i in vertexes:
                for j in vertexes:
                    current_matrix[(i, j)] = 0.0
                    for k in vertexes:
                        current_matrix[(i, j)] += st_matrix.get((i, k), 0.0) * st_matrix.get((k, j), 0.0)

            st_matrix = current_matrix.copy()

        return st_matrix

    def get_results(self, short=True):
        final = self.count_final_probabilities()
        items = list(final.items())
        items.sort()

        for item in items:
            if not short or (item[0][1] == -1 and -1 < item[0][0] < self.level):
                print('P{0} = {1}'.format(item[0], item[1]))

    def plot(self, graph):
        nx.draw_networkx(graph, node_color='aqua')
        plt.show()


gr = NSRGraph()
gr.get_results()