import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms as alg
import numpy as np
import math
from functools import reduce

dir = "input"


def find_gcd(list):
    x = reduce(math.gcd, list)
    return x


class NSRGraph:
    def __init__(self, edges, sinks, k):
        self.G = nx.parse_edgelist(edges, create_using=nx.DiGraph, nodetype=int, data=(('probability', float),))

        self.edges = [(edge[0], edge[1]) for edge in self.G.edges]
        self.vertexes = set(np.array(self.edges).flat)

        self.n = len(self.vertexes)
        self.m = len(self.edges)

        self.stochastic_matrix = self._get_stochastic_matrix()
        self.original_stochastic_matrix = self.stochastic_matrix.copy()
        self.sinks = sinks
        self._apply_sinks()

        self.k = k

        self.level = int(math.pow(10, math.ceil(math.log(max(self.vertexes), 10))))

        self.a_G = self._get_auxiliary_graph()
        self.a_stochastic_matrix = self._get_auxiliary_stochastic_matrix()

    @classmethod
    def from_files(cls):
        edges = open(dir + "/graph.txt").readlines()

        pairs = [line.split() for line in open(dir + "/sinks.txt").readlines()]
        sinks = {int(pair[0]): float(pair[1]) for pair in pairs}

        k = int(open(dir + "/k.txt").readline())

        return cls(edges, sinks, k)

    @classmethod
    def from_params(cls, edges, sinks, k):
        return cls(edges, sinks, k)

    def _get_stochastic_matrix(self):
        st_matrix = {(u, v): d['probability'] for (u, v, d) in self.G.edges(data=True)}
        return st_matrix

    def _apply_sinks(self):
        for sink in self.sinks:
            for j in self.vertexes:
                p = self.stochastic_matrix.get((sink, j), 0.0)
                if p:
                    self.stochastic_matrix[(sink, j)] *= (1 - self.sinks[sink])
                    self.stochastic_matrix[(sink, -1)] = self.sinks[sink]

        self.stochastic_matrix[(-1, -1)] = 1.0

    def get_distribution(self):
        pairs = [line.split() for line in open(dir + "/dist.txt").readlines()]
        dist = {int(pair[0]): int(pair[1]) for pair in pairs}
        return dist

    def _get_auxiliary_graph(self):
        edges = []

        for e in self.edges:
            for i in range(self.k):
                edges.append((i * self.level + e[0], (i + 1) % self.k * self.level + e[1]))

        a_G = nx.DiGraph(edges)
        return a_G

    def _get_edge(self, l, i, j):
        return l * self.level + i, (l + 1) % self.k * self.level + j

    def _get_auxiliary_stochastic_matrix(self):
        st_matrix = {}

        for l in range(1, self.k):
            for i in self.vertexes:
                for j in self.vertexes:
                    st_matrix[(self._get_edge(l, i, j))] = self.original_stochastic_matrix.get((i, j), 0.0)

        for i in self.vertexes:
            for j in self.vertexes:
                st_matrix[(self._get_edge(0, i, j))] = self.stochastic_matrix.get((i, j), 0.0)

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

    def _print_probability(self, v, p):
        print('P{0} = {1}'.format((v, -1), p))

    def count_experimental_results(self, short=True):
        final = self.count_final_probabilities()
        items = list(final.items())
        items.sort()

        results = {}

        for item in items:
            if not short or (item[0][1] == -1 and -1 < item[0][0] < self.level):
                results[item[0][0]] = item[1]

        return results

    def get_experimental_results(self):
        print('Experimental results:')
        results = self.count_experimental_results()

        for (v, p) in results.items():
            self._print_probability(v, p)

    def find_gcd_of_cycles_lengths_and_k(self):
        list = [self.k] + [len(c) for c in alg.cycles.cycle_basis(self.G.to_undirected())]

        for e in self.edges:
            if (e[1], e[0]) in self.edges:
                list += [2]
                break

        return find_gcd(list)

    def count_theoretical_results(self):
        if self.find_gcd_of_cycles_lengths_and_k() == 1:
            results = {v: 1.0 for v in self.vertexes}
        else:
            rest = self.vertexes.copy()

            for s in self.sinks:
                temp = rest.copy()
                paths = nx.shortest_path_length(self.G, target=s)
                for v in rest:
                    path_length = paths[v]
                    if path_length % self.k == 0:
                        temp.remove(v)
                rest = temp.copy()

            results = {v: float(v not in rest) for v in self.vertexes}

        return results

    def get_theoretical_results(self):
        if not alg.components.is_strongly_connected(self.G):
            print('Graph is not strongly connected')
            return

        print('Theoretical results:')
        results = self.count_theoretical_results()

        for (v, p) in results.items():
            self._print_probability(v, p)

    def plot(self, graph):
        nx.draw_networkx(graph, node_color='aqua')
        plt.show()


gr = NSRGraph.from_files()
gr.get_experimental_results()
print()
gr.get_theoretical_results()



