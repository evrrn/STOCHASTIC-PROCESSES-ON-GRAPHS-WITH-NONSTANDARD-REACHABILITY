import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from NSRGraph import NSRGraph


dir = "input_2"


def plot_bars(x, z, steps):
    x = x[:]
    x[-1] = 7
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    yrange = range(steps)
    for y in yrange:
        ax.bar(x, z[y], y, zdir='y', alpha=0.5)

    ax.set_yticks(yrange)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()


def plot_sink(quantity):
    plt.plot(quantity)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()


gr = NSRGraph.from_files()
# gr.plot(gr.G)
# gr.plot(gr.a_G)

steps = 4
vertexes = [1, 2, 3, 4, 5, 6, -1]

for _ in range(1):
    quantity = 20
    temp = []
    for i in range(1, 6):
        q = random.randint(0, quantity)
        temp.append(q)
        quantity -= q
    temp.append(quantity)
    random.shuffle(temp)

    for v in range(1, 7):
        gr.dist[v] = temp[v - 1]

    print(gr.dist)

    states = [gr.get_next_distribution(i) for i in range(steps)]
    plot_bars(vertexes, [[states[i][v] for v in vertexes] for i in range(steps)], steps)

    quantity = [gr.get_next_distribution(i)[-1] for i in range(steps * 3)]
    plot_sink(quantity)
