import random
from NSRGraph import NSRGraph

dir = "input_2"

gr = NSRGraph.from_files()
#gr.plot(gr.G)

for _ in range(3):

    gr.plot(gr.a_G)
    quantity = 20
    dist = {}

    temp = []
    for i in range(1, 6):
        q = random.randint(0, quantity)
        temp.append(q)
        quantity -= q
    temp.append(quantity)
    random.shuffle(temp)

    gr.dist = {v: temp[v - 1] for v in range(1, 7)}

    print(gr.dist)
    print(gr.get_next_distribution(4))

