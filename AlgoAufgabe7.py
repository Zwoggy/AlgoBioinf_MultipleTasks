"""
Aufgabe 7: Implementieren Sie den de Bruijn-Algorithmus entsprechend des Original-Papers (siehe Anhang).

"""


import networkx as nx
import matplotlib.pyplot as plt
true_sequence = 'ATGGCGTGCAATGGCGT'
reads = ['ATGGCGT', 'GGCGTGC', 'CGTGCAA', 'TGCAATG', 'CAATGGC', 'ATGGCGT']


def get_debruijn_edges_from_kmers(kmers):
    """
    Every possible (k-1)mer (n-1 suffix and prefix of kmers) is assigned
    to a node, and we connect one node to another if the (k-1)mer overlaps
    another. Nodes are (k-1)mers, edges are kmers.
    """
    # store edges as tuples in a set
    edges = set()
    k_mere_length_should_be = 5
    len_kmeres = len(kmers[0]) - k_mere_length_should_be
    # compare each (k-1)mer
    for k1 in kmers:
        for k2 in kmers:
            if k1 != k2:
                # if they overlap then add to edges
                if k1[len_kmeres:] == k2[:-len_kmeres]:
                    edges.add((k1[:-len_kmeres], k2[:-len_kmeres]))
                if k1[:-len_kmeres] == k2[len_kmeres:]:
                    edges.add((k2[:-len_kmeres], k1[:-len_kmeres]))

    return edges


def build_debruijn_graph(edges):
    """
    Construct a De Bruijn graph.

    Arguments:
    edges (list): List of sequence reads

    Returns:
    G (networkx.DiGraph): A directed multigraph representing the De Bruijn graph
    """
    G = nx.MultiDiGraph()
    for edge in edges:
        for i in range(len(edge) - 1):
            e = edge[i:i + 2]
            if G.has_edge(e[0], e[1]):
                G[e[0]][e[1]][0]['weight'] += 1
            else:
                G.add_edge(e[0], e[1], weight = 1)

    return G


def plot_debruijn_graph(G):
    """
    Draw a De Bruijn graph.

    Arguments:
    G (networkx.DiGraph): A directed multigraph representing the De Bruijn graph
    """
    pos = nx.spring_layout(G)
    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels = labels)
    nx.draw(G, pos, with_labels = True, arrows = True)
    plt.show()


if __name__ == "__main__":
    edges = get_debruijn_edges_from_kmers(reads)
    G = build_debruijn_graph(edges)
    plot_debruijn_graph(G)

