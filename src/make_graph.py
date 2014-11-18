import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

pairs = [("openmm", "yank"), ("cmake", "openmm"), ("cuda", "openmm"), ("fftw3f", "openmm"), ("swig", "openmm"), ("mpi4py", "yank"), ("netcdf4", "yank"),
("numpy", "mdtraj"), ("scipy", "mdtraj"), # ("mdtraj", "mixtape"), ("sklearn", "mixtape"), 
("sphinx-bibtex", "openmm"), ("mdtraj", "yank"), ("ambermini", "yank")]

pairs = [(b, a) for (a, b) in pairs]

all_nodes = list(np.unique(pairs))

hard_nodes_list = [
all_nodes, 
["cmake", "cuda", "sphinx-bibtex", "ambermini", "openmm", "mdtraj"],
[]
]




graph = nx.DiGraph(pairs)
positions = nx.spring_layout(graph)

for (k, hard_nodes) in enumerate(hard_nodes_list):
    plt.figure()
    easy_nodes = list(np.setdiff1d(all_nodes, hard_nodes))
    nx.draw(graph, pos=positions, alpha=0.0)
    nx.draw_networkx_nodes(graph, positions, nodelist=easy_nodes, node_color="b", alpha=0.5)
    nx.draw_networkx_nodes(graph, positions, nodelist=hard_nodes, node_color="r", alpha=0.5)

    nx.draw_networkx_edges(graph, positions, alpha=0.5)

    nx.draw_networkx_labels(graph, positions, font_size=18)
    
    plt.savefig("./figures/dependencies%d.png" % k)

