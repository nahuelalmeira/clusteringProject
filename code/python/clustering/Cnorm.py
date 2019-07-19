import os
import sys
import numpy as np
import igraph as ig
import networkx as nx
import pickle
from clustering_methods import compute_C_values

overwrite = False
if 'overwrite' in sys.argv:
    overwrite = True

networks_dir = '../../../networks/real'
data_dir = '../../../data/clustering'

networks = os.listdir(networks_dir)

for network in networks:
    if '_simple_gcc' in network:
        networks.remove(network[:-11])

pickle_file = 'Cnorm_data.pickle'
full_pickle_file = os.path.join(data_dir, pickle_file)
if os.path.isfile(full_pickle_file) and not overwrite:
    with open(full_pickle_file, 'rb') as f:
        data = pickle.load(f)
else:
    data = {}

networks.remove('dolphins')
for network in networks:
    print(network)
    if network in data:
        continue
    full_input_file = os.path.join(networks_dir, network, network + '.txt')
    #g = nx.read_edgelist(full_input_file)
    g = ig.Graph().Read_Edgelist(full_input_file, directed=False)
    C_values = compute_C_values(g, samples=10, package='ig', greedy_max_it=10000, verbose=True)
    data[network] = C_values
    print(C_values)
    with open(full_pickle_file, 'wb') as f:
        pickle.dump(data, f)

print(data)