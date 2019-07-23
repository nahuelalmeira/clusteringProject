import os
import sys
import pandas as pd
import numpy as np
import pickle
import pathlib
import igraph as ig
from triangleNet_methods import triangleGraph

## Auxiliar methods
def protocolName(mode, init_mu, step, samples, transit, decorr, rampe):
    
    protocol = '{}_{}_muInit{:.6f}_step{:.6f}_samples{:06}_transit{:06d}_decorr{:06d}'.format(mode, rampe, init_mu, step, samples, transit, decorr)
    
    return protocol
get_mu = lambda name: float(name.split('mu')[1].split('_')[0])
get_sample = lambda name: int(name.split('sample')[1].split('.')[0])


## Netwoks to compute 
net_dir = '../../../networks/real'
networks = os.listdir(net_dir)
for network in networks:
    if '_simple_gcc' in network:
        networks.remove(network[:-11])

## Protocol description
mode = 'maxC'
min_mu = 0
max_mu = 10
step = 0.1
samples = 100
decorr = 10
transit = 100
seed = 0

rampes = ['annealing', 'cooling']
seed_dir = 'seed{:05}'.format(seed)

sizes = []
for network in networks:
    full_input_file = os.path.join(net_dir, network, network + '.txt')
    try:
        g = ig.Graph().Read_Edgelist(full_input_file, directed=False)
    except:
        print('ERROR: Could not read file', network)
        continue
    N = g.vcount()
    M = g.ecount()
    sizes.append([network, N, M])
    
networks = list(zip(*sorted(sizes, key=lambda x: x[2])))[0]
for network in networks:
    print(network)
    for rampe in rampes:

        if rampe == 'annealing':
            init_mu = min_mu
        else:
            init_mu = max_mu
        protocol = protocolName(mode, init_mu, step, samples, transit, decorr, rampe)
        input_dir = os.path.join(net_dir, network, protocol, seed_dir)
        if not os.path.isdir(input_dir):
            print(input_dir)

        output_dir = os.path.join(net_dir, network, 'tr_' + protocol, seed_dir)
        pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        files = sorted(os.listdir(input_dir))  
        for f in files:
            if 'mu' not in f:
                continue
            mu = get_mu(f)
            sample = get_sample(f)
            if sample == 1:
                print(mu)
            full_input_path = os.path.join(input_dir, f)
            full_output_path = os.path.join(output_dir, 'tr_' + f)
            if os.path.isfile(full_output_path):
                continue
            g = ig.Graph().Read_Edgelist(full_input_path, directed=False)
            tr_g = triangleGraph(g)
            tr_g.write_edgelist(full_output_path)