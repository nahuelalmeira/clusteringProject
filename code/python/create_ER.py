import os
import sys
import pathlib
import igraph as ig

N = int(sys.argv[1])
k = float(sys.argv[2])
min_seed = int(sys.argv[3])
max_seed = int(sys.argv[4])

overwrite = False
if 'overwrite' in sys.argv:
    overwrite = True

seeds = range(min_seed, max_seed)
p = k/N

model_nets_dir = '../../networks/model'
ER_dir = 'ER'
ER_dir_k = ER_dir + '_k{:.2f}'.format(k)
ER_dir_k_N = ER_dir_k + '_N{:d}'.format(N)
base_output_dir = os.path.join(model_nets_dir, ER_dir, ER_dir_k, ER_dir_k_N)

for seed in seeds:
    print(seed)
    network = ER_dir_k_N + '_{:05}'.format(seed)
    full_output_dir = os.path.join(base_output_dir, network)
    full_output_path = os.path.join(full_output_dir, network + '.txt')
    pathlib.Path(full_output_dir).mkdir(parents=True, exist_ok=True)
    g = ig.Graph().Erdos_Renyi(N, p)
    g.write_edgelist(full_output_path)

    gcc = g.components(mode='weak').giant()
    network = ER_dir_k_N + '_{:05}'.format(seed) + '_gcc'
    full_output_dir = os.path.join(base_output_dir, network)
    full_output_path = os.path.join(full_output_dir, network + '.txt')
    pathlib.Path(full_output_dir).mkdir(parents=True, exist_ok=True)
    gcc.write_edgelist(full_output_path)
