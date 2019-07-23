import os
import sys
import pathlib
import igraph as ig

N = int(sys.argv[1])
m = int(sys.argv[2])
min_seed = int(sys.argv[3])
max_seed = int(sys.argv[4])

overwrite = False
if 'overwrite' in sys.argv:
    overwrite = True

seeds = range(min_seed, max_seed)

model_nets_dir = '../../networks/model'
BA_dir = 'BA'
BA_dir_m = BA_dir + '_m{:d}'.format(m)
BA_dir_m_N = BA_dir_m + '_N{:d}'.format(N)
base_output_dir = os.path.join(model_nets_dir, BA_dir, BA_dir_m, BA_dir_m_N)

for seed in seeds:
    print(seed)
    network = BA_dir_m_N + '_{:05}'.format(seed)
    full_output_dir = os.path.join(base_output_dir, network)
    full_output_path = os.path.join(full_output_dir, network + '.txt')
    pathlib.Path(full_output_dir).mkdir(parents=True, exist_ok=True)
    g = ig.Graph().Barabasi(N, m)
    g.write_edgelist(full_output_path)
