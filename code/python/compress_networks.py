import os
import sys
import tarfile

def protocolName(mode, init_mu, step, samples, transit, decorr, rampe):
    """
    Old version.
    """
    
    protocol = '{}_{}_muInit{:.6f}_step{:.6f}_samples{:06}_transit{:06d}_decorr{:06d}'.format(mode, rampe, init_mu, step, samples, transit, decorr)
    return protocol

def protocolNameNew(mode, init_mu, step, samples, transit, decorr, rampe, time_scaling):
    """
    Returns the name of the directory associated to a given MC protocol.
    """
    
    protocol = '{}_{}_muInit{:.6f}_step{:.6f}_samples{:06}_transit{:06d}_decorr{:06d}_timeScaling{:.6f}'.format(mode, rampe, init_mu, step, 
                                                                                                                samples, transit, decorr, time_scaling)
    return protocol

## Methods for getting the value of mu and sample number from the file names
get_mu = lambda name: float(name.split('mu')[1].split('_')[0])
get_sample = lambda name: int(name.split('sample')[1].split('.')[0])

model_nets_dir = '../../networks/model'

mode = 'maxC'
min_mu = 0
max_mu = 10
step = 0.1
samples = 100
seed = 0
time_scaling = 1

fast = False
decorr = 300
transit = 3000

m = 4
#alpha_values = np.arange(0, 1.1, 0.5)
alpha_values = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
if time_scaling > 1:
    seed = 1
else:
    seeds = 0
net_seed = 0
N_values = [500, 1000, 2000, 4000]

verbose = False
if 'verbose' in sys.argv:
    verbose = True

overwrite = False

HM_dir = 'HM'
HM_dir_m = 'HM_m{:d}'.format(m)
networks = []
network_dirs = []
for N in N_values:
    HM_dir_N = HM_dir_m + '_N{:d}'.format(N)
    for alpha in alpha_values:
        HM_dir_N_alpha = HM_dir_N + '_alpha{:03d}'.format(int(100*alpha))
        base_input_dir = os.path.join(model_nets_dir, HM_dir, HM_dir_m, HM_dir_N, HM_dir_N_alpha)
        network = HM_dir_N_alpha + '_{:05}_gcc'.format(net_seed)
        full_input_dir = os.path.join(base_input_dir, network)
        networks.append(network)
        network_dirs.append(full_input_dir)

rampes = ['annealing', 'cooling']
seed_dir = 'seed{:05}'.format(seed)

for network, net_dir in zip(networks, network_dirs):
    print(network)

    for rampe in rampes:

        if rampe == 'annealing':
            init_mu = min_mu
        else:
            init_mu = max_mu
        if time_scaling == 0:
            protocol = protocolName(mode, init_mu, step, samples, transit, decorr, rampe)
        else:
            protocol = protocolNameNew(mode, init_mu, step, samples, transit, decorr, rampe, time_scaling)
            
        input_dir = os.path.join(net_dir, protocol, seed_dir)
        if not os.path.isdir(input_dir):
            print(input_dir)
            continue

        files = sorted(os.listdir(input_dir))  
        data = []
        for input_name in files:
            if '.txt' not in input_name:
                continue
            
            net_name = input_name[:-4]
            full_input_name = os.path.join(input_dir, input_name)

            if verbose:
                print(input_name)

            ## Compress network file
            tar_input_name = net_name + '.tar.gz'
            full_tar_input_name = os.path.join(input_dir, tar_input_name)
            tar = tarfile.open(full_tar_input_name, 'w:gz')
            tar.add(full_input_name, arcname=input_name)
            tar.close()

            ## Remove network file
            os.remove(full_input_name)