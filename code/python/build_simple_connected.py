import os
import igraph as ig
import pathlib
networks_dir = '../../networks/real'

networks = os.listdir(networks_dir)
print(networks)

original_networks = [net for net in networks if '_simple_gcc' not in net]
gcc_networks = networks = [net for net in networks if '_simple_gcc' in net]

for network in original_networks:

    print(network)

    network_simple_gcc = network + '_simple_gcc' 
    if network_simple_gcc in gcc_networks:
        continue

    full_input_file = os.path.join(networks_dir, network, network + '.txt')
    try:
        g = ig.Graph().Read_Edgelist(full_input_file, directed=False)
    except:
        print('    ERROR: Could not read network', network)
        continue

    save = False
    if not g.is_simple():
        print("    Simplifying")
        g.simplify()
        save = True
    if not g.is_connected():
        print("    Taking giant")
        g = g.components(mode='weak').giant()
        save = True

    if save:
        print('    Writing simple connected network')
        output_dir = os.path.join(networks_dir, network_simple_gcc)
        pathlib.Path(output_dir).mkdir(exist_ok=True)
        output_file = os.path.join(output_dir, network_simple_gcc + '.txt')
        g.write_edgelist(output_file)

    