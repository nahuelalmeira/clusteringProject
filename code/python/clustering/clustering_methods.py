import numpy as np
import igraph as ig
import networkx as nx

def clusteredGraph(degSeq, maxIt=100, verbose=False):
   
    if not ig.is_graphical_degree_sequence(degSeq):
        print('Not a graphical degree sequence')
        return None
    
    degSeq = sorted(map(int, degSeq), reverse=True)
    #print(dispStubs)
    N = len(degSeq)
    
    currentDegSeq = [0]*N
    
    edgelist = []
    forbidden = []
    s = 0
    t = 1
    it = 1

    while True:
        
        edge = (s, t)
        
        available = True
        if currentDegSeq[s] == degSeq[s]:
            available = False     
        if currentDegSeq[t] == degSeq[t]:
            available = False
        if edge in forbidden:
            available = False
        
        if available:
            edgelist.append(edge)
            currentDegSeq[s] += 1
            currentDegSeq[t] += 1
            forbidden = []
            
        if s == N-2 and t == N-1:            
            if 2*len(edgelist) == sum(degSeq):
                break

            last_edge = edgelist.pop()
            forbidden.append(last_edge)
            s, t = edgelist[-1]
            for node in last_edge:
                currentDegSeq[node] -= 1
            it += 1
            #if verbose:
            #    print('Greedy, it', it)
               
        else:
            
            if t < N-1:
                t += 1
            else:
                if s < N-2:
                    s += 1
                    t = s + 1
        

        if it == maxIt:
            print('Max iterations reached:')
            print('    Expected M:', sum(degSeq)//2)
            print('    Obtained M:', len(edgelist))
            break

    if verbose:
        print('Greedy converged after', it, 'iterations.')

    return sorted(edgelist, key=lambda x: (x[0], x[1]))

def get_C_greedy(degSeq, maxIt=100, package='nx', verbose=False):

    edgelist = clusteredGraph(degSeq, maxIt=maxIt, verbose=verbose)
    if package == 'nx':
        g = nx.Graph()
        g.add_edges_from(edgelist)
        C_greedy = nx.transitivity(g)
        Cws_greedy = nx.average_clustering(g)
    elif package == 'ig':
        g = ig.Graph()
        N = np.max(edgelist) + 1
        g.add_vertices(N)
        g.add_edges(edgelist)
        C_greedy = g.transitivity_undirected(mode='zero')
        Cws_greedy = g.transitivity_avglocal_undirected(mode='zero')

    else:
        print('ERROR: Package "{}" not supported'.format(package))
    return C_greedy, Cws_greedy

def get_C_unc(degSeq):
    """
    Newman's formula for expected clustering coefficient
    of an uncorrelated network with given degree distribution.
    NOTE: For very heterogeneous network, it could give values
    greater than one.
    """
    N = len(degSeq)
    degSeq = np.array(degSeq)
    meank = np.mean(degSeq)
    meank2 = np.mean(degSeq**2)
    randC = (meank2 - meank)**2 / (N*meank**3)
    return randC

def get_C_rand(edgelist, swaps=100, n_samples=100):

    g = nx.Graph()
    g.add_edges_from(edgelist)
    M = g.number_of_edges()
    
    ### Perform swaps for randomizing
    nswap = swaps*M
    nx.double_edge_swap(g, nswap=nswap, max_tries=10*nswap)
    
    ## Compute values of C
    C_values = np.zeros(n_samples)
    C = nx.transitivity(g)
    C_values[0] = C
    nswap = int(nswap/n_samples)
    for i in range(1, n_samples):
        nx.double_edge_swap(g, nswap=nswap, max_tries=10*nswap)
        C = nx.transitivity(g)
        C_values[i] = C

    if n_samples == 1:
        return C
    return C_values

def get_C_rand_CM(degSeq, samples=10, package='nx', verbose=False):
    C_values = np.zeros(samples)
    Cws_values = np.zeros(samples)
    for i in range(samples):
        if verbose:
            print('rand CM, it', i)
        if package == 'nx':
            g = nx.configuration_model(degSeq, seed=i)
            g = nx.Graph(g)
            C = nx.transitivity(g)
            Cws = nx.average_clustering(g)
        elif package == 'ig':
            g = ig.Graph().Degree_Sequence(degSeq, method='vl')
            C = g.transitivity_undirected(mode='zero')
            Cws = g.transitivity_avglocal_undirected(mode='zero')
        else:
            print('ERROR: Package "{}" not supported'.format(package))
        C_values[i] = C
        Cws_values[i] = Cws

    return C_values, Cws_values

def get_C_Havel_Hakimi(degSeq):
    g = nx.havel_hakimi_graph(degSeq)
    C_HH = nx.transitivity(g)
    return C_HH

def compute_C_values(g, samples=10, package='nx', greedy_max_it=100, verbose=False):

    if package == 'nx':
        degSeq = list(dict(nx.degree(g)).values())
        C = nx.transitivity(g)
    elif package == 'ig':
        degSeq = g.degree()
        C = g.transitivity_undirected(mode='zero')
        Cws = g.transitivity_avglocal_undirected(mode='zero')
    else:
        print('ERROR: Package "{}" not supported'.format(package))
    
    if verbose:
        print('Computing C_greedy')
    C_greedy, Cws_greedy = get_C_greedy(degSeq, maxIt=greedy_max_it, package=package, verbose=verbose)
    C_rand, Cws_rand = get_C_rand_CM(degSeq, samples, package=package, verbose=verbose)
    
    C_norm = (C - C_rand.mean()) / (C_greedy - C_rand.mean())
    Cws_norm = (Cws - Cws_rand.mean()) / (Cws_greedy - Cws_rand.mean())
    
    return {'C': [C, C_rand.mean(), C_greedy, C_norm], 'Cws': [Cws, Cws_rand.mean(), Cws_greedy, Cws_norm]}