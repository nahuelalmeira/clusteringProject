import numpy as np
import igraph as ig
import networkx as nx

def clusteredGraph(degSeq):
   
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
    #print('Iteration', it)
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
            it += 1
            #print('Iteration', it)
            
            last_edge = edgelist.pop()
            forbidden.append(last_edge)
            s, t = edgelist[-1]
            for node in last_edge:
                currentDegSeq[node] -= 1
                
        else:
            
            if t < N-1:
                t += 1
            else:
                if s < N-2:
                    s += 1
                    t = s + 1
    
    return sorted(edgelist, key=lambda x: (x[0], x[1]))

def get_C_greedy(degSeq):

    edgelist = clusteredGraph(degSeq)
    g = nx.Graph()
    g.add_edges_from(edgelist)
    C_greedy = nx.transitivity(g)
    return C_greedy

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

def get_C_Havel_Hakimi(degSeq):
    g = nx.havel_hakimi_graph(degSeq)
    C_HH = nx.transitivity(g)
    return C_HH

