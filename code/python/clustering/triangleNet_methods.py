import igraph as ig
import numpy as np
from itertools import combinations

def getEdgeMultiplicity(G):
   """
   Compute edge multiplicity of graph G and return
   it as a list.
   """
   m_values = []
   for e in G.es():
       source, target = e.tuple
       nn_source = set(G.vs[source].neighbors())
       nn_target = set(G.vs[target].neighbors())
       nn_intersection = nn_source.intersection(nn_target)
       m = len(nn_intersection)
       m_values.append(m)
   return m_values

def triangleGraph(G):
   """
   Build a "Triangle graph", whose nodes are the triangles
   of G and whose links correspond to triangles in G that
   share an edge.
   """
   edge_list = set([])
   tr_set = set([])
   for e in G.es():
       source, target = e.tuple
       nn_source = set(G.neighbors(source))
       nn_target = set(G.neighbors(target))
       nn_intersection = nn_source.intersection(nn_target)
       if nn_intersection:
           local_trs = []
           for neighbor in nn_intersection:
               tr = '-'.join(map(str, sorted([source, target, neighbor])))
               tr_set.add(tr)
               local_trs.append(tr)
           edge_list = edge_list.union(set(combinations(local_trs, 2)))
   d = {edge: i for i, edge in enumerate(tr_set)}
   edge_list_idx = []
   for edge in edge_list:
       edge_list_idx.append((d[edge[0]], d[edge[1]]))
   N = len(tr_set)
   G = ig.Graph()
   G.add_vertices(range(N))
   G.add_edges(edge_list_idx)
   return G