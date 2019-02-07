import networkx as nx
import numpy as np
import scipy as sp
from sklearn.metrics.pairwise import cosine_similarity
scipy.linalg.blas.sgemm(...) #matrix-matrix
scipy.linalg.blas.sgemv #matrix-vector

if globalNum == 1 & t= False:
    print()

adj_numpy  = np.array([[0,6,0,0,0],[6,0,0,0,0],[0,0,0,3,0],[0,0,3,0,0],[0,0,0,0,0]])

adj_numpy  = np.array([[0,6,0,0,0],[6,0,0,0,0],[0,0,0,3,0],[0,0,3,0,0],[0,0,0,0,0]])

adj_numpy  = np.array([[0,6,0,0,4],[6,0,0,0,3],[0,0,0,3,3],[0,0,3,0,3],[4,3,3,3,0]])


adj_numpy  = np.array([[0,6,4,0,0,0,0,1],[6,0,5,0,0,0,0,1],[4,5,0,0,0,0,0,1],[0,0,0,0,7,3,1,1],[0,0,0,7,0,4,0,1],[0,0,0,3,4,0,0,1],[4,0,0,1,0,0,0,1],[1,1,1,1,1,1,1,0]])

adj_numpy  = np.array([[0,6,4,2,0,0,0,1],[6,0,5,0,0,0,0,1],[4,5,0,0,0,0,0,1],[2,0,0,0,7,3,0,1],[0,0,0,7,0,4,0,1],[0,0,0,3,4,0,0,1],[0,0,0,0,0,0,0,1],[1,1,1,1,1,1,1,0]])

G = nx.from_numpy_matrix(adj_numpy)
g = G
sub_graph = g.subgraph(g.neighbors(7))
adj_array = nx.adjacency_matrix(sub_graph).toarray()
edges_to_remove = np.argwhere(cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) <= egRem_threshold)
sub_graph = nx.Graph(sub_graph)
sub_graph.remove_edges_from(edges_to_remove)
sub_graph = nx.freeze(sub_graph)
cos_components = list(nx.connected_components(sub_graph))
cos_components.sort(key=len, reverse=True)
multi_node_components = [i for i in cos_components if len(i) > 1]
return u, {v: i for i, vs in enumerate(multi_node_components) for v in vs}


components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)).union(set([u])))))
{v: i for i, vs in enumerate(components2) for v in vs}
u, {v: i for i, vs in enumerate(components) for v in vs}
u, {v: i for i, vs in enumerate(components) for v in vs if vs!=u}
for i, vs in enumerate(components)
u, {v: i for i, vs in enumerate(components) for v in vs if v!=u}

cos = cosine_similarity(adj.dot(adj))
new_adj = np.multiply((cos > 0.8), adj.toarray())

edges_to_remove = np.where( new_adj != adj.toarray() )

edges_to_remove = np.where( new_adj != adj.toarray() )

a = np.matrix('1 2 3; 5 4 2; 6 2 1')
np.all(a.dot(a) , np.dot(a,a))

sg = g.subgraph(set(g.neighbors(u)).union(set([u])))
nx.adjacency_matrix(sg).toarray()
sg_adj = nx.adjacency_matrix(sg)
sg_adj.toarray()
sg_adj = nx.adjacency_matrix(G)
set(g.neighbors(u)).union(set([u]))



####### NEW CODE TO MIX CURR + sqCOS

	    sub_graph = g.subgraph(g.neighbors(u))
            nodes_count = len(sub_graph)
            edges_count = sub_graph.number_of_edges()
            if edges_count == 0 or nodes_count == 0:
                components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)).union(set([u])))))
                return u, {v: i for i, vs in enumerate(components) for v in vs if v != u}
            # [ current
            cut_vertices = set(nx.articulation_points(sub_graph))
            components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)) - cut_vertices)))
            components.sort(key=len, reverse=True)
            multicomp_1 = [i for i in components if len(i) > 1]
            # current end ]
	    for sc in multicomp_1:
	        subcomp = sub_graph.subgraph(set(sc))
	        print(nx.adjacency_matrix(subcomp).toarray())
                adj_array = nx.adjacency_matrix(subcomp).toarray()
                sq_cos = cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array))
                edges_to_remove = np.argwhere(sq_cos < egRem_threshold)
                subcomp2 = nx.Graph(subcomp)
                subcomp2.remove_edges_from(edges_to_remove)
                subcomp2 = nx.freeze(subcomp2)
                cos_components = list(nx.connected_components(subcomp2))
                cos_components.sort(key=len, reverse=True)
                multi_node_components = [i for i in cos_components if len(i) > 1]
                if len(multi_node_components_temp) == len(multi_node_components):
                    global globalNum
                    globalNum = globalNum + 1
                return u, {v: i for i, vs in enumerate(multi_node_components) for v in vs}


