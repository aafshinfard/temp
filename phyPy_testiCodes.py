import networkx as nx
import numpy as np
import community as louvain
from networkx.algorithms import community as nxcommunity
import scipy as sp
from sklearn.metrics.pairwise import cosine_similarity
scipy.linalg.blas.sgemm(...) #matrix-matrix
scipy.linalg.blas.sgemv #matrix-vector

if globalNum == 1 & t= False:
    print()

adj_numpy  = np.array([[0,6,0,0,0],[6,0,0,0,0],[0,0,0,3,0],[0,0,3,0,0],[0,0,0,0,0]])

adj_numpy  = np.array([[0,6,0,0,0],[6,0,0,0,0],[0,0,0,3,0],[0,0,3,0,0],[0,0,0,0,0]])

adj_numpy  = np.array([[0,6,0,0,4],[6,0,0,0,3],[0,0,0,3,3],[0,0,3,0,3],[4,3,3,3,0]])


adj_numpy  = np.array([[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1],[1,1,1,1,1,1,1,0]])
adj_numpy  = np.array([[0,6,4,0,0,0,0,1],[6,0,5,0,0,0,0,1],[4,5,0,0,0,0,0,1],[0,0,0,0,7,3,1,1],[0,0,0,7,0,4,0,1],[0,0,0,3,4,0,0,1],[4,0,0,1,0,0,0,1],[1,1,1,1,1,1,1,0]])

adj_numpy  = np.array([[0,6,4,2,0,0,1],[6,0,5,0,0,0,1],[4,5,0,0,0,0,1],[2,0,0,0,7,3,1],[0,0,0,7,0,4,1],[0,0,0,3,4,0,1],[1,1,1,1,1,1,0]])

adj_numpy  = np.array([[0,6,4,2,0,0,0,1],[6,0,5,0,0,0,0,1],[4,5,0,0,0,0,0,1],[2,0,0,0,7,3,0,1],[0,0,0,7,0,4,0,1],[0,0,0,3,4,0,0,1],[0,0,0,0,0,0,0,1],[1,1,1,1,1,1,1,0]])
adj_numpy  = np.array([
[0,1,1,1,0,0,0,0,0,0,0,0,1,0],
[1,0,1,1,0,0,0,0,0,0,0,0,1,0],
[1,1,0,1,1,0,0,0,0,0,0,0,1,0],
[1,1,1,0,0,1,0,0,0,0,0,0,1,0],
[0,0,1,0,0,1,1,1,0,0,0,0,1,0],
[0,0,0,1,1,0,1,1,0,0,0,0,1,1],
[0,0,0,0,1,1,0,1,0,0,0,0,1,0],
[0,0,0,0,1,1,1,0,0,0,0,0,1,1],
[0,0,0,0,0,0,0,0,0,1,0,1,1,1],
[0,0,0,0,0,0,0,0,1,0,1,1,1,1],
[0,0,0,0,0,0,0,0,0,1,0,1,1,0],
[0,0,0,0,0,0,0,0,1,1,1,0,1,0],
[1,1,1,1,1,1,1,1,1,1,1,1,0,1],
[0,0,0,0,0,1,0,1,1,1,0,0,1,0]])
sum(abs(adj_numpy-adj_numpy.T))

u = 12
threshold=0.75
G = nx.from_numpy_matrix(adj_numpy)
g = G
sub_graph = g.subgraph(g.neighbors(7))
adj_array = nx.adjacency_matrix(sub_graph).toarray()
cos = cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array))
edges_to_remove = np.argwhere(cos <= threshold)
sub_graph = nx.Graph(sub_graph)
sub_graph.remove_edges_from(edges_to_remove)
sub_graph = nx.freeze(sub_graph)
cos_components = list(nx.connected_components(sub_graph))
cos_components.sort(key=len, reverse=True)
multi_node_components = [i for i in cos_components if len(i) > 1]
return u, {v: i for i, vs in enumerate(multi_node_components) for v in vs}


for comp in components:
    if len(comp) > 1:
        sub_graph = g.subgraph(g.neighbors(comp))
adj_array = nx.adjacency_matrix(sub_graph).toarray()
new_adj = np.multiply(
            cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) >= 0.75, adj_array)
edges_to_remove = np.argwhere(new_adj != adj_array)
barcode_dict = dict(zip(range(len(sub_graph)), list(sub_graph.nodes)))
edges_to_remove_barcode = [(barcode_dict[i], barcode_dict[j]) for i, j in edges_to_remove]
sub_graph_copy = nx.Graph(sub_graph)
sub_graph_copy.remove_edges_from(edges_to_remove_barcode)
cos_components = list(nx.connected_components(sub_graph_copy))
cos_components.sort(key=len, reverse=True)
for com in cos_components:
            if len(com) > 1:
                communities.append(com)
        return

import community as cm
partition = cm.best_partition(sub_graph)
for com in set(partition.values()) :
    list_nodes = [nodes for nodes in partition.keys()
                                if partition[nodes] == com]
	if len(list_nodes)>1 :
		print({com: v for v in list_nodes})


from networkx.algorithms import community
a=community.greedy_modularity_communities(sub_graph)




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
            multicomp = [i for i in components if len(i) > 1]
            # current end ]
	    comul_multisubcomps =[]
	    for sc in multicomp:
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
                multisubcomp = [i for i in cos_components if len(i) > 1]
	        comul_multisubcomps = comul_multisubcomps + multisubcomp

            return u, {v: i for i, vs in enumerate(comul_multisubcomps) for v in vs}

### updated:
            sub_graph = g.subgraph(g.neighbors(u))
            nodes_count = len(sub_graph)
            # edges_count = sub_graph.number_of_edges()
            if nodes_count == 0:  # or edges_count == 0:
                #components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)).union(set([u])))))
                components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)) )))
                return u, {v: i for i, vs in enumerate(components) for v in vs if v != u}
            # adj_array = nx.adjacency_matrix(sub_graph).toarray()
            # sq_cos = cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array))
            # edges_to_remove = np.argwhere(sq_cos < threshold)
            
            adj_array = nx.adjacency_matrix(sub_graph).toarray()
            adj2 = adj_array
            adj2[np.argwhere(cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) < threshold)] = 0
            edges_to_remove = np.argwhere(adj2 != adj_array)
            edges_to_remove2 = [(list(sub_graph.nodes)[x], list(sub_graph.nodes)[y]) for x, y in edges_to_remove]
            sub_graph2 = nx.Graph(sub_graph)
            sub_graph2.remove_edges_from(edges_to_remove2)



            adj_array = nx.adjacency_matrix(sub_graph).toarray()
            adj2 = nx.adjacency_matrix(sub_graph).toarray()
            adj2[np.argwhere(cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) < threshold)] = 0
            edges_to_remove = np.argwhere(adj2 != adj_array)
            edges_to_remove2 = [(list(sub_graph.nodes)[x], list(sub_graph.nodes)[y]) for x, y in edges_to_remove]
            sub_graph2 = nx.Graph(sub_graph)
            sub_graph2.remove_edges_from(edges_to_remove2)
            
            adj_array = nx.adjacency_matrix(sub_graph).toarray()
            cos = cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array))
            row_names = list(x.nodes)
            column_names = list(x.nodes)
            cos_df = pd.DataFrame(cos, columns=column_names, index=row_names)
            edges_to_remove = cos_df[cos_df < threshold].stack().index.tolist()


	    # sub_graph2 = nx.freeze(sub_graph2)
            # removed = np.argwhere(nx.adjacency_matrix(sub_graph2).toarray() != adj.toarray())
            cos_components = list(nx.connected_components(sub_graph2))
            cos_components.sort(key=len, reverse=True)
            multi_node_components = [i for i in cos_components if len(i) > 1]
            #print(
            #    int(timeit.default_timer() - t0),
            #    "\nWorking", "with the one","with # comps equal to:", multi_node_components.__len__(),
            #    "\n with comps equal to:", multi_node_components,
            #    "\n and edges to remove:", edges_to_remove,
            #    "\n and # edges to remove:", edges_to_remove.__len__(),
            #    "\n and edges removed then:", removed,
            #    "\n cosine of squared:", cosine_similarity(np.dot(adj, adj)), "\nthreshold:", threshold,
            #    "\n and nodes list:", list(sub_graph2.nodes),
            #    file=sys.stderr)
            return u, {v: i for i, vs in enumerate(multi_node_components) for v in vs}




            adj_array = nx.adjacency_matrix(sub_graph).toarray()
            new_adj = np.multiply(cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) >= threshold , adj_array)
            row_names = list(sub_graph.nodes)
            column_names = row_names
            new_adj_df = pd.DataFrame(new_adj!=adj_array, columns=column_names, index=row_names)
            edges_to_remove2 = new_adj_df[new_adj_df == True].stack().index.tolist()


new_adj = np.multiply((cos > 0.8), adj.toarray())
