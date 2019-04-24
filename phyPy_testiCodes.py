import networkx as nx
import numpy as np
import community as louvain
from networkx.algorithms import community as nxcommunity
import scipy as sp
from sklearn.metrics.pairwise import cosine_similarity


def split_subgraph_into_chunks_randomly(node_set, max_size=200):
    "Split the subgraph into chunks for faster processing. Return chunks."
    chunks_count = 1
    if len(node_set) > max_size:
        chunks_count = 1 + int(len(node_set) / 200)
    ys = list(node_set)
    random.shuffle(ys)
    size, leftover = divmod(len(node_set), chunks_count)
    chunks = [ys[0 + size * i: size * (i + 1)] for i in list(range(chunks_count))]
    edge = size * chunks_count
    for i in list(range(leftover)):
        chunks[i % chunks_count].append(ys[edge + i])
    chunk_sets = [set() for _ in range(len(chunks))]
    for i, c in zip(range(len(chunks)), chunks):
        chunk_sets[i].update(set(c))
    return chunk_sets


def community_detection_cosine_of_squared(g, node_set):
    """Square the adjacency matrix and then use cosine similarity to detect communities. Return communities."""
    import scipy as sp
    import numpy as np
    from sklearn.metrics.pairwise import cosine_similarity

    communities = []
    if len(node_set) > 1:
        adj_array = nx.adjacency_matrix(g.subgraph(node_set)).toarray()
        new_adj = np.multiply(
            cosine_similarity(
                sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) >= 0.75, adj_array)
        edges_to_remove = np.argwhere(new_adj != adj_array)
        barcode_dict = dict(zip(range(len(node_set)), list(node_set)))
        edges_to_remove_barcode = [(barcode_dict[i], barcode_dict[j])
                                   for i, j in edges_to_remove]
        sub_graph_copy = nx.Graph(g.subgraph(node_set))
        sub_graph_copy.remove_edges_from(edges_to_remove_barcode)
        cos_components = list(nx.connected_components(sub_graph_copy))
        cos_components.sort(key=len, reverse=True)
        for com in cos_components:
            if len(com) > 1:
                communities.append(com)
    return communities


def community_detection_louvain(g, node_set, init_communities=False):
    """Apply Louvain community detection on a single component. Return communities."""
    import community as louvain

    if len(node_set) > 1:
        if not init_communities:
            partition = louvain.best_partition(g.subgraph(node_set))
        else:
            partition = louvain.best_partition(g.subgraph(node_set), init_communities)
    #     communities = []
    #     for com in set(partition.values()):
    #         community_nodes = {nodes for nodes in partition.keys() if partition[nodes] == com}
    #         if len(community_nodes) > 1:
    #             communities.append(community_nodes)
    # return communities
    return [{nodes for nodes in partition.keys() if partition[nodes] == com} for com in set(partition.values())]


def community_detection_k_clique(g, node_set, k=3):
    """Apply k-clique community detection. Return communities."""
    from networkx.algorithms import community as nxcommunity

    if len(node_set) > 1:
        return [set(i) for i in nxcommunity.k_clique_communities(g.subgraph(node_set), k)]
    return []


def community_detection_biconnected_components(g, node_set):
    """Separate bi-connected components. Return components."""
    cut_vertices = set(nx.articulation_points(g.subgraph(node_set)))
    components = list(nx.connected_components(g.subgraph(node_set - cut_vertices)))
    components.sort(key=len, reverse=True)
    return components


def split_subgraph_into_chunks_randomly(node_set, max_size=200):
    "Split the subgraph into chunks for faster processing. Return chunks."
    chunks_count = 1
    if len(node_set) > max_size:
        chunks_count = 1 + int(len(node_set) / 200)
    node_list = list(node_set)
    random.shuffle(node_list)
    size, leftover = divmod(len(node_set), chunks_count)
    chunks = [node_list[0 + size * i: size * (i + 1)] for i in list(range(chunks_count))]
    edge = size * chunks_count
    for i in list(range(leftover)):
        chunks[i % chunks_count].append(node_list[edge + i])
    chunk_sets = [set() for _ in range(len(chunks))]
    for i, c in zip(range(len(chunks)), chunks):
        chunk_sets[i].update(set(c))
    return chunk_sets

def merge_communities(g, communities, node_set=0, strategy=1, mode=1):
    """Merge communities if appropriate. """
    if len(communities) == 1:
        return communities
    if strategy == 1:  # Merge ad-hoc
        merge_network = nx.Graph()
        for i, j in enumerate(communities):
            merge_network.add_node(i)
        for i, com1 in enumerate(communities):
            for k, com2 in enumerate(communities):
                if i < k:
                    if mode == 1:  # disjoint input communities.
                        if nx.number_of_edges(
                                g.subgraph(com1.union(com2))) - \
                                nx.number_of_edges(g.subgraph(com1)) - \
                                nx.number_of_edges(g.subgraph(com2)) > 8:
                            merge_network.add_edge(i, k)
                    else:  # overlapping input communities.
                        if nx.number_of_edges(
                                g.subgraph(com1.union(com2))) - \
                                len(set(g.subgraph(com1).edges()).union(
                                    set(g.subgraph(com2).edges()))) \
                                > 8:
                            merge_network.add_edge(i, k)
        # res = []
        # for i in list(nx.connected_components(merge_network)):
        #     subset = set()
        #     for j in list(i):
        #         for barcode in list(communities[j]):
        #             subset.add(barcode)
        #     res.append(subset)
        # return res
        return [{barcode for j in i for barcode in communities[j]}
                for i in nx.connected_components(merge_network)]
    # Merge by Initializing Louvain with the communities
    return community_detection_louvain(g, node_set, communities)

communities = []
for bi_connected_component in community_detection_biconnected_components(g, set(g.neighbors(u))):
    communities += [merge_communities(g, sub_communities, bi_connected_component, strategy=1)
                    for chunk in split_subgraph_into_chunks_randomly(bi_connected_component, max_size=50)
                    for sub_communities in community_detection_k_clique(g, chunk, 3)]
communities = []
for bi_connected_component in community_detection_biconnected_components(g, set(g.neighbors(u))):
    sub_communities = [community_detection_k_clique(g, chunk, 3)
                       for chunk in split_subgraph_into_chunks_randomly(bi_connected_component, max_size=50)]
    print(sub_communities)
    communities += merge_communities(g, sub_communities[0], bi_connected_component, strategy=1)

communities = []
for bi_connected_component in community_detection_biconnected_components(g, set(g.neighbors(u))):
    sub_communities = [cluster
                       for chunk in split_subgraph_into_chunks_randomly(bi_connected_component, max_size=50)
                       for cluster in community_detection_k_clique(g, chunk, 3)]
    print(sub_communities)
    communities += merge_communities(g, sub_communities, bi_connected_component, strategy=1)

communities = []
for bi_connected_component in community_detection_biconnected_components(g, set(g.neighbors(u))):
    communities += merge_communities(g, [cluster
                       for chunk in split_subgraph_into_chunks_randomly(bi_connected_component, max_size=50)
                       for cluster in community_detection_k_clique(g, chunk, 3)], bi_connected_component, strategy=1)

communities = [merge_communities(g, [cluster
                       for chunk in split_subgraph_into_chunks_randomly(bi_connected_component, max_size=50)
                       for cluster in community_detection_k_clique(g, chunk, 3)], bi_connected_component, strategy=1)
               for bi_connected_component in community_detection_biconnected_components(g, set(g.neighbors(u)))]

communities = []
for bi_connected_component in community_detection_biconnected_components(g, set(g.neighbors(u))):
    sub_communities = []
    for chunk in split_subgraph_into_chunks_randomly(bi_connected_component, max_size=50):
        sub_communities += community_detection_k_clique(g, chunk, 3)
        # for chunk_community in Physlr.community_detection_k_clique(g, chunk, 3):
        #     sub_communities.append(chunk_community)
    print(sub_communities)
    communities += merge_communities(g, sub_communities, bi_connected_component, strategy=1)

def main(self):
    adj_numpy = np.array([
        [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1],
        [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0]])
    sum(abs(adj_numpy - adj_numpy.T))
    u = 12
    g = nx.from_numpy_matrix(adj_numpy)
    sub_graph = g.subgraph(g.neighbors(u))
    node_set = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13}
    node_set = set(g.neighbors(12))
    a = community_detection_cosine_of_squared(g, node_set)
    b = community_detection_louvain(g, node_set)
    c = community_detection_k_clique(g, node_set)
    d = community_detection_biconnected_components(g, node_set)

    address="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/TCGGGCAAGAGCACCA.tsv"

    eg_list = [(1,2)]


merge_network = nx.Graph()
for i in enumerate(communities):
    merge_network.add_node(i)


def wrap_up(mst, messages, sender, receiver):
    neighbors_to_wrapup = list(set(mst.neighbors(sender)) - {receiver})
    if not neighbors_to_wrapup:
        messages[(receiver, sender)] = 1
    messages[(receiver, sender)] = 1 + sum([messages[(sender, neighbor)] for neighbor in neighbors_to_wrapup])



dfs = list(nx.dfs_edges(sub_mst))
stack = [dfs[0][0]]
messages = dict()
# Gather
for edge in dfs:
    # edge[0]
    # edge[1]
    while stack[-1] != edge[0]:
        wrap_up(sub_mst, messages, stack.pop(), stack[-1])
    stack.append(edge[1])
while len(stack) != 1:
    wrap_up(sub_mst, messages, stack.pop(), stack[-1])
# Distribute
for edge in dfs:
    wrap_up(sub_mst, messages, edge[0], edge[1])
#stack.pop()


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
