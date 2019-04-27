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
    sub_mst = nx.maximum_spanning_tree(g.subgraph(node_set), weight="n")
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
    # stack.pop()


merge_network = nx.Graph()
for i in enumerate(communities):
    merge_network.add_node(i)


def wrap_up(mst, messages, sender, receiver):
    neighbors_to_wrapup = list(set(mst.neighbors(sender)) - {receiver})
    if not neighbors_to_wrapup:
        messages[(receiver, sender)] = 1
    messages[(receiver, sender)] = 1 + sum([messages[(sender, neighbor)] for neighbor in neighbors_to_wrapup])


    def detect_and_cut_junctions_of_tree(g, messages, pruning_threshold=50):
        """"
        detect the junction in the tree, and split the tree from that point
        """
        g = g.copy()
        set_of_all_junctions = {node
                                for node in list(g.nodes)
                                if g.degree(node) > 3}
        for candidate in set_of_all_junctions:
            candidate_messages = [messages[(candidate, e)] for e in g.neighbors(candidate)]
        for node_to_remove in set_of_nodes_for_pruning:
            g.remove_node(node_to_remove)
        return g
    

def prune_branches_of_trees(g, messages, pruning_threshold=20):
    """"Determine the backbones of the maximum spanning trees
            and remove branches smaller than branch_size."""
    g = g.copy()
    # list_of_edges_for_pruning = [(node, neighbor)
    #                              for node in list(g.nodes)
    #                              for neighbor in g.neighbors(node)
    #                              if messages[(node, neighbor)] < pruning_threshold]
    set_of_nodes_for_pruning = {neighbor
                                 for node in list(g.nodes)
                                 for neighbor in g.neighbors(node)
                                 if messages[(node, neighbor)] < pruning_threshold}
    #new_edge_list = [(i, j) for i, j in list_of_edges_for_pruning if i not in [k for t, k in list_of_edges_for_pruning]]
    #list_of_edges_for_pruning = [(i, j) for i, j in list_of_edges_for_pruning if i not in [k for t, k in list_of_edges_for_pruning]]
    print(set_of_nodes_for_pruning)
    for node_to_remove in set_of_nodes_for_pruning:
        g.remove_node(node_to_remove)
    # for edge in list_of_edges_for_pruning:
    #     stack = [edge[1]]
    #     while stack:
    #         node_to_delete = stack.pop()
    #         for new_node in g.neighbors(node_to_delete):
    #             if new_node != edge[0]:
    #                 stack.append(new_node)
    #         g.remove_node(node_to_delete)
    return g


pruned_sub_mst = prune_branches_of_trees(sub_mst, messages, 4)
nx.adj_matrix(pruned_sub_mst).toarray()