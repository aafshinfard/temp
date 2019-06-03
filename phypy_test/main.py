import networkx as nx
import numpy as np
import community as louvain
from networkx.algorithms import community as nxcommunity
import scipy as sp
from sklearn.metrics.pairwise import cosine_similarity

mport argparse
import itertools
import multiprocessing
import os
import random
import re
import statistics
import sys
import timeit
from collections import Counter


setup = '''
from networkx.algorithms import community as nxcommunity
import networkx as nx
import numpy as np
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
'''
print(min(timeit.Timer('a = g[u].keys()', setup=setup).repeat(7, 1000000)))
print(min(timeit.Timer('a = set(g.neighbors(u))', setup=setup).repeat(7, 1000000)))




def detect_communities_biconnected_components(g, node_set):
    """Separate bi-connected components. Return components."""
    cut_vertices = set(nx.articulation_points(g.subgraph(node_set)))
    components = list(nx.connected_components(g.subgraph(node_set - cut_vertices)))
    return components


def detect_communities_k_clique(g, node_set, k=3):
    """Apply k-clique community detection. Return communities."""
    return list(nx.algorithms.community.k_clique_communities(g.subgraph(node_set), k))


def detect_communities_louvain(g, node_set, init_communities=None):
    """Apply Louvain community detection on a single component. Return communities."""
    if len(node_set) < 2:
        return []
    import community as louvain
    partition = louvain.best_partition(g.subgraph(node_set), init_communities)
    return [{node for node in partition.keys() if partition[node] == com}
            for com in set(partition.values())]


def detect_communities_cosine_of_squared(g, node_set, squaring=True, threshold=0.75):
    """
    Square the adjacency matrix and then use cosine similarity to detect communities.
    Return communities.
    """
    import scipy as sp
    import numpy as np
    from sklearn.metrics.pairwise import cosine_similarity

    communities = []
    if len(node_set) > 1:
        adj_array = nx.adjacency_matrix(g.subgraph(node_set)).toarray()
        if squaring:
            new_adj = np.multiply(
                cosine_similarity(
                    sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) >= threshold, adj_array)
        else:
            new_adj = np.multiply(cosine_similarity(adj_array) >= threshold, adj_array)
        edges_to_remove = np.argwhere(new_adj != adj_array)
        barcode_dict = dict(zip(range(len(node_set)), list(node_set)))
        edges_to_remove_barcode = [(barcode_dict[i], barcode_dict[j])
                                   for i, j in edges_to_remove]
        sub_graph_copy = nx.Graph(g.subgraph(node_set))
        sub_graph_copy.remove_edges_from(edges_to_remove_barcode)
        cos_components = list(nx.connected_components(sub_graph_copy))
        for com in cos_components:
            communities.append(com)
    return communities


def determine_molecules_bc_k_cliques(g, u):
    """
    Assign the neighbours of this vertex to molecules by
    applying k-cliques community detection for each bi-connected component.
    """
    return [community
            for bi_connected_component in
            detect_communities_biconnected_components(g, set(g.neighbors(u)))
            for community in
            detect_communities_k_clique(g, bi_connected_component)]


def determine_molecules_bc_louvain(g, u):
    """
    Assign the neighbours of this vertex to molecules by
    applying louvain for each bi-connected component.
    """
    return [community
            for bi_connected_component in
            detect_communities_biconnected_components(g, set(g.neighbors(u)))
            for community in
            detect_communities_louvain(g, bi_connected_component)]


def determine_molecules_bc_cosine_of_squared(g, u):
    """
    Assign the neighbours of this vertex to molecules by
    applying cosine of squared of the adjacency matrix for each bi-connected component.
    """
    return [community
            for bi_connected_component in
            detect_communities_biconnected_components(g, set(g.neighbors(u)))
            for community in
            detect_communities_cosine_of_squared(g, bi_connected_component)]


def determine_molecules_consensus(g, u):
    """
    Assign the neighbours of this vertex to molecules
    by Applying a stack of different approaches.
    """
    communities = []
    communities2 = []
    for bi_connected_component in detect_communities_biconnected_components(g, set(g.neighbors(u))):
        for comset in detect_communities_k_clique(g, bi_connected_component):
            communities.extend(comset)
    for community in communities:
        for comset in detect_communities_cosine_of_squared(g, community, squaring=False, threshold=0.5):
            communities2.extend(comset)
    communities = []
    for community2 in communities2:
        for comset in detect_communities_cosine_of_squared(g, community2):
            communities.extend(comset)
    return communities
    # communities2 = []
    # for community in communities:
    #     communities2.extend(detect_communities_louvain(g, community))
    # return communities2
    # return [community2
    #         for bi_connected_component in
    #         detect_communities_biconnected_components(g, set(g.neighbors(u)))
    #         for community in
    #         detect_communities_k_clique(g, bi_connected_component)
    #         for community2 in
    #         detect_communities_cosine_of_squared(g, community, squaring=False)
    #         # for community3 in
    #         # detect_communities_cosine_of_squared(g, community2)
    #         # for community4 in
    #         # detect_communities_louvain(g, community3)
    #         ]


def partition_subgraph_into_bins_randomly(node_set, max_size=40):
    """
    Partition the subgraph into bins randomly for faster processing. Return bins.
    Warning: This function is not deterministic.
    """
    bins_count = 1 + len(node_set) // max_size
    node_list = list(node_set)
    random.shuffle(node_list)
    size, leftover = divmod(len(node_set), bins_count)
    bins = [node_list[0 + size * i: size * (i + 1)] for i in range(bins_count)]
    edge = size * bins_count
    for i in range(leftover):
        bins[i % bins_count].append(node_list[edge + i])
    print(bins)
    return [set(x) for x in bins]
    # bin_sets = [set() for _ in bins]
    # for i, c in enumerate(bins):
    #     bin_sets[i].update(c)
    #     print(c)
    # return bin_sets

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
        # for chunk_community in community_detection_k_clique(g, chunk, 3):
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


def detect_and_cut_junctions_of_tree(gcomponent, messages, junction_threshold=200):
    """"
    detect the junctions in the tree, and split the tree from those points.
    """
    gcomponent = gcomponent.copy()
    set_of_all_junctions = {node
                            for node in list(gcomponent.nodes)
                            if gcomponent.degree(node) > 3}
    nodes_to_remove = []
    for candidate in set_of_all_junctions:
        candidate_messages = [messages[(candidate, e)] for e in gcomponent.neighbors(candidate)]
        candidate_messages.sort()
        if candidate_messages[-3] > junction_threshold:
            nodes_to_remove.append(candidate )
    for node_to_remove in nodes_to_remove:
        gcomponent.remove_node(node_to_remove)
    return gcomponent


def determine_safer_backbones_of_trees(g, pruning_threshold):
    """"Determine the backbones of the maximum spanning trees
            and remove branches smaller than branch_size."""
    paths = []
    for component in nx.connected_components(g):
        gcomponent = g.subgraph(component)
        messages = determine_reachability_of_tree_by_message_passing(gcomponent)
        gcomponents = detect_and_cut_junctions_of_tree(gcomponent, messages, pruning_threshold)
        for component2 in nx.connected_components(gcomponents):
            gcomponent2 = g.subgraph(component2)
            u, v, _ = diameter_of_tree(gcomponent2, weight="n")
            path = nx.shortest_path(gcomponent2, u, v, weight="n")
            paths.append(path)
    paths.sort(key=len, reverse=True)
    return paths


def determine_safer_backbones(g, pruning_threshold=100):
    """"Determine the backbones of the graph
    with ambiguous nodes being removed """
    g = g.copy()
    backbones = []
    while not nx.is_empty(g):
        gmst = g
        paths = determine_safer_backbones_of_trees(gmst, pruning_threshold)
        backbones.extend(paths)
        vertices = [u for path in paths for u in path]
        neighbors = [v for u in vertices for v in g.neighbors(u)]
        g.remove_nodes_from(vertices)
        g.remove_nodes_from(neighbors)
    backbones.sort(key=len, reverse=True)
    return backbones



paths = []
for component2 in nx.connected_components(gcomponents):
    gcomponent2 = sub_mst.subgraph(component2)
    u, v, _ = diameter_of_tree(gcomponent2, weight="n")
    path = nx.shortest_path(gcomponent2, u, v, weight="n")
    paths.append(path)
paths.sort(key=len, reverse=True)


def wrap_up_messages_and_pass(mst, messages, sender, receiver):
    """Wrap up all incoming messages to this node (sender) except the one from receiver;
    and set (pass) the message from sender to receiver."""
    neighbors_to_wrap_up = list(set(mst.neighbors(sender)) - {receiver})
    if not neighbors_to_wrap_up:
        messages[(receiver, sender)] = 1
    messages[(receiver, sender)] = 1 + sum([messages[(sender, neighbor)] for neighbor in neighbors_to_wrap_up])


def determine_reachability_by_message_passing(mst):
    """Using message passing, determine for each edge of each vertex
    the number of vertices of the tree reachable from that vertex through that edge."""
    dfs = list(nx.dfs_edges(mst))
    if not dfs:
        return dict()
    stack = [dfs[0][0]]
    messages = dict()
    # Gather
    for edge in dfs:
        while stack[-1] != edge[0]:
            wrap_up_messages_and_pass(mst, messages, stack.pop(), stack[-1])
        stack.append(edge[1])
    while len(stack) != 1:
        wrap_up_messages_and_pass(mst, messages, stack.pop(), stack[-1])
    # Distribute
    for edge in dfs:
        wrap_up_messages_and_pass(mst, messages, edge[0], edge[1])
    return messages


def detect_junctions_of_tree(gcomponent, messages, junction_threshold):
    """"
    detect the junctions in the tree, and return.
    """
    gcomponent = gcomponent.copy()
    set_of_all_junctions = {node
                            for node in list(gcomponent.nodes)
                            if gcomponent.degree(node) > 2}
    nodes_to_remove = []
    for candidate in set_of_all_junctions:
        candidate_messages = [messages[(candidate, neighbor)]
                              for neighbor in gcomponent.neighbors(candidate)]
        candidate_messages.sort()
        if candidate_messages[-3] > junction_threshold:
            nodes_to_remove.append(candidate)
    return nodes_to_remove


def determine_junctions_of_trees(g, junction_threshold):
    """"
    Determine the backbones of the maximum spanning trees
    and remove junctions over junction threshold.
    """
    junctions = []
    for component in nx.connected_components(g):
        gcomponent = g.subgraph(component)
        messages = determine_reachability_by_message_passing(gcomponent)
        new_junctions = \
            detect_junctions_of_tree(gcomponent, messages, junction_threshold)
        junctions.extend(new_junctions)
    return junctions


def determine_junctions(g, junction_threshold=100):
    """"
    Determine the backbones of the graph with ambiguous nodes being removed
    """
    gmst = nx.maximum_spanning_tree(g, weight="n")
    junctions = determine_junctions_of_trees(gmst, junction_threshold)
    return junctions

def detect_communities_biconnected_components(g, node_set):
    """Separate bi-connected components. Return components."""
    cut_vertices = set(nx.articulation_points(g.subgraph(node_set)))
    components = list(nx.connected_components(g.subgraph(node_set - cut_vertices)))
    return components

def detect_communities_k_clique(g, node_set, k=3):
    """Apply k-clique community detection. Return communities."""
    return list(nx.algorithms.community.k_clique_communities(g.subgraph(node_set), k))

def determine_molecules_bc_k_cliques(g, u):
    """
    Assign the neighbours of this vertex to molecules by
    applying k-cliques community detection for each bi-connected component.
    """
    return [community
            for bi_connected_component in
            detect_communities_biconnected_components(g, set(g.neighbors(u)))
            for community in
            detect_communities_k_clique(g, bi_connected_component)]

def determine_molecules_consensus(g, node_set):
    """
    Assign the neighbours of this vertex to molecules
    by Applying a queue of different algorithms on top of each other.
    """
    from collections import deque

    communities = [node_set]
    communities_final = []

    alg_list = deque(["bc", "k3", "bc"])
    while alg_list:
        communities_final = []
        print(communities)
        algoritm = alg_list.popleft()
        if algoritm == "bc":
            for component in communities:
                communities_final.extend(
                    detect_communities_biconnected_components(g, component))
        elif algoritm == "k3":
            for component in communities:
                communities_final.extend(
                    detect_communities_k_clique(g, component))
        elif algoritm == "cos":
            for component in communities:
                communities_final.extend(
                    detect_communities_cosine_of_squared(
                        g, component, squaring=False, threshold=0.4))
        elif algoritm == "sqCos":
            for component in communities:
                communities_final.extend(
                    detect_communities_cosine_of_squared(g, component))
        communities = communities_final
    return communities_final

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

communities =[]
for bi_connected_component in detect_communities_biconnected_components(g, g[u].keys()):
    # print(bi_connected_component)
    communities.extend(detect_communities_cosine_of_squared(g, bi_connected_component, squaring=False, threshold=0.5))
print(communities)
    # for comset in detect_communities_biconnected_components(g, bi_connected_component):
    #    communities.extend(comset)