import argparse
import itertools
import multiprocessing
import threading
import os
import random
import re
import statistics
import sys
import timeit
from collections import Counter
from itertools import repeat

import networkx as nx
import tqdm
import numpy as np
import scipy as sp
#import pandas as pd
#from scipy.linalg import get_blas_funcs

from physlr.minimerize import minimerize
from physlr.read_fasta import read_fasta
from sklearn.metrics.pairwise import cosine_similarity
from scipy.linalg import get_blas_funcs

gemm = get_blas_funcs("gemm", [X, Y])


@staticmethod
def read_graph(filenames):
    "Read a graph in either GraphViz or TSV format."
    print(int(timeit.default_timer() - t0), "Reading", *filenames, file=sys.stderr)
    read_gv = False
    g = nx.Graph()
    for filename in filenames:
        with open(filename) as fin:
            c = fin.read(1)
            if c == "s":
                g = Physlr.read_graphviz(g, filename)
                read_gv = True
            elif c == "U":
                g = Physlr.read_tsv(g, filename)
            else:
                print("Unexpected graph format", c + fin.readline(), file=sys.stderr)
                sys.exit(1)
    print(int(timeit.default_timer() - t0), "Read", *filenames, file=sys.stderr)
    if read_gv:
        print(int(timeit.default_timer() - t0), "Sorting the vertices", file=sys.stderr)
        g = Physlr.sort_vertices(g)
        print(int(timeit.default_timer() - t0), "Sorted the vertices", file=sys.stderr)
    return g

@staticmethod
def determine_molecules(g, u):
    "Assign the neighbours of this vertex to molecules."
    if "AAACACCCAACTGCTA" not in u:
        components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)).union(set([u])))))
        return u, {v: i for i, vs in enumerate(components) for v in vs if v != u}
    strategy = 3  # 1:current, 2:current modified, 3:sqCos
    if strategy == 1:  # Physlr's current version (gitMaster)
        cut_vertices = set(nx.articulation_points(g.subgraph(g.neighbors(u))))
        components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)) - cut_vertices)))
        components.sort(key=len, reverse=True)
        return u, {v: i for i, vs in enumerate(components) for v in vs}
    if strategy == 2:  # Physlr's modified current version ! (gitModif)
        cut_vertices = set(nx.articulation_points(g.subgraph(g.neighbors(u))))
        components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)) - cut_vertices)))
        components.sort(key=len, reverse=True)
        return u, {v: i for i, vs in enumerate(components) if len(vs) > 1 for v in vs}
    if strategy == 3:
        sub_graph = g.subgraph(g.neighbors(u))
        nodes_count = len(sub_graph)
        print(
            int(timeit.default_timer() - t0),
            "\nEntering strategy 3 for node ", u," with num of neighs=", nodes_count,
             file=sys.stderr)
        # edges_count = sub_graph.number_of_edges()
        if nodes_count == 0:  # or edges_count == 0:
            components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)).union(set([u])))))
            # components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)))))
            return u, {v: i for i, vs in enumerate(components) for v in vs if v != u}
        adj = nx.adjacency_matrix(sub_graph)
        # edges_to_remove = np.argwhere(
        #    cosine_similarity(np.dot(adj, adj)) < threshold)
        # adj_array = nx.adjacency_matrix(sub_graph).toarray()
        # edges_to_remove = np.argwhere(cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) < threshold)
        # barcode_dict = dict(zip(range(nodes_count), list(sub_graph.nodes)))
        # edges_to_remove_barcode = [(barcode_dict[i], barcode_dict[j]) for i, j in edges_to_remove]
        # sub_graph2 = nx.Graph(sub_graph)
        # sub_graph2.remove_edges_from(edges_to_remove2)
        # sub_graph2.remove_edges_from(edges_to_remove_barcode)
        # cos_components = list(nx.connected_components(sub_graph2))
        # adj_array = nx.adjacency_matrix(sub_graph).toarray()
        edges_to_keep = np.argwhere(
            cosine_similarity(gemm(1, adj, adj)) >= threshold)
            #cosine_similarity(np.dot(adj, adj)) >= threshold)
        # edges_to_keep = np.argwhere(
        #    cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) >= threshold)
        barcode_dict = dict(zip(range(nodes_count), list(sub_graph.nodes)))
        edges_to_keep_barcode = {(barcode_dict[i], barcode_dict[j]) for i, j in edges_to_keep}
        cos_components = list(nx.connected_components(sub_graph.edge_subgraph(edges_to_keep_barcode)))
        cos_components.sort(key=len, reverse=True)
        #removed = np.argwhere(nx.adjacency_matrix(sub_graph2).toarray() != adj_array)
        #removed = np.argwhere(nx.adjacency_matrix(sub_graph2).toarray() != adj.toarray())
        multi_node_components = [i for i in cos_components if len(i) > 1]
        singletons = [i for i in cos_components if len(i) == 1]
        # return u, {v: i for i, vs in enumerate(multi_node_components) for v in vs}
        print(
            int(timeit.default_timer() - t0),
            "\nWorking", "with the one", "with # comps equal to:", multi_node_components.__len__(),
            "\n with comps equal to:", multi_node_components,
            # "\n and edges to remove:", edges_to_remove2,
            "\n and # edges to keep:", edges_to_keep.__len__(),
            "\n and # removed singletons:", singletons.__len__(),
            # "\n and edges removed then:", removed,
            # "\n and nodes list:", list(sub_graph2.nodes),
            file=sys.stderr)
       return u, {v: i for i, vs in enumerate(cos_components) if len(vs) > 1 for v in vs}

   if strategy == 4:
       sub_graph = g.subgraph(g.neighbors(u))
       nodes_count = len(sub_graph)
       # edges_count = sub_graph.number_of_edges()
       if nodes_count == 0:  # or edges_count == 0:
          components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)).union(set([u])))))
           # components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)))))
           return u, {v: i for i, vs in enumerate(components) for v in vs if v != u}
       # [ current
       cut_vertices = set(nx.articulation_points(sub_graph))
       components = list(nx.connected_components(g.subgraph(set(g.neighbors(u)) - cut_vertices)))
       components.sort(key=len, reverse=True)
       multicomp = [i for i in components if len(i) > 1]
       # current end ]
       comul_multisubcomps = []
       for sc in multicomp:
           subcomp = sub_graph.subgraph(set(sc))
           # print(nx.adjacency_matrix(subcomp).toarray())
           adj_array = nx.adjacency_matrix(subcomp).toarray()
           new_adj = np.multiply(
               cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) >= threshold, adj_array)
           row_names = list(subcomp.nodes)
           column_names = row_names
           new_adj_df = pd.DataFrame(new_adj != adj_array, columns=column_names, index=row_names)
           edges_to_remove2 = new_adj_df[new_adj_df == True].stack().index.tolist()
           subcomp2 = nx.Graph(subcomp)
           subcomp2.remove_edges_from(edges_to_remove2)
           # subcomp2 = nx.freeze(subcomp2)
           cos_components = list(nx.connected_components(subcomp2))
           cos_components.sort(key=len, reverse=True)
           multisubcomp = [i for i in cos_components if len(i) > 1]
           comul_multisubcomps = comul_multisubcomps + multisubcomp
       return u, {v: i for i, vs in enumerate(comul_multisubcomps) for v in vs}
           # sq_cos = cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array))
           # edges_to_remove = np.argwhere(sq_cos < threshold)
           # subcomp2 = nx.Graph(subcomp)
           # subcomp2.remove_edges_from(edges_to_remove)
           # subcomp2 = nx.freeze(subcomp2)
           # cos_components = list(nx.connected_components(subcomp2))
           # cos_components.sort(key=len, reverse=True)
           # removed = np.argwhere(nx.adjacency_matrix(sub_graph2).toarray() != adj_array)
           # adj = nx.adjacency_matrix(sub_graph)
           # cos = cosine_similarity(adj.dot(adj))
           # cos = cosine_similarity(adj)
           # cos = cosine_similarity(np.matmul(adj, adj))
           # new_adj = np.multiply((cos > egRem_threshold), adj.toarray())
           # edges_to_remove = np.argwhere(new_adj != adj.toarray())
           # cos = cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array))
           # edges_to_remove = np.argwhere(cos <= egRem_threshold)
           # edges_to_remove = np.argwhere(cosine_similarity(sp.linalg.blas.sgemm(1.0, adj_array, adj_array)) <= egRem_threshold)
           # len_comps = [len(i) for i in components2]
           # if len([i for i in len_comps if i > 1]) > 1:
           # numberOfMols = numberOfMols + len([i for i in len_comps if i == 1])
           # single_node_components = [i for i in components2 if len(i) == 1]
           # if len(components2) == 1:
           #     neighbor_stats.append(stat_tuple(nodes_count, edges_count))
           # if len(components2) > 1:
           #    neighbor_stats_multicomp.append(stat_tuple(nodes_count, edges_count))

def main():
    "Run Physlr."
    #Physlr().main()
    g = read_graph("/projects/btl/aafshinfard/projects/physler-dev/data/")
    


adj_numpy  = np.array([[0,6,4,2,0,0,0,1],[6,0,5,0,0,0,0,1],[4,5,0,0,0,0,0,1],[2,0,0,0,7,3,0,1],[0,0,0,7,0,4,0,1],[0,0,0,3,4,0,0,1],[0,0,0,0,0,0,0,1],[1,1,1,1,1,1,1,0]])
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
