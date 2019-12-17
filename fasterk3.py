import networkx as nx
import numpy as np
import scipy as sp
from sklearn.metrics.pairwise import cosine_similarity
g=nx.Graph()
edges=[(0,1),(0,2),(1,2),(1,3),(2,3),(3,4),(3,5),(4,5),(2,6),(5,6)]
g=nx.Graph()
edges=[(0,1),(0,2),(1,2),(1,3),(1,4),(2,3),(3,4)]
g.add_edges_from(edges)
g = g.to_undirected()
adj_array = nx.adjacency_matrix(g).toarray()
squared_adj = sp.linalg.blas.sgemm(1.0, adj_array, adj_array)
candidates3 = np.where((adj_array > 0) & (squared_adj > 0))
candidates = candidates3
connectors = []
for i, j in zip(candidates[0],candidates[1]):
	if i < j:	
		print("for i, j:", i, j)
		connectors.append([sorted([i,j, k]) for k in np.where((adj_array[i] > 0) & (adj_array[j] > 0))[0]])
perc_edges =[]
clique_graph = nx.Graph()
for connector in connectors:
	for clique in connector:
		clique = tuple(clique)
		if clique not in clique_graph.nodes():
			clique_graph.add_node(clique)
		if len(connector) > 1:
			for clique2 in connector:
				clique2 = tuple(clique2)
				if not (clique2 == clique):
					perc_edges.append((clique, clique2))
clique_graph.add_edges_from(perc_edges)





membership_dict = defaultdict(list)
    for clique in cliques:
        for node in clique:
            membership_dict[node].append(clique)

#cliques = nx.find_cliques(G)
cliques
cliques = [frozenset(c) for c in cliques if len(c) >= 3]
membership_dict = defaultdict(list)
for clique in cliques:
	for node in clique:
		membership_dict[node].append(clique)













for i, connector in enumerate(connectors):
	for j, connector2 in enumerate(connectors):
		if i < j:
			
