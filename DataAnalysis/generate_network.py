import networkx as nx
from math import log
from collections import defaultdict
from os import sys
sys.path.append("..")


def generate_network():
	TFs = defaultdict(lambda: 1)
	TFs_list = []
	Regions_list = []
	Properties_list = []

	network = nx.Graph()

	file = open("graph_data/TF_list.txt", "r")
	for line in file:
		updated_line = line.replace('\t', '_').strip('\n')
		updated_line_fixed = updated_line + "_" + str(TFs[updated_line])
		TFs[updated_line] += 1
		TFs_list.append(updated_line_fixed)

	for tf in TFs_list:
		network.add_node(tf)

	file = open("graph_data/Region_list.txt", "r")
	for line in file:
		Regions_list.append(line.strip('\n'))

	for region in Regions_list:
		network.add_node(region)

	file = open("graph_data/Property_list.txt", "r")
	for line in file:
		fixed_line = "_".join(line.split()[1:])
		Properties_list.append(fixed_line.strip('\n'))

	for property in Properties_list:
		network.add_node(property)

	file = open("graph_data/TFs_to_Properties.txt", "r")
	i = 0
	for line in file:
		edges = line.split()

		for j in range(0, len(edges)):
			if edges[j] == '1':

				network.add_edge(TFs_list[i], Properties_list[j])
		i+=1

	file = open("graph_data/Regions_to_TFs.txt", "r")
	i = 0
	for line in file:
		edges = line.split()
		for j in range(0, len(edges)):
			if edges[j] == '1':
				network.add_edge(TFs_list[j], Regions_list[i])
		i += 1

	file = open("graph_data/Regions_to_Properties.txt", "r")
	i = 0
	for line in file:
		edges = line.split()
		for j in range(0, len(edges)):
			if edges[j] == '1':
				network.add_edge(Regions_list[i], Properties_list[j])
		i += 1
	# Preprocessing
	
	# Remove TFs connected to two properties or less
	cnt = 0	
	new_TFs = []
	for TF in TFs_list:
		nbhd = network.neighbors(TF)
		x = list(set(Properties_list) & set(nbhd))
		if len(x) <= 1:
			network.remove_node(TF)
			cnt +=1
		else:
			new_TFs.append(TF)
	TFs_list = new_TFs
	print "Number of TFs Removed: " + str(cnt)
	print "Length of TF list: " + str(len(TFs_list))
	cnt = 0
	for region in Regions_list:
		nbhd = network.neighbors(region)
		x = list(set(TFs_list) & set(nbhd))
		x2 = list(set(Properties_list) & set(nbhd))
		if len(x2) <= 1 and len(x) <= 0:
			network.remove_node(region)
			cnt +=1
	print "Number of Regions Removed: " + str(cnt)

	print "Number of edges: " + str(len(network.edges()))
	print "Number of nodes: " + str(len(network.nodes()))
	return network, TFs_list, Regions_list, Properties_list

