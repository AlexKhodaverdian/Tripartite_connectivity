from generate_network import generate_network
#from os import sys
#sys.path.append("..")
from ILP_Solver.ILP_solver import *



network, TFs_list, Regions_list, Properties_list = generate_network()


beta = 5
gamma = 10
theta = 20
tau = 3
alpha = .5
max_duplicate_tfs = 2
forced_nodes_full = ['POU3F1',
'POU3F2',
'SOX2',
'SOX1',
'PAX6',
'OTX2',
'LHX2',
'NEUROG1',
'NEUROG2',
'NEUROD2',
'SP8',
'IRX3',
'SOX10',
'PKNOX2',
'HHEX',
'LMX1A',
'BARHL1',
'LHX5',
'NR2F2',
'DMBX1',
'MEIS2',
'OTX1',
'SOX21',
'FOXB1',
'SOX5',
'MEIS3',
'HOMEZ',
'TCF3',
'TCF4',
'ZIC1',
'ZIC2',
'ZIC3',
'ZIC4',
'ZIC5']

model, node_variables_TF, node_variables_Region = generate_anat_model(network, TFs_list, Regions_list, Properties_list, alpha, beta, gamma, theta, tau, max_duplicate_tfs, forced_nodes_full)
node_variables_TF_result, node_variables_Region_result = solve_anat_instance(model, network, node_variables_TF, node_variables_Region,time_limit=300)

file = open("Results/TFs_used.txt","w")

for node, value in node_variables_TF_result.items():
	if value == 1:
		file.write(node)
		file.write('\n')
file.close()

file = open("Results/Regions_used.txt","w")

for node, value in node_variables_Region_result.items():
	if value == 1:
		file.write(node)
		file.write('\n')
file.close()

