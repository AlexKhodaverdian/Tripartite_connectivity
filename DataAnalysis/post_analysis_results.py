from generate_network import generate_network
import csv
#from os import sys
#sys.path.append("..")
from ILP_Solver.ILP_solver import *
import numpy as np

network, TFs_list, Regions_list, Properties_list = generate_network()
TFs_used= []
Regions_used = []
for TF in open("TFs_used.txt", "r"):
	TFs_used.append(TF.strip())
print TFs_used
for region in open("Regions_used.txt", "r"):
	Regions_used.append(region.strip())
print Regions_used
TF_TF_column = ["Transcription Factors"]
TF_Region_column = ["Regions"]
TF_Property_column = ["Properties"]
total = 0
for node in TFs_used:
	Region_neighbors_for_TFs = list(set(Regions_used) & set(network.neighbors(node)))
	total += len(Region_neighbors_for_TFs)
	Property_neighbors_for_TFs = list(set(Properties_list) & set(network.neighbors(node)))
	TF_TF_column.append(node)
	TF_Region_column += Region_neighbors_for_TFs
	TF_Property_column += Property_neighbors_for_TFs

	TF_TF_column += [node for _ in range(0, max(len(Region_neighbors_for_TFs), len(Property_neighbors_for_TFs)) - 1)]
	TF_Region_column += ["" for _ in range(0, max(len(Region_neighbors_for_TFs), len(Property_neighbors_for_TFs)) - len(Region_neighbors_for_TFs))]
	TF_Property_column += ["" for _ in range(0, max(len(Region_neighbors_for_TFs), len(Property_neighbors_for_TFs)) - len(Property_neighbors_for_TFs))]
	print len(TF_TF_column), len(TF_Region_column),len( TF_Property_column)
	

Region_TF_column = ["Transcription Factors"]
Region_Region_column = ["Regions"]
Region_Property_column = ["Properties"]
for node in Regions_used:
	TF_neighbors_for_Regions = list(set(TFs_used) & set(network.neighbors(node)))
	Property_neighbors_for_Regions = list(set(Properties_list) & set(network.neighbors(node)))
	Region_Region_column.append(node)
	Region_TF_column += TF_neighbors_for_Regions
	Region_Property_column += Property_neighbors_for_Regions

	Region_Region_column += [node for _ in range(0, max(len(TF_neighbors_for_Regions), len(Property_neighbors_for_Regions)) - 1)]
	Region_TF_column += ["" for _ in range(0, max(len(TF_neighbors_for_Regions), len(Property_neighbors_for_Regions)) - len(TF_neighbors_for_Regions))]
	Region_Property_column += ["" for _ in range(0, max(len(TF_neighbors_for_Regions), len(Property_neighbors_for_Regions)) - len(Property_neighbors_for_Regions))]
	print len(Region_TF_column), len(Region_Region_column),len( Region_Property_column)


Property_TF_column = ["Transcription Factors"]
Property_Region_column = ["Regions"]
Property_Property_column = ["Properties"]
for node in Properties_list:
	TF_neighbors_for_Property = list(set(TFs_used) & set(network.neighbors(node)))
	Region_neighbors_for_Property = list(set(Regions_used) & set(network.neighbors(node)))
	Property_Property_column.append(node)
	Property_TF_column += TF_neighbors_for_Property
	Property_Region_column += Region_neighbors_for_Property

	Property_Property_column += [node for _ in range(0, max(len(TF_neighbors_for_Property), len(Region_neighbors_for_Property)) - 1)]
	Property_TF_column += ["" for _ in range(0, max(len(TF_neighbors_for_Property), len(Region_neighbors_for_Property)) - len(TF_neighbors_for_Property))]
	Property_Region_column += ["" for _ in range(0, max(len(TF_neighbors_for_Property), len(Region_neighbors_for_Property)) - len(Region_neighbors_for_Property))]
	print len(Property_TF_column), len(Property_Region_column),len( Property_Property_column)

TFs = zip(TF_TF_column, TF_Region_column, TF_Property_column)
Regions = zip(Region_Region_column, Region_TF_column, Region_Property_column)
Properties = zip(Property_Property_column, Property_TF_column, Property_Region_column)

avg_deg_all = []
for TF in TFs_list:
	avg_deg_all.append(len(set(network.neighbors(TF)) & set(Regions_list)))

avg_deg_TFs_chosen = []
for TF in TFs_used:
	avg_deg_TFs_chosen.append(len(set(network.neighbors(TF)) & set(Regions_list)))


with open('TF_detailed_results.csv', 'w') as fp:
	a = csv.writer(fp, delimiter=',')	
	a.writerows(TFs)

with open('Regions_detailed_results.csv', 'w') as fp:
	a = csv.writer(fp, delimiter=',')	
	a.writerows(Regions)

with open('Properties_detailed_results.csv', 'w') as fp:
	a = csv.writer(fp, delimiter=',')	
	a.writerows(Properties)

with open("Meta_stats.txt", "w") as fp:
	fp.write("Total Cost: " + str(total + len(Regions_used)))
	fp.write("\n")

	fp.write("Number of Edges Between TFs and Regions: " + str(total))
	fp.write("\n")

	fp.write("Number of Regions Used: " + str(len(Regions_used)))
	fp.write("\n")

	fp.write("Number of TFs Used: " + str(len(TFs_used)))
	fp.write("\n")
	fp.write("\n")


	fp.write("============= SUMMARY STATISTICS ================")
	fp.write("\n")

	fp.write("Mean TF Degree (wrt Regions) For All TFs: " + str(np.mean(avg_deg_all)))
	fp.write('\n')

	fp.write("Mean TF Degree (wrt Regions) For TFs Chosen: " + str(np.mean(avg_deg_TFs_chosen)))
	fp.write('\n')

	fp.write("Median TF Degree (wrt Regions) For All TFs: " + str(np.median(avg_deg_all)))
	fp.write('\n')

	fp.write("Median TF Degree (wrt Regions) For TFs Chosen: " + str(np.median(avg_deg_TFs_chosen)))
	fp.write('\n')
