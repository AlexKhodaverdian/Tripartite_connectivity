import numpy as np
import matplotlib.pyplot as plt

from generate_network import generate_network
from ILP_Solver.ILP_solver import *

alpha = 4
gamma = 9
network, TFs_list, Regions_list, Properties_list = generate_network()

deg_TF_Properties = []
deg_TF_Regions = []

for node in TFs_list:
	deg_TF_Properties.append(len(list(set(network.neighbors(node)) & set(Properties_list))))
	deg_TF_Regions.append(len(list(set(network.neighbors(node)) & set(Regions_list))))

plt.scatter(deg_TF_Properties, deg_TF_Regions)
plt.title("Degree Pair Plot for TFs")
plt.xlabel("Number of Property Neighbors")
plt.ylabel("Number of Region Neighbors")
plt.savefig("Results/Plots/TFs_PairPlot_All_Data.png")
plt.clf()

deg_Properties_TF = []
deg_Properties_Regions = []

for node in Properties_list:
	deg_Properties_TF.append(len(list(set(network.neighbors(node)) & set(TFs_list))))
	deg_Properties_Regions.append(len(list(set(network.neighbors(node)) & set(Regions_list))))

plt.scatter(deg_Properties_TF, deg_Properties_Regions)
plt.title("Degree Pair Plot for Properties")
plt.xlabel("Number of TF Neighbors")
plt.ylabel("Number of Region Neighbors")
plt.savefig("Results/Plots/Properties_PairPlot_All_Data.png")
plt.clf()

deg_Region_TF = []
deg_Region_Property = []

for node in Regions_list:
	deg_Region_TF.append(len(list(set(network.neighbors(node)) & set(TFs_list))))
	deg_Region_Property.append(len(list(set(network.neighbors(node)) & set(Properties_list))))

plt.scatter(deg_Region_TF, deg_Region_Property)
plt.title("Degree Pair Plot for Regions")
plt.xlabel("Number of TF Neighbors")
plt.ylabel("Number of Property Neighbors")
plt.savefig("Results/Plots/Region_PairPlot_All_Data.png")
plt.clf()

used_TFs = []
used_Regions = []
for line in open("Results/TFs_used.txt"):
	used_TFs.append(line.strip())
for line in open("Results/Regions_used.txt"):
	used_Regions.append(line.strip())

deg_TF_Properties2 = []
deg_TF_Regions2 = []

for node in used_TFs:
	deg_TF_Properties2.append(len(list(set(network.neighbors(node)) & set(Properties_list))))
	deg_TF_Regions2.append(len(list(set(network.neighbors(node)) & set(used_Regions))))

plt.scatter(deg_TF_Properties2, deg_TF_Regions2)
plt.title("Degree Pair Plot for TFs Used In Solution")
plt.xlabel("Number of Property Neighbors")
plt.ylabel("Number of Region Neighbors Used In Solution")
plt.savefig("Results/Plots/TFs_PairPlot_Used_TFs_Vs_Used_Regions.png")
plt.clf()

deg_Region_TF2 = []
deg_Region_Property2 = []

for node in used_Regions:
	deg_Region_TF2.append(len(list(set(network.neighbors(node)) & set(used_TFs))))
	deg_Region_Property2.append(len(list(set(network.neighbors(node)) & set(Properties_list))))

plt.scatter(deg_Region_TF2, deg_Region_Property2)
plt.title("Degree Pair Plot for Regions Used In Solution")
plt.xlabel("Number of TF Neighbors Used In Solution")
plt.ylabel("Number of Property Neighbors")
plt.savefig("Results/Plots/Region_PairPlot_Used_Regions_Vs_Used_TFs.png")
plt.clf()


deg_TF_Properties3 = []
deg_TF_Regions3 = []

for node in used_TFs:
	deg_TF_Properties3.append(len(list(set(network.neighbors(node)) & set(Properties_list))))
	deg_TF_Regions3.append(len(list(set(network.neighbors(node)) & set(Regions_list))))

plt.scatter(deg_TF_Properties3, deg_TF_Regions3)
plt.title("Degree Pair Plot for TFs Used In Solution")
plt.xlabel("Number of Property Neighbors")
plt.ylabel("Number of Region Neighbors")
plt.savefig("Results/Plots/TFs_PairPlot_Used_TFs_Vs_All_Regions.png")
plt.clf()

deg_Region_TF3 = []
deg_Region_Property3 = []

for node in used_Regions:
	deg_Region_TF3.append(len(list(set(network.neighbors(node)) & set(TFs_list))))
	deg_Region_Property3.append(len(list(set(network.neighbors(node)) & set(Properties_list))))

plt.scatter(deg_Region_TF3, deg_Region_Property3)
plt.title("Degree Pair Plot for Regions Used In Solution")
plt.xlabel("Number of TF Neighbors")
plt.ylabel("Number of Property Neighbors")
plt.savefig("Results/Plots/Region_PairPlot_Used_Regions_Vs_All_TFs.png")
plt.clf()


plt.scatter(deg_TF_Properties, deg_TF_Regions, c='b', label="All TFs")
plt.scatter(deg_TF_Properties3, deg_TF_Regions3, c='r', label="TFs Used In Sol")
plt.title("Degree Pair Plot for TFs")
plt.xlabel("Number of Property Neighbors")
plt.ylabel("Number of Region Neighbors")
plt.legend()
plt.savefig("Results/Plots/TFs_PairPlot_Mixed.png")
plt.clf()


plt.scatter(deg_Region_TF, deg_Region_Property, c='b', label="All Regions")
plt.scatter(deg_Region_TF3, deg_Region_Property3,c='r', label="Regions Used In Sol")
plt.title("Degree Pair Plot for Regions Used")
plt.xlabel("Number of TF Neighbors")
plt.ylabel("Number of Property Neighbors")
plt.savefig("Results/Plots/Region_PairPlot_Mixed.png")
plt.legend()
plt.clf()
