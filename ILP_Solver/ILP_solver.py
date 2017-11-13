from gurobipy import *
import networkx
import time as python_time
from os import sys
sys.path.append("..")
from utils import execution_time, print_edges_in_graph
from collections import defaultdict

epsilon = 0.000000001



def solve_anat_instance(model, graph, node_variables_TF, node_variables_Region, time_limit=1000000000):
	"""
	Given a tripartite-connectivity problem instance, returns the corresponding transcription factors and regions used in the optimal solution found.
	If time_limit is hit, and the optimal solution has not been found, the current optimal solution will be returned.

	:param model: a Gurobi model to be optimized
	:param graph: a directed graph with attribute 'weight' on all edges
	:param node_variables_TF: a dictionary of gurobi variables corresponding to transcription factors
	:param node_variables_Region: a dictionary of gurobi variables corresponding to regions

	:return: Two dictionaries, corresponding to transcription factors and regions, with 0 or 1 values for respective entries, indicating
			 whether the corresponding
	"""

	start_time = python_time.time()

	# Model Params
	model.params.Threads = 1
	model.params.TimeLimit = time_limit


	# SOLVE AND RECOVER SOLUTION

	print('-----------------------------------------------------------------------')
	model.optimize()

	# Recover minimal subgraph
	node_variables_TF_result, node_variables_Region_result = retreive_and_print_nodes(model, node_variables_TF, node_variables_Region)


	end_time = python_time.time()
	days, hours, minutes, seconds = execution_time(start_time, end_time)
	print('Tripartite-connectivity solving took %s days, %s hours, %s minutes, %s seconds' % (days, hours, minutes, seconds))

	# Return solution iff found

	return node_variables_TF_result, node_variables_Region_result if (model.status == GRB.status.OPTIMAL or model.status == GRB.status.TIME_LIMIT) else None




def generate_anat_model(graph, TFs_list, Regions_list, Properties_list, beta, gamma, theta, tau, max_duplicate_tfs = 2, max_percentage_penalized_edges = .2, minimum_percentage_appearance = .4, forced_nodes_full = [], special_genes = [], extra_penalized_edges = [], ):
	"""
	Generate Gurobi model to solve tripartite-connectivity problem.

	:param graph: a directed graph with attribute 'weight' on all edges
	:param TFs_list: List of all transcription factors
	:param Regions_list: List of all regions
	:param Properties_list: List of all properties
	:param beta: Number of TFs a property must be connected with
	:param gamma: Number of Regions a property must be connected with
	:param theta: Number of TFs each Region must be connected to
	:param tau: Number of Regions each TF must be connected to
	:param max_duplicate_tfs: Maximum Number of duplicates of same TF allowed in solution
	:param max_percentage_penalized_edges: Max percentage of penalized edges to allow in the solution
	:param minimum_percentage_appearance: Min percentage of forced regions
	:param forced_nodes_full: TFs to force to be used in the solution
	:param special_genes 11 genes which we wanted to see Regions for
	:param extra_penalized_edges: Edges which were not statistically significant

	:return: a Gurobi model pertaining to the Tripartite-connectivity instance, mappings from TFs to TF gurobi variables,
			 and Regions to Region gurobi variables
	"""


	# Source get +len(destination) sourceflow, destinations get -1, other nodes 0

	# Create empty optimization model
	model = Model('anat')



	# Create variables d_{uv}
	#Flow for edges
	node_variables_TF = {}
	for u in TFs_list:
		#node_variables_TF[u] = model.addVar(vtype=GRB.BINARY, name='theta_tf_%s' % (u,))
		node_variables_TF[u] = model.addVar(vtype=GRB.BINARY, name='theta_tf_%s' % (u,))
		#node_variables_TF[u] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, name='theta_tf_%s' % (u,))
	node_variables_Region = {}
	for u in Regions_list:
		#node_variables_Region[u] = model.addVar(vtype=GRB.BINARY, name='theta_region_%s' % (u,))
		node_variables_Region[u] = model.addVar(vtype=GRB.BINARY, name='theta_region_%s' % (u,))
		#node_variables_Region[u] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, name='theta_region_%s' % (u,))

	edges = [(u,v) for u in TFs_list for v in Regions_list if graph.has_edge(u,v)]
	edges_variables = {}
	for u,v in edges:
		edges_variables[u,v] = model.addVar(vtype=GRB.BINARY, name='edge_%s_%s' % (u,v))
		#edges_variables[u,v] = model.addVar(vtype=GRB.CONTINUOUS, lb=0, ub = 1, name='edge_%s_%s' % (u,v))

	model.update()

	#Constrained vs unconstrained edges(Ex: max constrained edges 20%)
	model.addConstr(quicksum(edges_variables[u,v] for u,v in edges_variables) >= max_percentage_penalized_edges * quicksum(edges_variables[u,v] for u,v in set(extra_penalized_edges) & set(edges_variables.keys())))

	# Constraint on number of touched regions
	for i in range(1, 5):
		s = set()
		for region in Regions_list:
			ss = set()
			for gene in set(graph.neighbors(region)):
				for geneus in special_genes:
					if geneus in gene:
						ss.add(geneus)
			if len(ss) == i:
				s.add(region)
		model.addConstr(quicksum(node_variables_Region[u] for u in s) >= minimum_percentage_appearance * len(s))

	# CONSTRAINTS
	unique_TFs = set()
	for u in TFs_list:
		l = u.index("_")
		unique_TFs.add(u[0:l])

        # At most TF appears twice
        for u in unique_TFs:
                model.addConstr(quicksum(node_variables_TF[x] for x in TFs_list if u in x) <= max_duplicate_tfs)

	# At least TF appears once
	for u in forced_nodes_full:
		if len([node_variables_TF[x] for x in TFs_list if u in x]) == 0:
			continue
		model.addConstr(quicksum(node_variables_TF[x] for x in TFs_list if u in x) >= 1)

	for u in TFs_list:
		model.addConstr(node_variables_TF[u] >= 0)


	for u in Regions_list:
		model.addConstr(node_variables_Region[u] >= 0)

	for u,v in edges:
		model.addConstr(edges_variables[u,v] >= node_variables_TF[u] + node_variables_Region[v] - 1)


	# TF Region 1 constraint (Ie with TFs)
	for node in TFs_list:
		Region_neighbors = list(set(graph.neighbors(node)) & set(Regions_list))

		model.addConstr(
			quicksum(node_variables_Region[u] for u in Region_neighbors) >= min(max(5,len(Region_neighbors)), theta) * node_variables_TF[node]
		)

	# TF Region 2 constraint (Ie with TFs)
	for node in Regions_list:
		TF_neighbors = list(set(graph.neighbors(node)) & set(TFs_list))

		model.addConstr(
			quicksum(node_variables_TF[u] for u in TF_neighbors) >= min(len(TF_neighbors), tau) * node_variables_Region[node]
		)

	# Beta constraint (Ie with TFs)
	for node in Properties_list:
		TF_neighbors = list(set(graph.neighbors(node)) & set(TFs_list))
		unique_TFs_beta = defaultdict(int)
		for u in TF_neighbors:
			l = u.index("_")
			unique_TFs_beta[u[0:l]] = min(max_duplicate_tfs, unique_TFs_beta[u[0:l]]+1)
		summation = 0
		for key,value in unique_TFs_beta.items():
			summation += value
		print summation, len(TF_neighbors)
		model.addConstr(
			quicksum(node_variables_TF[u] for u in TF_neighbors) >= min(summation, beta)
		)

	# Gamma constraint (Ie with Regions)
	for node in Properties_list:
		Regions_neighbors = list(set(graph.neighbors(node)) & set(Regions_list))
		model.addConstr(
			quicksum(node_variables_Region[u] for u in Regions_neighbors) >= min(len(Regions_neighbors), gamma)
		)

	# OBJECTIVE
	# Minimize total path weight
	edges = [(u,v) for u in TFs_list for v in Regions_list if graph.has_edge(u,v)]
	objective_expression = quicksum(node_variables_Region[u] for u in Regions_list) + quicksum(
		edges_variables[u,v] for u, v in edges)
	model.setObjective(objective_expression, GRB.MINIMIZE)

	return model, node_variables_TF, node_variables_Region



def retreive_and_print_nodes(model, node_variables_TF, node_variables_Region):
	"""

	:param model: an optimized gurobi model
	:param node_variables_TF: a dictionary mapping original TFs to corresponding gurobi variables
	:param node_variables_Region: a dictionary mapping original Regions to corresponding gurobi variables
	:param detailed_output: flag which when True will print the edges in the optimal subgraph

	:return: Two dictionaries, corresponding to transcription factors and regions, with 0 or 1 values for respective entries, indicating
			 whether the corresponding
	"""

	node_variables_TF_result, node_variables_Region_result = {},{}
	if model.status == GRB.status.OPTIMAL or model.status == GRB.status.TIME_LIMIT:
		node_variables_TF_result = model.getAttr('x', node_variables_TF)
		node_variables_Region_result = model.getAttr('x', node_variables_Region)
	return node_variables_TF_result, node_variables_Region_result

