##Ford–Fulkerson algorithm applied to find the autocatalytic cycles
##
import sys
import os
import glob
import shutil
import os.path
from pathlib import Path
sys.setrecursionlimit(20**6) 
from collections import deque
import operator


def only_spontaneous_reaction(f,b1,b2):
	t = f + b1 + b2
	#print(t)
	total_spontaneous = True
	total_energy = 0
	reaction_list =[]
	energy_list = []
	spontaneidad = True
	for st in t:
		if st[0:-2].isnumeric():
			reaction_list.append(st)
			#if spontan[st[:-2]] == "NoValue":
			if spontan[st] == "NoValue":
				energy_list.append("NoValue")
			else:
				#energy_list.append(float(spontan[st[0:-2]]))
				#cheking = float(spontan[st[0:-2]])<cout_energy
				energy_list.append(float(spontan[st]))
				cheking = float(spontan[st])<cout_energy				
				spontaneidad = spontaneidad  and cheking
	if "NoValue" in energy_list:
		total_energy = "NoValue"
		total_spontaneous = "NoValue"
	else:
		total_energy = sum(energy_list)
		if total_energy<cout_energy and spontaneidad:
			total_spontaneous = "Yes"
		else:
			total_spontaneous = "Not"
	return total_spontaneous,total_energy,reaction_list,energy_list

class Edge:
    def __init__(self, u, v, r):
        self.u = u
        self.v = v
        self.r = r

    def __str__(self) -> str:
        return "(" + self.u + "," + self.v + ")"

    __repr__ = __str__


def make_edges(data):
	edges = []
	for u, v, r, id, energy in data:
		edges.append(Edge(u, v, r))
	return edges			

def make_graph(edges):
    graph = dict()
    for e in edges:
        if e.u in graph:
            graph[e.u].append(e.v)
        else:
            graph[e.u] = [e.v]
    return graph


def get_ids(G):
    n = 0
    p = dict()
    for v in G.keys():
        if not v in p:
            p[v] = n
            n += 1
        for u in G[v]:
            if not u in p:
                p[u] = n
                n += 1
    return p

def get_canonical(G, p):
    adj = dict()
    for i in range(len(p)):
        adj[i] = []
    for v in G.keys():
        for u in G[v]:
            adj[p[v]].append(p[u])
    return adj




def augmenting_path(s, t):
    Q = deque()
    Q.append(s)
    parent = dict()
    parent[s] = -1
    while len(Q) > 0:
        q = Q.popleft()
        if q == t:
            return True, parent
        for u in adj[q]:
            if flo[(q, u)] >= cap[(q, u)]:
                continue
            if u not in parent:
                parent[u] = q
                Q.append(u)
    return False, dict()

def get_flow(s, t):
    for key in flo.keys():
        flo[key] = 0
    flow = 0
    while flow < 2:
        has_flow, parent = augmenting_path(s, t)
        if has_flow:
            flow += 1
        else:
            break
        v = 1*t
        while parent[v] >= 0:
            flo[(parent[v], v)] += 1
            flo[(v, parent[v])] -= 1
            v = parent[v]

    edges = []
    for key, f in flo.items():
        if f > 0:
            for i in range(f):
                edges.append(key)

    if flow == 2:
        leng = max(len(disjoint_paths(s, t, edges)[0]),len(disjoint_paths(s, t, edges)[1]))
        return True, edges,leng
    else:
        return False, [],0

def get_path(x, parent):
    path = []
    while x != -1:
        path.append(x)
        x = parent[x]
    path.reverse()
    return path

def disjoint_paths(s, t, edges):
    rho = dict()
    for u, v in edges:
        rho[u] = []
        rho[v] = []

    for u, v in edges:
        rho[u].append(v)

    w = rho[s][0]
    path1 = [s, w]
    while w != t:
        w = rho[w][0]
        path1.append(w)

    w = rho[s][1]
    path2 = [s, w]
    while w != t:
        w = rho[w][0]
        path2.append(w)
    return [path1, path2]

def translate(i):
	return prev[i] 
	
def bfs(src, na):
    print("### Star")
    Q = deque()
    Q.append(src)
    parent = dict()
    parent[src] = -1
    name_cicles = na
    while len(Q) > 0:
    	q = Q.popleft()
    	print("### Sub ###")
    	co = 0
    	for u in adj[q]:
    		co +=1
    		#print(len(adj[q])) #Print the expected iteration 
    		#print(co) #Print the actual iteration
    		if u not in parent:
    			parent[u] = q
    			Q.append(u)
    			if prev[u][:-2].isnumeric():
    				two_paths, edges,le_1 = get_flow(u, src)
    				if two_paths:
    					smile_list = list(map(translate, get_path(u, parent)))
    					path_two_a = []
    					for path_ad in disjoint_paths(u, src, edges):
	    						path_two_a.append(list(map(translate, path_ad))[1:-1])
	    				print(only_spontaneous_reaction(smile_list,path_two_a[0],path_two_a[1]))
	    				spo,energia,reac_li,ener_li = only_spontaneous_reaction(smile_list,path_two_a[0],path_two_a[1])
	    				if spo=="Yes":
	    					with open(dir_results+"/"+"AutoCatCycles_"+name_cicles+".txt", 'a') as type_file:
	    						smile_list = list(map(translate, get_path(u, parent)))
	    						type_file.write(smile_list[0]+"\t")
	    						type_file.write(smile_list[-1]+"\t")
	    						type_file.write(",".join(smile_list)+"\t")
	    						type_file.write(str(len(smile_list)))
	    						print(list(map(translate, get_path(u, parent))))
	    						path_two = []
	    						for path_ad in disjoint_paths(u, src, edges):
	    							type_file.write("\t")
	    							smile_list2 = list(map(translate, path_ad))[1:-1]
	    							type_file.write(",".join(smile_list2))
	    							type_file.write("\t")
	    							type_file.write(str(len(smile_list2)))
	    							path_two.append(list(map(translate, path_ad)))
	    						print(path_two)
	    						type_file.write("\t"+str(energia)+"\t"+str(reac_li)+"\t"+str(ener_li)+"\t"+"Strictly_Sponstaneous")
	    						type_file.write("\n")

	    				if spo == "Not":
	    					with open(dir_results+"/"+"AutoCatCycles_"+name_cicles+".txt", 'a') as type_file:
	    						smile_list = list(map(translate, get_path(u, parent)))
	    						type_file.write(smile_list[0]+"\t")
	    						type_file.write(smile_list[-1]+"\t")
	    						type_file.write(",".join(smile_list)+"\t")
	    						type_file.write(str(len(smile_list)))
	    						print(list(map(translate, get_path(u, parent))))
	    						path_two = []
	    						for path_ad in disjoint_paths(u, src, edges):
	    							type_file.write("\t")
	    							smile_list2 = list(map(translate, path_ad))[1:-1]
	    							type_file.write(",".join(smile_list2))
	    							type_file.write("\t")
	    							type_file.write(str(len(smile_list2)))
	    							path_two.append(list(map(translate, path_ad)))
	    						print(path_two)
	    						type_file.write("\t"+str(energia)+"\t"+str(reac_li)+"\t"+str(ener_li)+"\t"+"Non_Strictly_Sponstaneous")
	    						type_file.write("\n")	    					

	    				if spo == "NoValue":
	    					with open(dir_results+"/"+"AutoCatCycles_"+name_cicles+".txt", 'a') as type_file:
	    						smile_list = list(map(translate, get_path(u, parent)))
	    						type_file.write(smile_list[0]+"\t")
	    						type_file.write(smile_list[-1]+"\t")
	    						type_file.write(",".join(smile_list)+"\t")
	    						type_file.write(str(len(smile_list)))
	    						print(list(map(translate, get_path(u, parent))))
	    						path_two = []
	    						for path_ad in disjoint_paths(u, src, edges):
	    							type_file.write("\t")
	    							smile_list2 = list(map(translate, path_ad))[1:-1]
	    							type_file.write(",".join(smile_list2))
	    							type_file.write("\t")
	    							type_file.write(str(len(smile_list2)))
	    							path_two.append(list(map(translate, path_ad)))
	    						print(path_two)
	    						type_file.write("\t"+str(energia)+"\t"+str(reac_li)+"\t"+str(ener_li)+"\t"+"NonValue")
	    						type_file.write("\n")






def search_autocatalysis_cycles (dir_cycle_results,spontaneous,data, energy_cut_off, erase_water,to_seach):
	
	global dir_results
	global cout_energy 
	global G
	global p
	global adj
	global prev
	global cap 
	global flo
	global spontan
	cout_energy = energy_cut_off
	spontan = spontaneous.copy()
	to_search = to_seach
	dir_results = dir_cycle_results[:-1]
	#filed = pathlib.Path(dir_results)
	filed = Path(dir_results)
	if filed.exists ():
		question = input("Do you want to erase the previos " + dir_results + " directory (yes/not):")
		if question == "yes" or question == "Yes" or question == "YES":
			shutil.rmtree(dir_results)
			os.mkdir(dir_results)
		else:
			print("Finishing the script")
			exit()	
	else:
		os.mkdir(dir_results)
		cout_energy = energy_cut_off

	G = make_graph(make_edges(data))
	p = get_ids(G)
	adj = get_canonical(G, p)
	n = len(p)
	
	if erase_water:
		## quitanto el O de los ciclos
		adj[p["O"]] = []
	prev = list(p.keys())

	cap = dict()
	flo = dict()
	for i in range(n):
		for u in adj[i]:
			cap[(i, u)] = 0
			cap[(u, i)] = 0
			flo[(i, u)] = 0
			flo[(u, i)] = 0

	for i in range(n):
		for u in adj[i]:
			cap[(i, u)] += 1
	
	if to_seach != "all":
		for a,b in p.items():
			if a == to_search:  ## Chose catalitic node
				mole_check = b
				print(b)
				nam = a
		print(prev[mole_check])
		if not prev[mole_check][:-2].isnumeric() and prev[mole_check] != "O": 
			bfs(mole_check,nam)
	
	else:
		for a,b in p.items():
			if not prev[b][:-2].isnumeric() and prev[b] != "O":
				bfs(b,a)


def get_feeder_incycle_count_dir(mol,plots_dir,dir_cycle_results,reac_dict,prod_dict):
	
	#reac_dict = reac_dict_5.copy()
	#prod_dict = prod_dict_5.copy()
	results_dir = dir_cycle_results[:-1]
	#print(results_dir )
	#print("Here")
	list_of_files = glob.glob(results_dir+"/AutoCatCycles*")
	#print(list_of_files)
	filed = Path(plots_dir)
	if filed.exists ():
		question = input("Do you want to erase the previos " + plots_dir + " directory (yes/not):")
		if question == "yes" or question == "Yes" or question == "YES":
			shutil.rmtree(plots_dir)
			os.mkdir(plots_dir)
	else:
		os.mkdir(plots_dir)


	as_feeder_global = {}
	as_consumer_global = {}
	for filesearch in list_of_files:
		print(filesearch)
	for file_name_out in list_of_files:
		len_root = len(results_dir) + len("/AutoCatCycles_")
		molecule =file_name_out[len_root:-4] 
		#print(file_name_out)
		#print(len_root)
		if mol != "all":
			if  molecule != mol:
				continue

		file_name = results_dir + "/AutoCatCycles_" + molecule + ".txt"

		filed = Path(file_name)
		if not filed.is_file():
			print("There isn't the file: ",file_name)
			exit()	
			
		cycles = open(file_name,'r')  ###### Name of the output of cycles you want to analyst
		lines = cycles.readlines()
		print("Here")
		count = 0
		
		for line in lines:
			count += 1
			#print(count)
			###########
			line=line.rstrip('\n')
			line=line.split('\t')
			cat_node = line[0]
			cat_reaction = line[1]
			branch_0 = line[2]
			len_0 = int(line[3])
			branch_1 = line[4]
			len_1 =int(line[5])
			branch_2 = line[6]
			len_2 = int(line[7])
			total_energy = line[8]
			reaction_list = line [9]
			reaction_energy_list = line [10]
			spontaneaty = line[11]
			if spontaneaty =="Non_Strictly_Sponstaneous" or spontaneaty == "NonValue":
				continue
			print(spontaneaty)
			import draw_cycles_14 as draw_cy
			import reading_mod_outs as re_mod
			
			############
			if branch_1 !="":
				b_p = branch_0+","+branch_1
			else:
				b_p = branch_0

			branch_principal = b_p.split(",")
			if branch_2 != "":
				b_s = cat_reaction+","+branch_2+","+cat_node
			else:
				b_s = cat_reaction+","+cat_node
			branch_secundary = b_s.split(",")
			###########
			principal_len_cycle = len_0+len_1
			secundary_len_cycle = len_2+2
			###########		
			for i in range(principal_len_cycle):		
				if branch_principal[i][:-2].isnumeric():
					reac_out, prod_out,water_count_reac,water_count_prod = draw_cy.get_mol_out_cycle(True,i,branch_principal[i],branch_0,branch_1,branch_2,reac_dict,prod_dict)
					for mol_feed in reac_out:
						if mol_feed in as_feeder_global:
							as_feeder_global[mol_feed] += 1
						else:
							as_feeder_global[mol_feed] = 1
					for mol_consumer in prod_out:
						if mol_consumer in as_consumer_global:
							as_consumer_global[mol_consumer] += 1
						else:
							as_consumer_global[mol_consumer] = 1
			for j in range(1,secundary_len_cycle):		
				if branch_secundary[j][:-2].isnumeric():
					reac_out, prod_out,water_count_reac,water_count_prod = draw_cy.get_mol_out_cycle(False,j,branch_secundary[j],branch_0,branch_1,branch_2,reac_dict,prod_dict)
					for mol_feed in reac_out:
						if mol_feed in as_feeder_global:
							as_feeder_global[mol_feed] += 1
						else:
							as_feeder_global[mol_feed] = 1
					for mol_consumer in prod_out:
						if mol_consumer in as_consumer_global:
							as_consumer_global[mol_consumer] += 1
						else:
							as_consumer_global[mol_consumer] = 1
	
	print("sorting")
	as_feeder_global  = dict( sorted(as_feeder_global.items(), key=operator.itemgetter(1),reverse=True))
	as_consumer_global  = dict( sorted(as_consumer_global.items(), key=operator.itemgetter(1),reverse=True))
	
	#as_feeder_global = dict(list(as_feeder_global)[0:10])
	#as_consumer_globa = dict(list(as_consumer_global)[0:10])
	print("feeders")
	for molecule, value in list(as_feeder_global.items())[0:10]:
		print(molecule, " : ", value)
		
	print("consumers")
	for molecule, value in list(as_consumer_global.items())[0:10]:
		print(molecule, " : ", value)	
	return as_feeder_global,as_consumer_global
	

if __name__ == "__main__":
	import reading_mod_outs as re_mod
	import plots_functions_2 as plots
	plots_dir = "Plots_dir/"
	dir_cycle_results ="AutoCycles-Files/"
	data_dir ="Data/"
	file_name = "rels_5.txt"
	subplot_dir = "Feeders_dir/"
	length = 10
	figure_caption_1 = "Feeders"
	figure_caption_2 = "Consumer"
	reac_dict,prod_dict = re_mod.get_reac_prod_dirs(data_dir,file_name)	
	global_feed_dir, global_consumer_dir = get_feeder_incycle_count_dir("all",plots_dir,dir_cycle_results,reac_dict,prod_dict)
	plots.plot_feeder_incycle(global_feed_dir,plots_dir,subplot_dir,length,figure_caption_1)
	plots.plot_feeder_incycle(global_consumer_dir,plots_dir,subplot_dir,length,figure_caption_2)
	
	#print(global_feed_dir)
	#print(global_consumer_dir)
	

