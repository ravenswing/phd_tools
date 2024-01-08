import sys
import string
import re
import math
import csv
import igraph
import random
import cairo
#import matplotlib as plt


def readFile(file):
	file=open(file)
	line=file.readline()
	
	text_file=""
	
	while(line):
		text_file+=line
		line=file.readline()
		
	file.close()
	
	return(text_file)
	

def matriuToString(matriu):
	espai=" "
	text=""
	
	for fila in matriu:
		text=text+espai.join(fila)
		text=text+"\n"
	
	return(text)

def tabularText(text):
	
	noms=[]
	valors=[]
	
	files=string.split(text,"\n")
	
	i=0
	while(i<len(files)):
		fila=files[i]
		if (fila != ""):
			nums=re.split("[\s]+",fila.strip())
			valors.append(nums)
		
		i=i+1
	
	return(valors)

def writeFile(filename,text):
	file=open(filename,"w+")
	
	file.write(text)

	file.close()

def random_color():
	#	a = float(random.randrange(0,255))
	#	b = float(random.randrange(0,255))
	#	c = float(random.randrange(0,255))


	#	a2=a/255
	#	b2=b/255
	#	c2=c/255

	a2=int(int(random.randrange(0,255)))
	b2=int(int(random.randrange(0,255)))
	c2=int(int(random.randrange(0,255)))

	color=[a2,b2,c2]
	
	return(color)

def get_color(index):
	colors=[]

	mediumseagreen=[60,179,113]
	dark_orange=[255,140,0]
	slate_blue=[106,90,205]
	crimson=[220,20,60]
	medium_orchid=[186,85,211]
	gold=[255,215,0]
	medium_turquoise=[72,209,204]
	deep_sky_blue=[0,191,255]
	teal=[0,128,128]
	deep_pink=[255,20,147]
	midnight_blue=[25,25,112]
	dark_sea_green=[143,188,143]
	orange_red=[255,69,0]
	green_yellow=[173,255,47]
	dark_blue=[0,0,139]
	olive=[128,128,0]
	plum=[221,160,221]
	slate_gray=[112,128,144]
	salmon=[250,128,114]
	purple=[128,0,128]
	black=[0,0,0]
	dark_red=[139,0,0]
	dark_violet=[148,0,211]
	khaki=[240,230,140]
	hot_pink=[255,105,180]
	aqua=[0,255,255]
	silver=[192,192,192]
	limegreen=[50,205,50]
	magenta=[255,0,255]
	lightslategray=[119,136,153]
	khaki=[240,230,140]
		
	colors.append(mediumseagreen)
	colors.append(dark_orange)
	colors.append(slate_blue)
	colors.append(crimson)
	colors.append(medium_orchid)
	colors.append(gold)
	colors.append(medium_turquoise)
	colors.append(plum)
	colors.append(teal)
	colors.append(deep_pink)
	colors.append(midnight_blue)
	colors.append(green_yellow)
	colors.append(deep_sky_blue)
	colors.append(purple)
	colors.append(dark_blue)
	colors.append(olive)
	colors.append(orange_red)
	colors.append(slate_gray)
	colors.append(salmon)
	colors.append(dark_sea_green)
	colors.append(black)
	colors.append(dark_red)
	colors.append(dark_violet)
	colors.append(khaki)
	colors.append(hot_pink)
	colors.append(aqua)
	colors.append(silver)
	colors.append(limegreen)
	colors.append(magenta)
	colors.append(lightslategray)
	colors.append(khaki)

	if (index<0 or index>len(colors)-1):
		return random_color()
	else:
		return colors[index]






def multiplyListValues(list,multiplicator):
	new_list=[]
	for value in list:
		new_list.append(value*multiplicator)

	return(new_list)

def findMaxMinValue(list):
	maximum_value=0
	minimum_value=0
	for value in list:
		if value > maximum_value:
			maximum_value=value
		if value < minimum_value:
			minimum_value=value
	return(maximum_value,minimum_value)


def meanVarStd(list):
	sum=0
	for value in list:
		sum=sum+value
#	print len(list)
	mean=sum/len(list)

	sum_var=0
	for value in list:
		diff=(value-mean)*(value-mean)
		sum_var=sum_var+diff
	variance=sum_var/len(list)

	std=math.sqrt(variance)

	return(mean,variance,std)


def usage():
	print "----------------------------------------------------------------------------------------------"
	print "                                   Network Analysis Usage:"
	print "----------------------------------------------------------------------------------------------"
	print "-md  : Specify distance matrix file"
	print "-mc : Specify correlation matrix file"
	print "-name: Specify the name for the image file containing the graph"
	print "-comm: Specify the algorithm for finding community structures: fastgreedy, edge_betweness, walktrap"
	print "----------------------------------------------------------------------------------------------"
	print "                                                                                              "
	print ""

def inputParser(arg):
	num_arg=len(arg)-1
	
	if(num_arg == 0) or (num_arg < 8):
		return(False,"")

	index=1
	while index <= num_arg:
		param=arg[index]
		
		if param=="-md":
			if num_arg>index:
				fitxer_distancies=str(arg[index+1])
				index=index+2
			else:
				print "ERROR. You should specify the distance matrix file"
				return(False,"")
		
		
		elif param=="-mc":
			if num_arg>index:
				fitxer_corr=str(arg[index+1])
				index=index+2
			else:
				print "ERROR. You should specify the correlation matrix file"
				return(False,"")

		elif param=="-name":
			if num_arg>index:
				figure_name=str(arg[index+1])
				index=index+2
			else:
				print "ERROR. You should specify the name for the image file"
				return(False,"")
					
		elif param=="-comm":
			if num_arg>index:
				community_method=str(arg[index+1])
				index=index+2
			else:
				print "ERROR. You should specify the algorithm for finding the community structure"
				return(False,"")
		


	return(True,fitxer_distancies,fitxer_corr,figure_name,community_method)


params=inputParser(sys.argv)

if(not params[0]):
	usage()
	exit()

fitxer_distancies = params[1]
fitxer_corr = params[2]
figure_name= params[3]
community_method=params[4]
#fitxer_distancies = 'distmat_9014.dat'
#fitxer_corr = 'corr_9014_CA.dat'
#fitxer_distancies = 'matriu.dat'
#fitxer_corr = 'matriu_corr.dat'

data_dist = readFile(fitxer_distancies)
table_dist=tabularText(data_dist)

data_corr = readFile(fitxer_corr)
table_corr=tabularText(data_corr)

#print table_dist
#print table_corr



fila=0
distance=6.0

# Miro quins residues estan a menys de 5 angstroms de mitjana i miro el valor de correlacio que els correspon. Creo un dictionary que conte els valors de dist, corr i numeros de residues
valors=dict()

valors={'nodes':[],'edges':[],'weights':[], 'residues':[]}

while fila<(len(table_dist)):
	
	columna=0
	
	valors["nodes"].append(fila)
	valors["residues"].append(str(fila+1))
	while columna<len(table_dist[fila]):
		if (float(table_dist[fila][columna])<float(distance) and (fila < columna)):
			valors["edges"].append([fila,columna])
			valors["weights"].append(math.fabs(-math.log10(math.fabs(float(table_corr[fila][columna])))))
		#			print fila,columna,table_corr[fila][columna],math.fabs(-math.log10(math.fabs(float(table_corr[fila][columna]))))
		columna=columna+1
	fila=fila+1

#print valors

# Faig un arxiu pel pymol per poder visualitzar sobre l'estructura els residues que estan en contacte

#i=0
#
#text="hide everything\n"
#text=text+"show cartoon\n"
#
#while (i<len(valors["residue_i"])):
#	if (valors["residue_i"][i] < valors["residue_j"][i]):
#		text=text+"show spheres,"+str(valors["residue_i"][i]+1)+"/ca \n"
#		text=text+"show spheres,"+str(valors["residue_j"][i]+1)+"/ca \n"
#	i=i+1
#
#
#writeFile("pymol_select_residues",text)


#print valors["weights"]

# Genero un graph amb la info recopilada

g = igraph.Graph()

# assigno dades al graph

g.add_vertices(valors["nodes"])
g.add_edges(valors["edges"])
g.es["weight"]=valors["weights"]

g.vs["residues"]=valors["residues"]


#print g.is_weighted()
#
#print g.es["weight"]
#print g.vs["name"]
#print g.vs["residues"]
#
#print g.get_adjacency(2,"weight",0,False)
#
#print g

print "******Starting to compute shortest paths****"

def computeShortestPaths(graph,nodes,weights,output_type):
	shortest_paths=[]

	for node_idx_origin in nodes:
	
		node_idx_destination=node_idx_origin+2
	
		while (node_idx_destination < len(nodes)):
			shortest_paths.append(graph.get_shortest_paths(node_idx_origin,node_idx_destination, weights, "OUT",output_type))
			node_idx_destination=node_idx_destination+1

	return(shortest_paths)


shortest_paths=computeShortestPaths(g,valors["nodes"],valors["weights"],"vpath")

sp_graph = igraph.Graph()

sp_graph.add_vertices(valors["nodes"])

edges_already_added = [[0 for col in range(len(valors["nodes"]))] for row in range(len(valors["nodes"]))]


list_edges_shortest_paths=[]
list_edges=[]



for path_vertices in shortest_paths:
	path_vertices=path_vertices[0]
	
	for i in range(len(path_vertices)-1):
		num_edges = edges_already_added[path_vertices[i]][path_vertices[i+1]]
		if( num_edges == 0):
			list_edges_shortest_paths.append([path_vertices[i],path_vertices[i+1]])
		edges_already_added[path_vertices[i]][path_vertices[i+1]] = num_edges + 1
		edges_already_added[path_vertices[i+1]][path_vertices[i]] = num_edges + 1




sp_graph.add_edges(list_edges_shortest_paths)


#Genero matriu pesos per representar els camins mes curts

list_weights_shortest_paths=[]
list_weights_shortest_paths_unweighted=[]
labels_weighted_nodes=[""] * len(valors["nodes"])
vertex_sizes=["0"] * len(valors["nodes"])


# Funcio per Buscar element maxim d una matriu

def GetMaximumValueMatrix(matrix):
	maximum_value=0
	
	for list in matrix:
		for value in list:
			if (value > maximum_value):
				maximum_value=value
	return(maximum_value)

maximum_weight=float(GetMaximumValueMatrix(edges_already_added))


#print "****************************************"
#print "Maxim weight trobat per tots els edges que composen matriu: "+ str(maximum_weight)


for edge in list_edges_shortest_paths:
	reps = edges_already_added[edge[0]][edge[1]]/maximum_weight
	list_weights_shortest_paths_unweighted.append(edges_already_added[edge[0]][edge[1]])
	if(reps > 0.3):
		list_weights_shortest_paths.append(reps*10)
		labels_weighted_nodes[edge[0]] = valors["residues"][edge[0]]
		labels_weighted_nodes[edge[1]] = valors["residues"][edge[1]]
		vertex_sizes[edge[0]] = "8"
		vertex_sizes[edge[1]] = "8"
	else:
		list_weights_shortest_paths.append(0)


print "************Shortest paths DONE*************"

#print "****************************************"
#print "Llista de pesos pel shortest paths unweighted:"
#print list_weights_shortest_paths_unweighted
#print "****************************************"
#print "Numero de pesos pel shortest paths unweighted:"
#print len(list_weights_shortest_paths_unweighted)
#
#print "****************************************"
#print "Llista de pesos pel shortest paths weigthed:"
#print list_weights_shortest_paths
#print "****************************************"
#print "Numero de pesos pel shortest paths weighted:"
#print len(list_weights_shortest_paths)
#
#
#print "****************************************"
#print "Llista de labels pel shortest paths :"
#print labels_weighted_nodes
#print "****************************************"
#print "Numero de labels pel shortest paths :"
#print len(labels_weighted_nodes)

#writeFile(list_weights_shortest_paths,figure_name+"weights_shortest_paths")
#writeFile(list_weights_shortest_paths_unweighted,figure_name+"weights_shortest_paths_unweighted")
#writeFile(labels_weighted_nodes,figure_name+"labels_shortest_paths")
#print list_weights_shortest_paths
#print list_weights_shortest_paths_unweighted
#print labels_weighted_nodes

############################################################


#Calculo el numero d'arestes mitjanes per totes les parelles de nodes del graf g
average_path_length=g.average_path_length()

print "****************************************"
print "Calculo el average path length:"
print average_path_length

#Calculo ara la contribucio de cada residue k a la network recalculant de nou el CPL despres de treure el node k de la network
betweenness=g.betweenness(weights=g.es["weight"])

print "****************************************"
print "Calculo ara el betweeness de cada vertex del network"
print betweenness


# Funcio que donada una llista de pesos per cada edge i els edges que formen part del cami, calcula la suma dels pesos dels edges que formen part del path

#def lengthPath(epath,weight_list):
#	length=0
#	for edge_ID in epath:
#		length=length+weight_list[edge_ID]
#
#	return(length)
#
#
#
#
#shortest_edge_paths=computeShortestPaths(g,g.vs["name"],g.es["weight"],"epath")
#
#print len(shortest_edge_paths)
#
#length_shortest_edge_paths=[]
#
#for epath in shortest_edge_paths:
#	length_shortest_edge_paths.append(lengthPath(epath[0],g.es["weight"]))
#
#
#CPL=meanVarStd(length_shortest_edge_paths)[0]
#
#print "****************************************"
#print "El valor del CPL es:"
#print CPL

## Aqui es va traient un node del graf i es recalcula el CPL
#
#CPL_rem_list=[]
#
#for index_to_remove in range(len(g.vs)):
#	new_indices_list=[]
#	if(index_to_remove>0):
#		new_indices_list = new_indices_list + list(xrange(index_to_remove))
#	if(index_to_remove<len(g.vs)-1):
#		new_indices_list = new_indices_list + list(xrange(index_to_remove+1,len(g.vs)))
#
#	subgraph = g.subgraph(new_indices_list)
#
#	shortest_edge_paths=computeShortestPaths(subgraph,subgraph.vs["name"],subgraph.es["weight"],"epath")
#
#	length_shortest_edge_paths=[]
#	for epath in shortest_edge_paths:
#		length_shortest_edge_paths.append(lengthPath(epath[0],subgraph.es["weight"]))
#
#
#	CPL_rem=meanVarStd(length_shortest_edge_paths)[0]
#	CPL_rem_list.append(CPL_rem)
##	print CPL_rem
#
#print "Llista que conte tots els CPL_rem:"
#print CPL_rem_list
#
#
#print("Number of weiths of original graph: " + str(len(g.es["weight"])))
#print("Number of weiths of subgraph: " + str(len(subgraph.es["weight"])))
#
#
##Calculo ara el score Zk:
#
#CPL_rem_mean,CPL_rem_variance,CPL_rem_std=meanVarStd(CPL_rem_list)
#
#score_zk=[]
#
#for CPL in CPL_rem_list:
#	zk=(CPL-CPL_rem_mean)/CPL_rem_std
#	score_zk.append(zk)
#
#print "Calculo els scores zk per cada residue"
#print score_zk
#print "El numero total d scores es: "+str(len(score_zk))
#
#print "Els residues son:"
#print valors["residues"]
#
#figure()
#title("Score Zk")
#xlabel('Residue'); ylabel('Score Zk')
#plt.plot(valors["residues"],score_zk, linewidth=0, marker='o', color='#48D1CC')
#savefig("./score_zk.png")
#
#


print "******Starting to find communities ****"

###############################################################

#A partir d'aqui apliquem els diferents metodes per trobar communitats : fastgreedy, edge_betweness, walktrap. Veure documentacio igraph per veure dif entre algoritmes

if (community_method == "edge_betweness"):
	VertexDendrogram = g.community_edge_betweenness(None, False, "weight")

if (community_method == "walktrap"):
	VertexDendrogram=g.community_walktrap( "weight", 4)

if (community_method == "fastgreedy"):
	VertexDendrogram=g.community_fastgreedy( "weight")

#vertexClustering = g.community_infomap("weight")
community_graph = VertexDendrogram.as_clustering(VertexDendrogram.optimal_count).cluster_graph(None,"sum")


print "**********Communities Found*********"

community_vertices_weights = []
community_vertices_names = []

list_community_vertices_names=[]


for graph in VertexDendrogram.as_clustering(VertexDendrogram.optimal_count).subgraphs():
	#	print graph.vs["residues"]
	community_vertices_weights.append(len(graph.vs["residues"]))
	community_vertices_names.append(' '.join(graph.vs["residues"]))
	list_community_vertices_names.append(graph.vs["residues"])


community_graph.vs["size"] = community_vertices_weights
#community_graph.vs["label"] = community_vertices_names




# Genero fitxer per pymol per veure les comunitats i assigno el mateix color del pymol als nodes del graph
i = 0

community_vertices_colors=[]
community_vertices_node_names=[]

text="hide everything\n"
text=text+"show ribbon\n"
color_python=[]

while (i< len(community_vertices_names)):
	name_replaced = re.sub('\s+', '+', community_vertices_names[i])
	text=text+"select node "+str(i)+", resi "+name_replaced+"\n"
	color_pymol=get_color(i)
	node_label=""
	string_color="rgb("
	el_count = 0
	for el in color_pymol:
		if (el_count < 2):
			string_color=string_color+str(el)+","
		else:
			string_color=string_color+str(el)
		el_count = el_count + 1
	string_color=string_color+")"
	community_vertices_colors.append(string_color)
	text=text+"set_color c"+str(i)+"="+str(color_pymol)+"\n"
	text=text+"color c"+str(i)+", (node_"+str(i)+")\n"
	text=text+"label first resn ala in node_"+str(i)+", "+str(i)+"\n"
	node_label=str(i)
	community_vertices_node_names.append(node_label)
	
	i=i+1


#print community_vertices_colors

writeFile("pymol_community_analysis_"+figure_name,text)

# Assigno colors al graph

community_graph.vs["color"] = community_vertices_colors
community_graph.vs["label"] = community_vertices_node_names
community_graph.vs["label_dist"]= 0
#community_graph.vs["label_size"]= 20
community_graph.vs["label_color"]= "white"
#community_graph.vs["bbox"] = (1000, 1000)
#community_graph.vs["margin"] = 200



visual_style = {}
visual_style["vertex_size"] = community_vertices_weights
visual_style["edge_width"] = multiplyListValues(community_graph.es["weight"],10)

#Graph community


igraph.plot(community_graph,figure_name+".png",**visual_style)

#Graph de shortest path

visual_style = {}
visual_style["vertex_size"] = vertex_sizes
visual_style["edge_width"] = list_weights_shortest_paths
visual_style["vertex_label"] = labels_weighted_nodes
visual_style["vertex_label_size"] = "12"
visual_style["vertex_label_dist"] = "2"
#visual_style["vertex_label_angle"] = "15"



colors_weighted_nodes=["rgb(255,255,255)"] * len(valors["nodes"])

# Assigno el mateix color als nodes de graph de shortest path que als de community

i=0
for label in labels_weighted_nodes:
	if (label != ""):
		node=0
		for list in list_community_vertices_names:
			for element in list:
				if (label == element):
					colors_weighted_nodes[i] =community_vertices_colors[node]
			node=node+1
	i=i+1



sp_graph.vs["color"] = colors_weighted_nodes


igraph.plot(sp_graph,figure_name+"_shortest_paths.png",**visual_style)


text=""
text=text+"hide everything\n"
text=text+"show cartoon\n"
text=text+"set stick_color, black\n"
text=text+"bg_color white\n"

text_output=""
text_output=text_output+" Residue1 - Residue2 - Weight Path\n"
group=""

for index in range(len(list_weights_shortest_paths)):
	edge = list_edges_shortest_paths[index]
	weight = list_weights_shortest_paths[index]

	if weight > 0 :
#		print str(edge[0]+1)+" "+str(edge[1]+1)+" "+str(weight)

		text_output=text_output+str(edge[0]+1)+" "+str(edge[1]+1)+" "+str(weight)+"\n"
		
		text=text+"create test"+str(index)+", name ca in resi "+str(edge[0]+1)+"+"+str(edge[1]+1)+"\n"
		text=text+"show spheres, test"+str(index)+"\n"
		text=text+"set sphere_scale,"+str(0.1+weight/10)+", test"+str(index)+"\n"
		text=text+"bond name ca in resi "+str(edge[0]+1)+", name ca in resi "+str(edge[1]+1)+"\n"
		text=text+"set stick_radius,"+str(weight/10)+", test"+str(index)+"\n"
		text=text+"show sticks, test"+str(index)+"\n"
		group=group+"test"+str(index)+" "

print "****************************************"

text=text+"group PATH, "+group

writeFile("pymol_shortest_path_"+figure_name,text)
writeFile("output_shortest_path_"+figure_name+".out",text_output)






