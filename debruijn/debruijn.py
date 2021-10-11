#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib as plt 
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "REN Yani"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "yani.ren@etu.u-paris.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()

def read_fastq(fastq_file):
    """Read a fastq file and yield the sequence.
    :Parameters: path of the fastq file format.
    """
    with open(fastq_file, "r") as filin:
        lines = filin.readlines() 
        for i in range(1, len(lines), 4): 
            yield lines[i].strip()
           
def cut_kmer(read, kmer_size):
    """For every read cut and return kmer.
    :Parameters: 
    	read : a sequence in the fastq file.
    	kmer_size : the size of kmer. 
    """
    for i in range(0, len(read),1): 
        if len (read[i:i+kmer_size]) == kmer_size: 
            yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Read the fastq file and cut every sequence in kmer in order to construct a dictionnary of kmer.
    :Parameters: 
    	fastq_file: path of fasta format file, file to read 
    	kmer_size : int, the size of kmer 
    :Return: 
    	kmer_dict : dictionnary of kmer
    """ 
    kmer_dict = {}
    # To read the fastq file
    read_output = list(read_fastq(fastq_file))

    # To cut every seq in k-mer 
    for read in read_output: 
        cut_kmer_output = cut_kmer(read, kmer_size)
    
    # To construct the kmer_dict 
    for kmer in list(cut_kmer_output): 
        if kmer not in kmer_dict: 
            kmer_dict[kmer]=1
        else : 
            kmer_dict[kmer]+=1
    return kmer_dict

def build_graph(kmer_dict):
    """To construct graph with prefix and suffix kmers as nodes, from kmer_dict.
    :Parameters: 
	kmer_dict : dictionnary of kmer
    :Return:
    	G : graph build using kmer_dict
    """
    G = nx.DiGraph()
    for kmer, poids in kmer_dict.items(): 
        G.add_edge(kmer[:-1], kmer[1:], weight=poids)
    return(G)

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Return graph without some path to remove.
    :Parameters: 
    	graph : given graph 
    	path_list : contains path to remove
    	deletry_entry_node : boolean that indicates if entry node will be removed.
    	delete_sink_node : boolean that indicates if sink node will be removed.
    """
    for path in path_list: 
    	if delete_entry_node == True and delete_sink_node == True: 
    	    graph.remove_nodes_from(path)
    	elif delete_entry_node == True and delete_sink_node == False: 
    	    graph.remove_nodes_from(path[:-1])
    	elif delete_entry_node == False and delete_sink_node == True:
    	    graph.remove_nodes_from(path[1:])
    	else:
    	    graph.remove_nodes_from(path[1:-1])
    return(graph)

def std(data):
    """Compute and return the standard deviance value of given data."""
    return(statistics.stdev(data))
    
def select_best_path(graph, path_list, path_length, weight_avg_list,delete_entry_node=False, delete_sink_node=False):
    """Browse paths in 'path_list' and select the best path based on the best average weigth, then on the best length value of the path. If average weigth and path length are same, the best path is choose randomly. Return a graph with best path. 
    :Parameters:
    	graph : graph to treat 
    	path_list : contains paths to compare for choosing the best path
    	path_length : contains the length of path
    	weigth_avg_list : contains the average of path
    	deletry_entry_node : boolean that indicates if entry node will be removed.
    	delete_sink_node : boolean that indicates if sink node will be removed.
    """
    best_path_id = 0
    for i in range(1, len(weight_avg_list)): 
        if std([weight_avg_list[best_path_id] , weight_avg_list[i]])>0:	
    	    if weight_avg_list[best_path_id]> weight_avg_list[i]: 
    	        best_path = path_list[best_path_id]
    	    elif weight_avg_list[best_path_id] < weight_avg_list[i]:
    	        best_path_id = i
    	        best_path = path_list[i]
        elif std([weight_avg_list[best_path_id] , weight_avg_list[i]])==0:
            if std([path_length[best_path_id], path_length[i]])>0:
                if path_length[best_path_id] > path_length[i]: 
                    best_path = path_list[best_path_id] 
                elif path_length[best_path_id] < path_length[i]: 
                    best_path_id = i 
                    best_path = path_list[i]                  
            elif std([path_length[best_path_id], path_length[i]])==0:
                randm = randint(0,1) 
                if randm==0: 
            	    best_path = path_list[best_path_id]
                else: 
            	    best_path_id = i 
            	    best_path = path_list[i]             	    
    path_list = [path for path in path_list if path!=best_path]        	    
    graph=remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    
    return(graph)       
           
             
def path_average_weight(graph, path):
    """Compute and return the average weight of a given path in a graph.
    :Parameters: 
        graph : graph to treat. 
   	path : a path of the graph that we want to compute average weight.
    """
    weight_avg=statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])
    return(weight_avg)


def solve_bubble(graph, ancestor_node, descendant_node):   
    """Determine all possible path in a graph between two nodes. 
    :Parameters: 
    	ancestor_node : node in common of descendant node.
    	descendant_node : successor node of ancestor node. 
    """
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    path_length = []
    weight_avg_list = []
    print("PATH LIST ",len(path_list), path_list)
    for path in path_list:
        print("PATH ", path, len(path))
        path_length.append(len(path))   
        weight_avg=path_average_weight(graph, path)
        weight_avg_list.append(weight_avg)      
    graph=select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node=False, delete_sink_node=False)
                     
    return(graph)
   
    	
def simplify_bubbles(graph):
    """Detect bubbles in graph. Make a recursive call of itself for remove bubbles and return a graph without bubbles.""" 
    bubble = False
    for node in list(graph.nodes):        
        predecessors_nodes = list(graph.predecessors(node))        
        if len(predecessors_nodes)>1: 
            for i in range(len(predecessors_nodes)-1):
                for j in range(i+1,len(predecessors_nodes)): 
                    node1=predecessors_nodes[i]
                    node2=predecessors_nodes[j]
                    ancestor_node = nx.lowest_common_ancestor(graph, node1, node2)
                 
                    if ancestor_node!=None:                
                        bubble = True                  
                        break   
            if bubble:           
    	        graph = simplify_bubbles(solve_bubble(graph, ancestor_node, node))
    	        break                   	
    return(graph) 
    
    
def solve_entry_tips(graph, starting_nodes):
    """ Detect entry tips and return a graph without entry tips.
    :Parameters: 
    	graph : given graph.
    	starting_nodes : nodes of given graph that don't have predecessors. 
    """
    for node in list(graph.nodes): 
        predecessors_nodes = list(graph.predecessors(node)) 
        node_start_list = []    
        if len(predecessors_nodes)>1: 
            for node_start in starting_nodes:
            	has_path = nx.has_path(graph, node_start, node)
            	if has_path == True: 
            	    node_start_list.append(node_start)
            
            if len(node_start_list)==2:               
                path_list = []
                for node_start in node_start_list:                                
                    path_list.extend(list(nx.all_simple_paths(graph, source=node_start, target=node)))                 
                if std([path_average_weight(graph, path_list[0]), path_average_weight(graph, path_list[1])])>0: 
                    if path_average_weight(graph, path_list[0])> path_average_weight(graph, path_list[1]):
                        best_path = path_list[0]
                    else:
                        best_path = path_list[1]
                elif std([path_average_weight(graph, path_list[0]), path_average_weight(graph, path_list[1])])==0:
                    print(type(path_list[0]))
                    if std([len(path_list[0]), len(path_list[1])])>0: 
                        if len(path_list[0]) > len(path_list[1]):
                            best_path = path_list[0]
                        elif len(path_list[0]) < len(path_list[1]):
                            best_path = path_list[1]
                        else:
                            best_path = random.choice([path_list[0], path_list[1]])               
                path_list = [path for path in path_list if path!=best_path]               
                graph = remove_paths(graph, path_list, delete_entry_node = True, delete_sink_node = False)
                break
    return(graph)     

def solve_out_tips(graph, ending_nodes):
    """Detect out tips and return a graph without out tips for a given graph.
    :Parameters: 
    	graph : given graph.
    	ending_nodes : nodes of given graph that don't have sucesssors.  
    """
    for node in list(graph.nodes): 
        successors_nodes = list(graph.successors(node)) 
        node_end_list = []    
        if len(successors_nodes)>1: 
            for node_end in ending_nodes:                
                has_path = nx.has_path(graph, node, node_end)
                if has_path == True: 
            	    node_end_list.append(node_end)            
            if len(node_end_list)==2:
                path_list = []
                for node_end in node_end_list:                                
                    path_list.extend(list(nx.all_simple_paths(graph, source=node, target=node_end)))                    
                path_length = []
                weigth_avg_list = []
                if std([path_average_weight(graph, path_list[0]), path_average_weight(graph, path_list[1])])>0: 
                    if path_average_weight(graph, path_list[0])> path_average_weight(graph, path_list[1]):
                        best_path = path_list[0]
                    else:
                        best_path = path_list[1]
                elif std([path_average_weight(graph, path_list[0]), path_average_weight(graph, path_list[1])])==0:
                    print(type(path_list[0]))
                    if std([len(path_list[0]), len(path_list[1])])>0: 
                        if len(path_list[0]) > len(path_list[1]):
                            best_path = path_list[0]
                        elif len(path_list[0]) < len(path_list[1]):
                            best_path = path_list[1]
                        else:
                            best_path = random.choice([path_list[0], path_list[1]])
                                    
                path_list = [path for path in path_list if path!=best_path]               
                graph = remove_paths(graph, path_list, delete_entry_node = False, delete_sink_node = True)
                break
    return(graph)

def get_starting_nodes(graph):
    """Search and return starting nodes for a given graph."""
    starting_nodes = list()
    # To search for each node if they have predecessors nodes, if not, they are considered as starting nodes. 
    for node in list(graph.nodes):
        node_pred = list(graph.predecessors(node))
        #print(node_pred)
        if len(node_pred)==0:            
            starting_nodes.append(node)
    return(starting_nodes)
            

def get_sink_nodes(graph):
    """Search and return ending nodes for a given graph."""
    sink_nodes = list()
    # To search for each node if they have predecessors nodes, if not, they are considered as input nodes. 
    for node in list(graph.nodes):
        node_next = list(graph.successors(node))
        #print(node_pred)
        if len(node_next)==0: 
            sink_nodes.append(node)
    return(sink_nodes)

def get_contigs(graph, starting_nodes, ending_nodes):
    """Search for all possible path and find contig. Return a list of found contigs.
    :Parameters: 
    	graph : graph to analyze
    	starting_nodes : nodes that don't have predecessors.
    	ending_nodes : nodes that don't have sucessors. 
    :Return: 
    	contigs_list : contains all contig in the graph. 
    """
    contigs_list = []   
    all_contig = [list(nx.all_simple_paths(graph, source=node_start, target=node_end)) for node_start in starting_nodes for node_end in ending_nodes if nx.has_path(graph, node_start, node_end)==True]
    for contig in all_contig:
        contig_deb = contig[0][0] 
        contig_fin = [contig[0][i][-1] for i in range(1,len(contig[0]))] 
        contig_final = contig[0][0]+"".join(contig_fin)      
        contigs_list.append((contig_final, len(contig_final)))    
  
    return(contigs_list) 

def save_contigs(contigs_list, output_file):
    """Write contig in a text file format.
    :Parameters: 
	contigs_list : contains all contig to write. 
	output_file : path of the output file.
    """
    with open(output_file,"w") as filout: 
        for i in range(len(contigs_list)):
            entete =">contig_{} len={}\n".format(i,int(contigs_list[i][1]))
            text = str(contigs_list[i][0])
            filout.write(entete)
            text = str(fill(text))      
            filout.write(text)

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
#    with open(graph_file, "wt") as save:
#            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    fastq_file = args.fastq_file
    kmer_size = args.kmer_size
    output_file = args.output_file
    graphimg_file = args.graphimg_file
    print("ARGUMENTS : ", args)
     
    # Construct kmer dictionnary 
    kmer_dict = build_kmer_dict(fastq_file, kmer_size)
   
    # Built graph 
    graph = build_graph(kmer_dict)
	
    # Browse graph for searching starting and ending nodes, to determine contig 
    starting_nodes = get_starting_nodes(graph)
    print("INPUT : ", starting_nodes)
    sink_nodes = get_sink_nodes(graph)
    print("OUPUT ", list(sink_nodes))
    
    # Solve bubbles if there are bubbles
    graph = simplify_bubbles(graph)
    
    # Solve tips if there are tips  
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, sink_nodes)
    
    # Write contigs in texte format file
    contigs_list = get_contigs(graph, starting_nodes, sink_nodes)
    print("CONTIGS : ",contigs_list)
    save_contigs(contigs_list, output_file)

    
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    #if args.graphimg_file:
        #draw_graph(graph, args.graphimg_file)
    #Save the graph in file
    #if args.graph_file:
        #save_graph(graph, args.graph_file)



if __name__ == '__main__':
    main()
