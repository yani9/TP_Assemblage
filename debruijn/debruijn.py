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

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
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
    with open(fastq_file, "r") as filin:
        lines = filin.readlines() 
        for i in range(1, len(lines), 4): 
            yield lines[i].strip()
           
def cut_kmer(read, kmer_size):
    for i in range(0, len(read),1): 
        if len (read[i:i+kmer_size]) == kmer_size: 
            yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    # To read the fastq file
    read_output = list(read_fastq(fastq_file))

    # To cut every seq in k-mer 
    for read in read_output : 
        cut_kmer_output = cut_kmer(read, kmer_size)
    
    # To construct the kmer_dict 
    for kmer in list(cut_kmer_output): 
        if kmer not in kmer_dict: 
            kmer_dict[kmer]=1
        else : 
            kmer_dict[kmer]+=1

    return kmer_dict

def build_graph(kmer_dict):
    """To construct graph with prefix and suffix kmers as nodes, from kmer_dict."""
    G = nx.DiGraph()
    for kmer, poids in kmer_dict.items(): 
        G.add_edge(kmer[:-1], kmer[1:], weight=poids)
    return(G)

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
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
    return(statistics.stdev(data))
    
def select_best_path(graph, path_list, path_length, weight_avg_list,delete_entry_node=False, delete_sink_node=False):
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
    weight_avg=statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])
    return(weight_avg)


def solve_bubble(graph, ancestor_node, descendant_node):   
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
    tip = False
    for node in list(graph.nodes): 
        predecessors_nodes = list(graph.predecessors(node)) 
        pred_node_start = []        
        if len(predecessors_nodes)>1:
            for pred_node in predecessors_nodes:
                while len(list(graph.predecessors(pred_node)))!=0: 
                    pred_node = graph.predecessors(pred_node)
                print("PRED ", pred_node)
                if pred_node in starting_nodes:
                    pred_node_start.append(pred_node)
            if len(pred_node_start)>1:
                print(pred_node_start, node)
                tip=True 
                break       
        if tip:
            for pred_node in pred_node_start :
                path_list = list(nx.all_simple_paths(graph, source=pred_node, target=node))
                #path_list = [path for path in path_list if path==node]
                path_length = []
                weigth_avg_list = []
                if str([path_average_weight(graph, path_list[0]), path_average_weight(graph, path_list[1])])>0: 
                    if path_average_weight(graph, path_list[0])> path_average_weight(graph, path_list[1]):
                        path_list = path_list[1]
                    elif path_average_weight(graph, path_list[0])<path_average_weight(graph, path_list[1]):
                        path_list = path_list[0]
                elif str(path_average_weight(graph, path_list[0]), path_average_weight(graph, path_list[1]))==0:
                    print(type(path_list[0]))
                    if std(len(path_list[0]), len(path_list[1]))>0: 
                        if len(path_list[0])>len(path_list[1]):
                            path_list = path_list[1]
                        elif len(path_list[0])<len(path_list[1]):
                             path_list = path_list[0]
                        else:
                            path_list = random.choice([path_list[0], path_list[1]])
                
                """
                for path in path_list:
                    path_length.append(len(path))
                    weight_avg = path_average_weight(graph, path)
                    weight_avg_list.append(weight_avg)
                """    
                graph = solve_entry_tips(remove_paths(graph, path_list, delete_entry_node = True,
                delete_sink_node = False))
                print(graph) 
    return(graph) 

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    node_input = list()
    #print("ALL NODES : ", list(graph.nodes))

    # To search for each node if they have predecessors nodes, if not, they are considered as input nodes. 
    for node in list(graph.nodes):
        node_pred = list(graph.predecessors(node))
        #print(node_pred)
        if len(node_pred)==0: 
            #print(node)
            node_input.append(node)

    return(node_input)
            

def get_sink_nodes(graph):
    node_output = list()
    #print("ALL NODES : ", list(graph.nodes))

    # To search for each node if they have predecessors nodes, if not, they are considered as input nodes. 
    for node in list(graph.nodes):
        node_next = list(graph.successors(node))
        #print(node_pred)
        if len(node_next)==0: 
            #print(node)
            node_output.append(node)

    return(node_output)

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs_list = []   
    all_contig = [list(nx.all_simple_paths(graph, source=node_start, target=node_end)) for node_start in starting_nodes for node_end in ending_nodes if nx.has_path(graph, node_start, node_end)==True]
    for contig in all_contig :
        contig_deb = contig[0][0] 
        contig_fin = [contig[0][i][-1] for i in range(1,len(contig[0]))] 
        contig_final = contig[0][0]+"".join(contig_fin)      
        contigs_list.append((contig_final, len(contig_final)))    
  
    return(contigs_list) 

def save_contigs(contigs_list, output_file):
    with open(output_file,"w") as filout: 
        for i in range(len(contigs_list)):
            entete =">contig_{} len={}\n".format(i,int(contigs_list[i][1]))
            text = str(contigs_list[i][0])
            filout.write(entete)
            text = str(fill(text)) 
            #print("To write : ",text)       
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
    #args = get_arguments()
    
    # To read fastqc file 
    #fastq_file  = "/home/sdv/m2bi/yren/Documents/Assemblage/debruijn-tp/data/eva71_hundred_reads.fq"
    #fastq_file = "/home/sdv/m2bi/yren/Documents/Assemblage/debruijn-tp/data/eva71_two_reads.fq"
    #fastq_file = "./data/eva71_two_reads.fq"
    fastq_file = "./data/eva71_hundred_reads.fq"
    read_output = list(read_fastq(fastq_file))
    print(read_output)

    # To construct kmers 
    kmer_size = 22
    #read = "TCAGAGCTCTAGAGTTGGTTCTGAGAGAGATCGGTTACTCGGAGGAGGCTGTGTCACTCATAGAAGGGATCAATCACACCCACCACGTGTACCGAAACAA"
    for read in read_output: 
        cut_kmer_output = cut_kmer(read, kmer_size)
#    print(list(cut_kmer_output), len(list(cut_kmer_output)))
    

    # To construct kmer dico 
    kmer_dict = build_kmer_dict(fastq_file, kmer_size)
#    print(kmer_dict)

    # To built graph 
    graph = build_graph(kmer_dict)
    #print(type(graph))

    # To browse graph for searching input and output nodes, to determine contig 
#    node_input = get_starting_nodes(graph)
    starting_nodes = get_starting_nodes(graph)
    print("INPUT : ", starting_nodes)
#    node_output = get_sink_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    print("OUPUT ", list(sink_nodes))

#    contigs_list = get_contigs(graph, node_input, node_output)
    contigs_list = get_contigs(graph, starting_nodes, sink_nodes)

    output_file = "save_contig_list.fasta"
    save_contigs(contigs_list, output_file)

    graph_1 = nx.DiGraph()
    graph_1.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    graph_1 = select_best_path(graph_1, [[1,2], [3,2]], [1, 1], [5, 10], delete_entry_node=True)
   
    graph_2 = nx.DiGraph()
    graph_2.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7), (7, 8)])
    graph_2 = select_best_path(graph_2, [[5, 6], [5, 7, 8]], [1, 2], [13, 10], delete_sink_node=True)
    graph_3 = nx.DiGraph()
    graph_3.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (2, 8), (8, 9),
                            (9, 5), (5, 6), (5, 7)])
    graph_3 = select_best_path(graph_3, [[2, 4, 5], [2, 8, 9, 5]],
                                         [1, 4], [13, 10])
                                         
                                         
    starting_nodes = get_starting_nodes(graph_3)
    sink_nodes = get_sink_nodes(graph_3)
    graph = simplify_bubbles(graph_3)
    print("RESULT ", list(graph)) 
    
    graph_1 = nx.DiGraph()
    graph_1.add_weighted_edges_from([(1, 2, 10), (3, 2, 2), (2, 4, 15), (4, 5, 15)])
    graph_1 = solve_entry_tips(graph_1, [1, 3])
    print("RESULT2 ", list(graph_1)) 
    
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
