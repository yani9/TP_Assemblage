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
    # To construct graph with prefix and suffix kmers as nodes, from kmer_dict.
    G = nx.DiGraph()
    for kmer, poids in kmer_dict.items(): 
        G.add_edge(kmer[:-1], kmer[1:], weight=poids)
  
    return G

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

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
    all_contig = list() 
    for node_start in starting_nodes: 
        for node_end in ending_nodes:
            contig = list(nx.all_simple_paths(graph, source=node_start, target=node_end))
            #all_contig.append(contig)
    return(contig, len(contig))

def save_contigs(contigs_list, output_file):
    pass


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
    fastq_file = "/home/sdv/m2bi/yren/Documents/Assemblage/debruijn-tp/data/eva71_two_reads.fq"
    read_output = list(read_fastq(fastq_file))
    print(read_output)

    # To construct kmers 
    kmer_size = 22
    #read = "TCAGAGCTCTAGAGTTGGTTCTGAGAGAGATCGGTTACTCGGAGGAGGCTGTGTCACTCATAGAAGGGATCAATCACACCCACCACGTGTACCGAAACAA"
    for read in read_output: 
        cut_kmer_output = cut_kmer(read, kmer_size)
    print(list(cut_kmer_output), len(list(cut_kmer_output)))
    

    # To construct kmer dico 
    kmer_dict = build_kmer_dict(fastq_file, kmer_size)
    print(kmer_dict)

    # To built graph 
    graph = build_graph(kmer_dict)
    #print(type(graph))

    # To browse graph for searching input and output nodes, to determine contig 
    node_input = get_starting_nodes(graph)
    print("INPUT : ", list(node_input))
    node_output = get_sink_nodes(graph)
    print("OUPUT ", list(node_output))

    contig = get_contigs(graph, node_input, node_output)
    print(contig)

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
