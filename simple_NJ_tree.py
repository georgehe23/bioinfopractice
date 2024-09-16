# Written for Assignment 2, Question 3a (phylogenetics) for
# Semester 2, 2024 in Algorithms for Bioinformatics (COMP90014)
# at the University of Melbourne
"""Setting up the environment:


create python virtual environment
> python3 -m venv NJ_venv
> source NJ_venv/bin/activate


make sure the system has an available Tkinter installed:
> sudo apt-get install python3-tk
this is used for visualisation from the commandline, 
the appropriate command to run this run and visualise 
after the environment and data is setup is:
> python simple_NJ_tree.py


load required modules
> pip install biopython==1.84 
graphviz==0.20.3 
matplotlib==39.2 
networkx==3.3 
numpy==2.1.0 
pydot==3.0.1
tk==0.1.0
"""


import Bio
import Bio.Align
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle

from networkx.drawing.nx_pydot import graphviz_layout
from collections import defaultdict

matplotlib.use('TkAgg')


"""
The following files are required to complete this assignment. 
- https://github.com/melbournebioinformatics/COMP90014_2024/raw/master/assignments/data/ebola_glycoproteins.fasta
- https://github.com/melbournebioinformatics/COMP90014_2024/raw/master/assignments/data/ebola_glycoproteins_MSA.fasta

Download these and place in a subfolder called /data.
"""


# setting up the data and environment
EBOLA_SEQS_FILEPATH = "data/ebola_glycoproteins.fasta"
# I might be lucky and be able to upload the data files to my GitHub, and if so, you can copy it from there.


""""
The following python files are required to complete this assignment. 
- https://github.com/melbournebioinformatics/COMP90014_2024/raw/master/assignments/src/A2utils.py
The relevant visualisation function has been extracted and is listed below:
"""
def draw_evo_tree(T):
    fig, ax = plt.subplots(figsize=(10, 10))
    pos = graphviz_layout(T, prog="dot")
    nx.draw_networkx_nodes(T, pos, node_color='#636363', node_size=1)
    labels = {n: T.nodes[n]['label'] for n in T.nodes() if T.out_degree(n) == 0}
    print(labels)
    prev_labels = {k: v for k, v in labels.items() if not v.startswith('CASE')}
    case_labels = {k: v for k, v in labels.items() if v.startswith('CASE')}
    text1 = nx.draw_networkx_labels(T, pos, labels=prev_labels, font_size=8, font_family="sans-serif", verticalalignment='top')
    text2 = nx.draw_networkx_labels(T, pos, labels=case_labels, font_size=8,  font_color="red", font_family="sans-serif", verticalalignment='top')
    for _, t in text1.items():
        t.set_rotation('vertical')
    for _, t in text2.items():
        t.set_rotation('vertical')
    cstyle = 'bar,angle=-180,fraction=0'
    cstyle = 'angle,angleA=-90,angleB=180,rad=0'
    nx.draw_networkx_edges(T, pos, width=1, edge_color='#636363', arrowstyle='-', node_size=0, connectionstyle=cstyle, min_source_margin=0, min_target_margin=0)
    plt.tight_layout()
    plt.subplots_adjust(bottom=-0.4)
    plt.show()



"""
Below I create my phylogenetic tree using neighbour joining.
I do a pairwise global protein sequence alignment because
I was analysing the evolution of Ebola sequences to find
the best vaccine for this assignment question I was working on. 
This took me countless days and nights and I tried many sad methods, 
most of which included indexing indexes on indexes on indexes. 
I finally managed to simplify it to this and I don't know why 
I spent so much time on my other attempts. 
:(
Enjoy
"""
# load data
def load_sequences(filepath: str):
    from Bio import SeqIO
    seqs = []
    for seq in SeqIO.parse(filepath, "fasta"):
        seqs.append((seq.name, str(seq.seq)))
    return seqs
seqs = load_sequences(EBOLA_SEQS_FILEPATH)    


# main function
def build_evo_tree(sequences):
    G = nx.DiGraph()  # init variables
    node_id = 0
    matrix_dict = preprocessing(seqs)

    G, node_id, node_dict = set_leaf_nodes(G, matrix_dict, node_id)  # init graph (leaves)
    dist_matrix = get_distance_matrix(matrix_dict)  # all seqs dict: {seq:[list of strain names]}

    while len(dist_matrix) > 2:  # do Neighbour Joining from leaves
        pair = closest_pair(dist_matrix) 
        dist_matrix, matrix_dict, node_dict, node_id, G = add_pair_parent(pair, dist_matrix, matrix_dict, node_dict, node_id, G)
    pair = closest_pair(dist_matrix)  # add root node
    dist_matrix, matrix_dict, node_dict, node_id, G = add_pair_parent(pair, dist_matrix, matrix_dict, node_dict, node_id, G)

    return G



# sort sequences so they are referencable across data structures by index
def preprocessing(seqs):  
    all_dict = {name:seq for name,seq in seqs}  # unpack and sort seqs
    all_dict = dict(sorted(all_dict.items()))
    
    seqs_dict = defaultdict(list)  # seqs_dict = {seq:[list of strain names]}
    for name,seq in all_dict.items():
        seqs_dict[seq] += [name]

    return seqs_dict # all sublisted names are sorted



def get_distance_matrix(seq_dict):    
    aligner = Bio.Align.PairwiseAligner()  # use aligner to make a distance matrix for protein sequences
    aligner.mode = 'global'  # global protein sequence alignment of the same region
    pam30 = Bio.Align.substitution_matrices.load('PAM30')  # using PAM as this is global alignment
    aligner.substitution_matrix = pam30  # using PAM30 due to expected close similarity

    n = len(seq_dict)
    dist_matrix = np.zeros((n, n))  # init dist matrix

    # unpack sorted sequence from sequence dict
    unpack_all = sorted(seq_dict.keys())
    for i in range(n):  # pairwise sequence comparison of protein seqs using PAM30 
        for j in range(i + 1, n):
            alignment = aligner.align(unpack_all[i], unpack_all[j])
            score = -alignment.score  # use negative alignment similarity score as distance
            dist_matrix[i, j] = dist_matrix[j, i] = score

    # normalise the distance matrix
    if np.max(dist_matrix) != np.min(dist_matrix):
        dist_matrix = (dist_matrix - np.min(dist_matrix)) / (np.max(dist_matrix) - np.min(dist_matrix))
    
    return dist_matrix  # returns sorted distance matrix, in same order as entries in seq_dict



def set_leaf_nodes(graph, seq_dict, node_id):  # init graph with leaf nodes
    node_dict = {}

    for seq in seq_dict:  # each unique sequence has a list of related strains in seq_dict's defaultdict(list)
        for strain_name in range(len(seq_dict[seq])):
            
            node_dict.update({node_id:seq_dict[seq][strain_name]})  # track added nodes in a node_dict
            graph.add_node(node_id, label=seq_dict[seq][strain_name])  # add leaves to graph
            node_id += 1  # keep unique node id

    return graph, node_id, node_dict



"""NEIGHBOUR JOINING ALGORITHM"""
# find average distance each node is from other nodes
def node_avg_dist(i, num_nodes, dist_matrix):
    return np.mean([dist_matrix[i][j] for j in range(num_nodes) if i != j])



# find closeness proximity and average distance from other nodes for each node pair
def calculate_Qij(i, j, num_nodes, dist_matrix):  # returns the Qij score for the node pair
    return dist_matrix[i][j] - node_avg_dist(i, num_nodes, dist_matrix) - node_avg_dist(j, num_nodes, dist_matrix)



# find starting node pair through smallest Qij
def closest_pair(dist_matrix):
    best_Qij = np.inf
    n = len(dist_matrix)  # go through matrix and get Qij score per node pair

    for i in range(n):  # iterate through pairs and check Qij score
        for j in range(i + 1, n):
            Qij = calculate_Qij(i, j, n, dist_matrix)

            if Qij < best_Qij:  # set best Qij
                best_Qij = Qij
                pair = (i,j)

    return pair  # best pair returned as indexes of dist_matrix



def remove_pair(pair, matrix_dict, dist_matrix):
    pair = sorted(pair, reverse=True)  # remove later row/column first to maintain indexing

    cleaned_dict = {a: b for x, (a, b) in enumerate(matrix_dict.items()) if x not in pair}
    dist_matrix = np.delete(dist_matrix, pair, axis=0)  # remove rows
    dist_matrix = np.delete(dist_matrix, pair, axis=1)  # remove columns

    return dist_matrix, cleaned_dict



def update_distances(i, j, k, dist_matrix, matrix_dict, node_dict, node_id):
    new_dist_matrix = np.zeros((k + 1, k + 1))  # init new matrix
    new_dist_matrix[:k, :k] = dist_matrix  # copy old matrix into top left, to preserve sort order
    
    for m in range(k):  # add k distances from other nodes
        new_dist_matrix[k][m] = 0.5 * (dist_matrix[i][m] + dist_matrix[j][m] - dist_matrix[i][j])
        new_dist_matrix[m][k] = new_dist_matrix[k][m]

    matrix_dict.update({k:[node_id]}) # define new parent node in matrix_dict
    node_dict.update({node_id:node_id}) # define new parent node in overall node_dict

    return new_dist_matrix, matrix_dict, node_dict



def add_pair_parent(pair, dist_matrix, matrix_dict, node_dict, node_id, graph):
    i, j = pair  # unpack best i,j indexes 
    k = len(dist_matrix)  # newest row/column index
    dist_matrix, matrix_dict, node_dict = update_distances(i, j, k, dist_matrix, matrix_dict, node_dict, node_id)

    # get new edge lengths from i,j to the new internal node
    new_i_edge = 0.5 * (dist_matrix[i][j] + node_avg_dist(i, k, dist_matrix) - node_avg_dist(j, k, dist_matrix))
    new_j_edge = dist_matrix[i][j] - new_i_edge
    new_edges = (new_i_edge, new_j_edge)

    graph.add_node(node_id, label=node_id)  # add new parent node to graph
    
    # get names of indexed seqs in matrix_dict, and use them to get their node index in node_dict
    for x in range(len(pair)):  # for i and j edges
        for matrix_index, seq in enumerate(matrix_dict):  # index matrix dict columns
            if matrix_index == pair[x]:  # proceed with only index-matched dict entries
                    
                for matrix_name in matrix_dict[seq]:  # get names from matrix_dict {seq:[names]}
                    for node_index in node_dict:  # find matched matrix_dict name in node_dict
                        dict_name = node_dict[node_index]  
                        
                        if dict_name == matrix_name:  # cross-ref dicts with names and get node index
                            graph.add_edge(node_id, node_index, len=new_edges[x])  # add graph edges

    dist_matrix, matrix_dict = remove_pair(pair, matrix_dict, dist_matrix)  # remove from matrix/dict
    node_id += 1  # keep track of node ids added

    return dist_matrix, matrix_dict, node_dict, node_id, graph




G = build_evo_tree(seqs)

seqs = load_sequences(EBOLA_SEQS_FILEPATH)
G = build_evo_tree(seqs)
draw_evo_tree(G)
