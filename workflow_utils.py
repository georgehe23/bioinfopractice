import subprocess
import networkx as nx
import matplotlib.pyplot as plt

from networkx.drawing.nx_pydot import graphviz_layout



def setup_custom_venv(requirements_file, custom_venv_name):
    script_path = 'bin/setup.sh'
    with open(script_path, 'r') as script:
        modified_script = script.read().replace('custom_venv', custom_venv_name).replace('requirements.txt', requirements_file)
    subprocess.run(['bash', '-c', modified_script], check=True, text=True)



def draw_evo_tree(T):
    """"
    Source: https://github.com/melbournebioinformatics/COMP90014_2024/raw/master/assignments/src/A2utils.py
    The relevant visualisation function has been extracted from the assignment utilities file.
    It plots a phylogenetic tree using the NetworkX package.
    Case samples are labelled in red and previous sample are labelled in black. 
    Samples and cases are grouped in common clades according to their last common ancestor. 

    Input: a phylogenetic tree in the form of a NetworkX graph.
    Output: an interactive window that displays the phylogenetic tree graph.
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    pos = graphviz_layout(T, prog="dot")
    nx.draw_networkx_nodes(T, pos, node_color='#636363', node_size=1)
    labels = {n: T.nodes[n]['label'] for n in T.nodes() if T.out_degree(n) == 0}
    
    # print(labels)
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