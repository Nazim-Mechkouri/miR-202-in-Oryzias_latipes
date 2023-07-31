import sys, os, re
import Bio
import ete3
from Bio import Phylo
from ete3 import Tree

newick_file = sys.argv[1]
# Load the original Newick tree file
original_tree = Phylo.read(newick_file, 'newick')
print(original_tree)

species_list_file = sys.argv[2]

# Read the Newick tree from file
t = Tree(newick_file)
#print(tree)
species_list = []
# Read the species list from file
with open(species_list_file, 'r') as f:
    for line in f:
        liner = line.split()[0]
        species_list.append(liner)

species_list = [element.replace("\n","") for element in species_list]
species_list = [element.replace("_",".") for element in species_list]
species_list = [element.replace(">","") for element in species_list]
species_set = set(species_list)

print(len(species_list),"	espèces totales de départ")

# Filter species_list to include only the species present in the newick tree
species_list = [species for species in species_list if species in t.get_leaf_names()]

print(len(species_list),"	espèces communes entre celles de départ et le fichier newick utilisé")

subtree = t.prune(species_list, preserve_branch_length=True)
t.write(subtree, outfile="small_cladogram.newick")

