# miR-202-in-Oryzias_latipes


This repository contains a set of Python and shell scripts designed to analyze and predict miR-202's interactions with targets. The scripts aim to process sequence data, create heatmaps, create newick files and sub-trees based on relationships between our species, and visualize the prediction of the interactions between the miRNA and the targeted gene.

## Getting Started

### Prerequisites

Make sure you have the following software and libraries installed on your system:

- Python 3.7 or later
- Biopython 1.79
- Pandas 1.0.0 to 1.3.0
- NumPy 1.19.0 to 1.21.0
- Matplotlib 3.3.0 to 3.4.2
- Seaborn 0.10.0 to 0.11.1
- argparse

You can use the provided `requirements.txt` file to install the necessary libraries using the following command:

pip install -r requirements.txt

Or in command lines : 
- pip install pandas
- pip install metaplot
- pip install numpy
- pip install Biopython
- pip install seaborn
- pip install ete3
- pip install argparse 


### Usage

The scripts can be executed from the command line. Below are the main commands along with their usage:

1. **script_launcher**:

./script_launcher GENE_NAME


Example: 

./script_launcher mst1



This script performs a series of steps, including creating files for each species, executing RNAduplex for each file, and producing a Newick file containing only the species of interest. The Newick file is then used as input for the next script.

2. **script_Dataframe_treatment**:

./script_Dataframe_treatment GENE_NAME FINAL_NEWICK_FILE


Example:

./script_Dataframe_treatment mst1 final_newick_for_mst.nwk


This script reads the merged_RNAduplex file, crosses the information with the given Newick file, and generates a dataframe containing data for common species. The dataframe is further processed to handle species with multiple miRNA sequences in their genomes (_2, _3, etc.) and create a HEATMAP of the predictions of the miRNA/target interactions.

Please note that these scripts are specifically designed for miRNA and target interaction analysis. THe newick files are either in the "genus_species" or "genus.species" formats (with no informations about the higher clades). For other applications, additional modifications may be necessary.

## Workflow Details

### script_launcher

1. **script_pair_candidate_and_orthologous.py**:
- Creates a file for each species, pairing each miRNA sequence with its orthologous sequence in that species.

2. **script_Aligner_RNAduplex.sh**:
- Executes RNAduplex for each file, creating additional files (_2, _3, etc.) for species with multiple miRNA sequences that interact with the target.

3. **merge.py**:
- Merges RNAduplex predictions for all species into a single file, "unfiltered_merged_RNAduplex".

4. **script_pre-treatment.py**:
- Basically, this script Processes the "unfiltered_merged_RNAduplex" file, adding extensions (_2, _3, etc.) for species with multiple miRNA sequences. Duplicates with identical miRNA sequences are removed, and the result is stored in "merged_RNAduplex".

5. **script_newick_creator.py**:
- Creates a Newick file containing only the species of interest, which can be rearranged as required.


### script_Dataframe_treatment

1. **script_Dataframe_creator.py**:
- Reads the "merged_RNAduplex" and Newick files to create a dataframe with information for common species. Note that here, the duplicates species do not have an index value yet since "Oryzias_latipes_2" for example is not found in the newick file, yet have to be given the same index as the original Oryzias_latipes.

2. **script_filter_dataframe.py**:
- Filters the dataframe, assigning the same index to species with extensions (_2, _3, etc.) as their original species. Also, creates a file with species names and their corresponding indexes.

3. **script_for_scripts.py**:
- Reads the filtered dataframe produced just above, and writes the new species indexes to "script_finalHeatmap.py", which is used to generate heatmaps.

4. **script_finalHeatmap.py**:
- Produces two heatmaps of miRNA-target interactions and target-miRNA interactions in SVG format, considering the phylogenetic relationships between species given in the previously produced newick file.

## Note

Please remember that these scripts are designed for specific miRNA and target interaction analysis. Adjustments and modifications may be required for other applications. If you encounter any issues or have questions, feel free to reach out for assistance. Happy analyzing!
