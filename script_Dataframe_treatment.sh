#!/bin/bash

#My arguments
gene_name=$1
newick_file=$2

#My variables
gene_check=$1
newick_check=$1
new_directory=$gene_name"_candidates_combined_orthologous"


#My arguments
gene_name=$1
newick_file=$2




#My conditions
if [[ $# -eq 0 ]] ; then
    echo 'Aucun argument donné : Entrez les deux arguments ($1=nom de gène valide	    $2=fichier_newick)'
    exit 0
fi

if [[ $# -eq 1 ]] ; then
    echo 'Vous avez entré un argument, ce script en demande deux. ($1=nom de gène valide	$2=fichier_newick)'
    exit 0
fi

if [[ -d "$new_directory" ]];
then
    echo "$new_directory  this directory exists."
    echo ""
    gene_check=1
    FILE_NEWICK= $newick_file
    if [[ -f "$newick_file" ]]; then
        echo "$newick_file     this Newick file exists."
        echo ""
        newick_check=1
    else 
        echo "$FILE_NEWICK There is no newick file in the given path."
    fi    
else
	echo "$DIR directory does not exist. Maybe an error il 'script_launcher' run ?"
fi




if [[ $gene_check -eq 1 ]] && [[ $newick_check -eq 1 ]] ; then

	cp ../scripts/script_Dataframe_creator.py $new_directory/RNAduplex_aligned/.
	cp ../scripts/script_filter_dataframe.py $new_directory/RNAduplex_aligned/.
	cp ../scripts/script_keep_only_duplicate_species.py $new_directory/RNAduplex_aligned/.
	cp ../scripts/script_for_scripts.py $new_directory/RNAduplex_aligned/.

	touch $new_directory/RNAduplex_aligned/Dataframe_file
	touch $new_directory/RNAduplex_aligned/filtered_Dataframe_file
	echo "Starting script_Dataframe_creator.py and script_filter_dataframe"
	echo ""
	cp $newick_file $new_directory/RNAduplex_aligned/$newick_file
	python3 $new_directory/RNAduplex_aligned/script_Dataframe_creator.py $new_directory/RNAduplex_aligned/merged_RNAduplex $new_directory/RNAduplex_aligned/$newick_file
	python3 $new_directory/RNAduplex_aligned/script_filter_dataframe.py Dataframe_file filtered_Dataframe_file
	python3 $new_directory/RNAduplex_aligned/script_keep_only_duplicate_species.py filtered_Dataframe_file
	mv filtered_Dataframe_file $new_directory/RNAduplex_aligned/
	echo "Finished script_Dataframe_creator.py and script_filter_dataframe"
	rm Dataframe_file
	python3 $new_directory/RNAduplex_aligned/script_for_scripts.py $new_directory/RNAduplex_aligned/filtered_Dataframe_file $new_directory/RNAduplex_aligned/script_finalHeatmap.py
	echo ""
	echo "Generating heatmap files !!"
	python3 $new_directory/RNAduplex_aligned/script_finalHeatmap.py $new_directory/RNAduplex_aligned/merged_RNAduplex $new_directory/RNAduplex_aligned/$newick_file
	echo "Finished all tasks sucesfully"
	echo ""
	echo ""
	mkdir output
	mkdir output/$gene_name'_outputs'
	mkdir output/$gene_name'_outputs'/Heatmaps_and_Trees

	mv heatmap_miRNA_matches.svg output/$gene_name'_outputs'/Heatmaps_and_Trees/.
	mv heatmap_DNA_matches.svg output/$gene_name'_outputs'/Heatmaps_and_Trees/.
	echo " Heatmap have been generated, and stored in output/'gene'_outputs/"
	
fi	
