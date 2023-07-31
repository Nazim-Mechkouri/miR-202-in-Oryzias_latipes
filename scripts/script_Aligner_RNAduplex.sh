#!/bin/bash

dire=$1

mkdir $dire/RNAduplex_aligned

for file in $dire/*; do
    #echo $file
    if [[ -f $file ]]; then
        p=$file
        base="${file%.*}" 
        base=$( echo ${base##*/} )
         
          # Extract the filename without extension
        extension="${file##*.}"  # Extract the file extension
	    #echo $base
	    #echo $extension
        # Create the new file


        new_file="${base}_Aligned.${extension}"  # Create the new filename

        echo -n "" > $dire/RNAduplex_aligned/$new_file

     	RNAduplex < $file > $dire/RNAduplex_aligned/$new_file
        
        #echo "Created file: $new_file"

    fi
done
