import os
import sys

output_file = "unfiltered_merged_RNAduplex"
directory = sys.argv[1]  # Replace with the directory path where your files are located

merged_data = []

for filename in os.listdir(directory):
    if filename.endswith(".fasta"):  # Assuming all files have the .txt extension
        print("I found a file")
        file_path = os.path.join(directory, filename)
        with open(file_path, "r") as file:
            print("I am reading the file")
            lines = file.readlines()
            
              # Check if the file has the required number of lines
            seq_id = lines[0].strip()
            print(seq_id)
            specie_name = lines[1].strip()
            print(specie_name)
            try:
              data = lines[2].strip()
            except IndexError:
              data = 'null'
            print(data)
            merged_data.append(f"{specie_name}\t{data}")

# Write merged data to the output file
with open(output_file, "w") as outfile:
    outfile.write("\n".join(merged_data))


