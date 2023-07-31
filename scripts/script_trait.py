import pandas as pd
import sys


def read_species_file(file_path):
    species_data = pd.read_csv(file_path, sep='\t', header=None, names=['Species', 'Data'])
    return species_data

def compare_species(species_data):
    species_data['Species'] = species_data['Species'].str.replace(r'_(\d+)$', '', regex=True)

    # Remove duplicates based on both 'Species' and 'Data' columns
    species_data.drop_duplicates(subset=['Species', 'Data'], keep='first', inplace=True)

    # Fill NaN values in 'Data' column with 'null'
    species_data['Data'].fillna('null', inplace=True)

if __name__ == "__main__":
    file_path = sys.argv[1]
    species_data = read_species_file(file_path)
    compare_species(species_data)

    # Write the output to a file without column headers
    output_file = "merged_RNAduplex"
    species_data.to_csv(output_file, sep='\t', header=False, index=False, mode='w')
