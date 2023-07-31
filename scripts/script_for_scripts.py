import argparse

def extract_species_indexes(tsv_file):
    species_indexes = {}
    with open(tsv_file, 'r') as f:
        for line in f:
            species_name, index = line.strip().split('\t')
            species_indexes[species_name] = index
    return species_indexes

def write_commands(species_indexes, text_file):
    with open(text_file, 'r') as f:
        lines = f.readlines()

    with open(text_file, 'w') as f:
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("# This line is purely esthetic, so that we put out specie of interest"):
                f.write(f"{lines[i]}\n")  # Write the description line
                for species, index in species_indexes.items():
                    f.write(f'merged_df.loc[merged_df["species"] == "{species}", "index"] = {index}\n')
            else:
                f.write(f"{lines[i]}")  # Write the line as it is
            i += 1

    print("Commands written to the text file successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract species and indexes from TSV file and write commands.")
    parser.add_argument("tsv_file", help="Path to the TSV file")
    parser.add_argument("text_file", help="Path to the text file")
    args = parser.parse_args()

    species_indexes = extract_species_indexes(args.tsv_file)
    write_commands(species_indexes, args.text_file)
