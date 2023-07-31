import sys

def update_species_names(species_dict, species):
    if species in species_dict:
        species_dict[species] += 1
        new_species_name = f"{species}_{species_dict[species]}"
    else:
        species_dict[species] = 1
        new_species_name = species
    return new_species_name


def process_file(input_file):
    species_dict = {}
    with open(input_file, "r+") as infile:
        lines = infile.readlines()
        infile.seek(0)  # Move the file pointer to the beginning of the file
        for line in lines:
            line = line.strip().replace(">", "")  # Delete ">" from the line
            if line.startswith("#"):  # Skip any commented lines
                infile.write(f"{line}\n")
                continue
            
            species_name, rest = line.split("\t", 1)
            new_species_name = update_species_names(species_dict, species_name)
            
            infile.write(f"{new_species_name}\t{rest}\n")
        infile.truncate()  # Remove any remaining content after the last line


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input_file")
        sys.exit(1)

    input_file_path = sys.argv[1]
    process_file(input_file_path)

