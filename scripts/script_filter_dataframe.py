import argparse

def fill_null_indexes(species_data):
    species_map = {}
    result = []

    # First, create a mapping of species names to their indexes
    for line in species_data:
        name, index = line.strip().split('\t')
        if "_2" not in name:
            species_map[name] = index

    # Then, process the data again to fill in the "null" indexes
    for line in species_data:
        name, index = line.strip().split('\t')
        if "_2" in name and index == "null":
            base_species = name.replace("_2", "")
            if base_species in species_map:
                index = species_map[base_species]
        result.append(f"{name}\t{index}")

    return result


def main(input_file, output_file):
    try:
        with open(input_file, 'r') as f:
            tsv_data = f.readlines()

        # Fill the "null" indexes
        updated_lines = fill_null_indexes(tsv_data)

        # Write the updated TSV data to the output file and remove "null" lines
        with open(output_file, 'w') as f:
            for line in updated_lines:
                name, index = line.strip().split('\t')
                if index != "null":
                    f.write(f"{name}\t{index}\n")
        print(f"Data written to {output_file} successfully.")
    except FileNotFoundError:
        print("Input file not found.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fill null indexes in TSV file and remove 'null' lines.")
    parser.add_argument("input_file", help="Path to the input TSV file")
    parser.add_argument("output_file", help="Path to the output TSV file")
    args = parser.parse_args()

    main(args.input_file, args.output_file)

