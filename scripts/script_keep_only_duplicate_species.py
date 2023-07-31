import argparse

def filter_species_by_extension(species_data):
    result = []
    extensions_to_keep = ["_2", "_3", "_4", "_5"]

    for line in species_data:
        name, index = line.strip().split('\t')
        if any(extension in name for extension in extensions_to_keep):
            result.append(f"{name}\t{index}")

    return result


def main(input_file):
    try:
        with open(input_file, 'r') as f:
            tsv_data = f.readlines()

        # Filter species by extension
        filtered_lines = filter_species_by_extension(tsv_data)

        # Write the filtered TSV data back to the same file
        with open(input_file, 'w') as f:
            f.write('\n'.join(filtered_lines))
        print(f"Filtered data written back to {input_file} successfully.")
    except FileNotFoundError:
        print("Input file not found.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter species by extension and update the file.")
    parser.add_argument("input_file", help="Path to the input TSV file")
    args = parser.parse_args()

    main(args.input_file)
