import sys, re, os

def remove_clade_names(newick_str):
    # Split the Newick string by commas and parentheses
    tokens = []
    current_token = ""
    for char in newick_str:
        if char in "(),":
            if current_token.strip():  # Ignore empty tokens (e.g., caused by double commas)
                tokens.append(current_token.strip())
            tokens.append(char)
            current_token = ""
        else:
            current_token += char
    if current_token.strip():
        tokens.append(current_token.strip())

    # Remove clade names and associated distances (e.g., "XXXX:0")
    cleaned_tokens = []
    for token in tokens:
        if ":" in token and token.split(":")[-1] == "0":
            continue  # Skip tokens with clade names and distance 0
        cleaned_tokens.append(token)

    # Reconstruct the Newick string without clade names
    cleaned_newick = "".join(cleaned_tokens)

    return cleaned_newick

# Sample Newick string

newick_file = sys.argv[1]


with open(newick_file, 'r') as tree_file:
    sample_newick = tree_file.read().strip()


cleaned_newick = remove_clade_names(sample_newick)
print(cleaned_newick)
