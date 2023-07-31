import sys
import os

candidates = sys.argv[1]
orthologous = sys.argv[2]

temp = sys.argv[1].split("/")[-1]
tmp = temp.split(".")[0]
name = tmp.replace("miR-202_in_", "")
print(name)

gene_name = sys.argv[2].split("/")[-1].split("_")[1]
print(gene_name)

# Create the output directory
output_dir = gene_name + '_candidates_combined_orthologous'
os.makedirs(output_dir, exist_ok=True)


def read_fasta(filename):
    sequences = {}
    with open(filename, 'r') as file:
        seq_id = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                seq_id = line[1:]
                sequences[seq_id] = []
            else:
                sequences[seq_id].append(line.replace('-', ''))
    return sequences


def echo_sequences(fasta1, fasta2):
    sequences1 = read_fasta(fasta1)
    sequences2 = read_fasta(fasta2)


    cpt = 0
    for i, (seq_id, seq) in enumerate(sequences1.items(), start=1):
        if i > 4:
            break # Keep only the first 4 sequences


        for seq_id_2, seq_2 in sequences2.items():
            if seq_id_2 == name:
                cpt+=1
                print(cpt)
                output_filename_2 = f'{output_dir}/{name}_combined_{i}_orthologous.fasta'
                with open(output_filename_2, 'w') as output_file_2:
                    output_file_2.write(f'>{seq_id}\n')
                    output_file_2.write('\n'.join(seq) + '\n')
                    output_file_2.write(f'>{seq_id_2}\n')
                    output_file_2.write('\n'.join(seq_2) + '\n')


echo_sequences(candidates, orthologous)