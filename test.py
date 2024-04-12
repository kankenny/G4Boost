from Bio import SeqIO

input_file = "test_sequence.fasta"

fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    print(f'Name: {name}, Sequence: {sequence}\n')
