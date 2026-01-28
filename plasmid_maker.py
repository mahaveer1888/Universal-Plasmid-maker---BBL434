import sys


def read_fasta(filename):
    seq = ""
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq



def read_design(filename):
    mcs_sequences = []
    antibiotic_sequences = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            dna_seq, name = [x.strip() for x in line.split(",")]

            if "cloning" in name.lower() or "mcs" in name.lower():
                mcs_sequences.append(dna_seq)
            else:
                antibiotic_sequences.append(dna_seq)

    return mcs_sequences, antibiotic_sequences



ORI_V = "TTGACATGCGTACGTTAGCTAGCTAGCGTACGTAGCTAGCTAGCTA"
REP_GENES = "ATGAAAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAA"
PRIMASE_HELICASE = "ATGGGCTACGATCGATCGTACGTAGCTAGCTAGCATCGATCG"



def build_plasmid(insert_seq, mcs_seqs, antibiotic_seqs):

    plasmid = ""

    # Default replication genes
    plasmid += ORI_V
    plasmid += REP_GENES
    plasmid += PRIMASE_HELICASE

    # Add all MCS sequences
    for mcs in mcs_seqs:
        plasmid += mcs

    # Insert gene
    plasmid += insert_seq

    # Add all antibiotic marker genes
    for ab in antibiotic_seqs:
        plasmid += ab

    return plasmid



def write_fasta(filename, header, sequence):
    with open(filename, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 70):
            f.write(sequence[i:i+70] + "\n")



if len(sys.argv) != 3:
    print("Usage: python plasmid_maker.py Input.Fa Design.txt")
    sys.exit(1)

input_fasta = sys.argv[1]
design_file = sys.argv[2]

insert_sequence = read_fasta(input_fasta)
mcs_list, antibiotic_list = read_design(design_file)

final_plasmid = build_plasmid(insert_sequence, mcs_list, antibiotic_list)

write_fasta("Output.Fa", "Universal_Plasmid", final_plasmid)

