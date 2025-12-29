#! /usr/bin/env python3

import argparse
import Bio.SeqIO

parser = argparse.ArgumentParser(description="Simple Prokaryotic GFF Annotation Tool")

parser.add_argument("--input","-i",required=True, help="Input Fasta File")
parser.add_argument("--min_len","-m", type=int, default=100, help="Minimum ORF length")
parser.add_argument("--translate","-t", type=int, default=0, help="Number of ORFs to translate")
parser.add_argument("--out","-o", default="output", help="Output file prefix")

args = parser.parse_args()

print("input:{0},min_len:{1},translate:{2},output:{3}".format(args.input,args.min_len,args.translate,args.out))

# Lists to store gff lines, GC results and translated protein sequences
gff_lines = []
gc_results = []
protein_results = []

# Count for ORFs
orf_count = 1

# Function to make the reverse complement of the sequence
def rev_comp(seq):
    comp = {"A":"T","T":"A","G":"C","C":"G"}
    reverse = ""
    for base in seq[::-1]:
        reverse += comp[base]
    return reverse

# Function to show codons
def show_codons(seq,frame):
    for pos in range(frame, len(seq) -2, 3):
        codon = seq[pos:pos+3]
        print("pos:{0},codon:{1}".format(pos,codon))

# Function to translate codons into protein
def start_codon(seq,frame):
    for pos in range(frame, len(seq)-2, 3):
        codon = seq[pos:pos+3]
        if codon == "ATG":
            print("Start at postion:{0},codon:{1}".format(pos,codon))

# Function to calculate GC content
def gc_content(seq):
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100

def translate_orf(seq):
    codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein = ""
    # Read codons in steps of 3 bases
    for idx in range(0,len(seq),3):
        codon = seq[idx:idx+3]
        # Translate full codons only
        if len(codon) == 3:
            protein += codon_table.get(codon,"X")
    return protein

# Function to find ORFs
def get_orfs(seq,frame,strand,min_len,record_id):
    global orf_count
    stop_codons = {"TAA","TAG","TGA"}

    # Go through the sequence in steps of 3
    for start_pos in range(frame, len(seq)-2,3):
        start_codon = seq[start_pos:start_pos+3]

        # Look for the start codon ATG
        if start_codon == "ATG":

            # After ATG, look forward for a stop codon
            for end_pos in range(start_pos + 3, len(seq)-2, 3):
                codon = seq[end_pos:end_pos + 3]
                if codon in stop_codons:
                    orf_len = end_pos + 3 - start_pos

                    # Check for the length requirement
                    if orf_len > min_len:
                        print("ORF positions:{0} to {1}".format(start_pos,end_pos))
                        print("Codons:{0} to {1}".format(start_codon,codon))

                        # Take the ORF sequence and print
                        orf_seq = seq[start_pos:end_pos+3]
                        print("ORF SEQUENCE:",orf_seq)

                        # Calculate GC content
                        gc_value = gc_content(orf_seq)
                        print("GC content:{0}%".format(gc_value))
                        gc_results.append("orf{0}:{1}%".format(orf_count,gc_value))

                        # Translate if -t is given
                        if args.translate > 0 and orf_count <= args.translate:
                            protein = translate_orf(orf_seq)
                            print("Protein:",protein)
                            protein_results.append(">orf{0}\n{1}".format(orf_count,protein))


                        # Create GFF annotation line
                        gff_line = "{0}\tORFfinder\tCDS\t{1}\t{2}\t.\t{3}\t0\tID=orf{4}".format(record_id, start_pos+1, end_pos+3, strand, orf_count)
                        gff_lines.append(gff_line)
                        orf_count += 1
                    break

# Read input fasta file
for record in Bio.SeqIO.parse(args.input,"fasta"):
    print("ID:",record.id)
    print("Seq:",record.seq)
    print("Length:",len(record.seq))

    seq = str(record.seq)
    rev_seq = rev_comp(seq)
    record_id = record.id
    break

# Find ORFs on + strand of frame 0,1,2
for frame in range(3):
    get_orfs(seq,frame,"+",args.min_len,record_id)

# Find ORFs on - strand of frames 0,1,2
for frame in range(3):
    get_orfs(rev_seq,frame,"-",args.min_len,record_id)
prefix = args.out
gff_file = prefix + ".gff"
gc_file = prefix + "_gc.txt"
protein_file = prefix + "_protein.txt"

# Write output GFF file
with open(gff_file,"w") as outputgff:
    outputgff.write("##gff-version 3\n")
    for line in gff_lines:
        outputgff.write(line + "\n")
print("GFF file written:", gff_file)


# Write GC content file
with open(gc_file, "w") as gcfile:
    for line in gc_results:
        gcfile.write(line + "\n")
print("GC content written:",gc_file)

# Write protein translations
with open(protein_file, "w") as pfile:
    for line in protein_results:
        pfile.write(line + "\n")
print("Protein sequence written:", protein_file)
