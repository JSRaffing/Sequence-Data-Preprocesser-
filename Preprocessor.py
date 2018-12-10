#!/usr/bin/python

import argparse
# all the options for this script
parser = argparse.ArgumentParser(description='practicing command line options.')
parser.add_argument('--in_file', action="store", type=str, required=True, help='Name of the inputted fasta file')
parser.add_argument('--out_file', action="store", type=str, required=True, help='Name of the trimmed reads file')
parser.add_argument('--unk_file', action="store", type=str, required=True, help='Name of the file with unprocessed reads')
parser.add_argument('--n_mismatch', nargs='?', type=int, default=0, action="store", help='The tolerance of mismatches')
parser.add_argument('--min_len', nargs='?', type=int, default=0, action="store", help='The minimum length of sequence')
parser.add_argument('--forward', action="store", type=str, required=True, help='The forward primer')
parser.add_argument('--reverse', action="store", type=str, required=True, help='The reverse primer')

args = parser.parse_args()

# creating variables of inputted options
infile = args.in_file
outfile = args.out_file
unkfile = args.unk_file
n_mismatch = args.n_mismatch
minlen = args.min_len
forwardp = args.forward
reversep = args.reverse


d = {} # creating dictionary for parsed fasta file
reasonswhy = {} # creating dictionary for reasons why sequence isn't parsed
from Bio import SeqIO
# parsing fasta file using BioPython
for seq_record in SeqIO.parse(infile, "fasta"):
    d[str(seq_record.id)] = str(seq_record.seq)

print(d)
Length_dict = len(d)  # counts number of processed sequences

# dictionary for what the sequence should return when comparing to primer
complement = {"A": "A", "G": "G", "C":"C", "T": "T", "M": ("A", "C", "M"), "N": ("A", "T", "G", "C")}

x = {}  # opening dictionary for forward direction slices
l = {}  # opening dictionary for sequences that pass min length requirement
j = {}  # opening dictionary for sequences that do not pass min length

# cutting into slices, by iterating over index and splitting by length of primer
for k in d:

    if len(d[k]) >= minlen:
        x[k] = [str(d[k][i:i+(len(forwardp))]) for i in range(len(forwardp))]
        l[k] = str(d[k])
    else:
        j[k] = [str(d[k])]
        reasonswhy[k] = "length is less than 150"  # appends to dictionary with reason
        b = open(unkfile, "a")  # opening file for untrimmed sequences
        b.write(">" + k + "\n" + d[k] + "\n")
print(x)
print(l)

z = {}  # dictionary that for each sequence has list of differences between slice and primer
# count mismatches between forward primer and each string in slice by
for k in x:
        z[k] = [sum(ch1 not in complement.get(ch2) for ch1, ch2 in zip(string, forwardp)) -1 for string in x[k]]
print(z)

# find lowest number of mismatches in list, and index of slice with lowest number of mismatches
h = {}  # dictionary with lowest number of mismatches
g = {}  # dictionary with index of slice that has lowest number of mismatches
for k, v in z.items():
    h[k] = min(v)
    g[k] = v.index(min(v))
print(h)
print(g)

# cut section out of string based on index of lowest mismatch

p = {}  # dictionary for sequences that meet limit of mismatches
n = {}  # dictionary for sequences that are over the limit of mismatches
for k in l:
    if h[k] <= n_mismatch:
        p[k] = l[k][g[k]+len(forwardp)::]
    else:
        n[k] = str(l[k])
        reasonswhy[k] = "more than 3 mismatches for forward primer"  # appends to dictionary with reasons
        b = open(unkfile, "a")  # writes to file of untrimmed sequences
        b.write(">" + k + "\n" + l[k] + "\n")
print(p)
print(n)

# section for trimming reverse primer now
# dictionary for complement nucleotides
romplement = {"A": "T", "G": "c", "C": "G", "T": "A", "M": ("A", "C"), "N": ("A", "T", "G", "C")}

x1 = {}  # dictionary of sequence reverse slices
j1 = {}  # dictionary
l1 = {}  # dictionary of sequences that pass through the reverse slicing

# cutting sequences into slices in reverse
for k in p:
    x1[k] = [p[k][::-1][i:i + (len(reversep))] for i in range(len(reversep))]
    l1[k] = p[k][::-1]
print(x1)
print(l1)

# count mismatches between reverse primer and complement of reverse slices
z1 = {}  # dictionary for number of mismatches of reverse primer and and complement of reverse slices
for k in x1:
        z1[k] = [sum(romplement[ch1] != ch2 for ch1, ch2 in zip(reversep, string)) - 1 for string in x1[k]]
print(z1)

# find slice with lowest number of mismatches for each sequence
h1 = {}  # dictionary of lowest number of mismatches for each sequence
g1 = {}  # dictionary of index of slice with lowest number of mismatches
for k in z1:
    h1[k] = min(z1[k])
    g1[k] = z1[k].index(min(z1[k]))
print(h1)
print(g1)

# Trimming slice out of the string while sending sequences that don't have a minimum to another dictionary and file
p1 = {}
n1 = {}
for k in l1:
        if h1[k] <= n_mismatch:
            p1[k] = l1[k][g1[k]+len(reversep)::][::-1]
            c = open(outfile, "a")  # send trimmed sequences to file
            c.write(">" + k + "\n" + p1[k] + "\n")
        else:
            n1[k] = str(l1[k])
            reasonswhy[k] = "more than 3 mismatches for reverse primer"
            b = open(unkfile, "a")  # file that contains sequences that don't meet mismatch limit
            b.write(">" + k + "\n" + l1[k] + "\n")
print(p1)
print(n1)
print("Almost done......")
trimmed_sequences = len(p1)  # number of reads trimmed
avg_trimmed_sequences = sum([len(x) for x in p1.values()])/trimmed_sequences  # average length of trimmed reads
avg_raw_sequences = sum([len(x) for x in d.values()])/Length_dict  # average length of raw reads
print(trimmed_sequences)
print(Length_dict)
print(avg_trimmed_sequences)
print(avg_raw_sequences)

jr = open("log.txt", "a")  # file with basic statistics
jr.write('Processed Reads: ' + repr(Length_dict) + '\n')
jr.write('Trimmed Reads: ' + repr(trimmed_sequences) + '\n')
jr.write('Average Read Length (Raw): ' + repr(avg_raw_sequences) + '\n')
jr.write('Average Read Length (Trimmed): ' + repr(avg_trimmed_sequences) + '\n' )
for k in reasonswhy:  # appends reasons why sequences were not processed
    jr.write(">" + k + '\n' + reasonswhy[k] + '\n')
