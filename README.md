# Sequence-Data-Preprocesser-
This program is a preprocesser that removes primers and adapters from sequences in a fasta file. The assumption is made that the adapters are at the ends of the sequences, attached to the primers. Biopython is used to parse the file.

## **Assumed Sequence Makeup:** <br/>
[adapters][forward primer] - [sequence] - [reverse primer][adapter]<br/>

## **Arguments:**<br/>
--in_file - Fasta file with sequences to be trimmed<br/>
--out_file - Name of trimmed sequence file<br/>
--unk_file - Name of file with sequences that did not meet requirements<br/>
--n_mismatch - Tolerance for mismatches<br/>
--min_len The minimum length of sequence to be trimmed<br/>
--forward Forward primer<br/>
--reverse Reverse primer<br/>

## **Required Programs**
Biopython is used to parse the fasta file. It can be downloaded following the instructions on https://biopython.org/wiki/Download

### **Author**
Jennien Raffington

