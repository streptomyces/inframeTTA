# inframeTTA.pl

### Examples

~~~ {.sh}

 perl code/inframeTTA.pl -out sco_plasmids_tta.txt -fnafn sco_plasmids_tta.fna \
 -faafn sco_plasmids_tta.faa sco_plasmids.gbk 

 perl code/inframeTTA.pl -out sco_plasmids_tta.txt -ntfas sco_plasmids_tta.fna \
 -protfas sco_plasmids_tta.faa sco_plasmids.gbk 

 perl code/inframeTTA.pl -ntfas sco_plasmids_tta.fna \
 -protfas sco_plasmids_tta.faa sco_plasmids.gbk 

 perl code/inframeTTA.pl -outfile sco_plasmids_tta.txt \
 -ntfas sco_plasmids_tta.fna sco_plasmids.gbk 

~~~

### Options


* -outfile: File name to which tabular (tab delimited) output is written.


* -fnafn|-ntfas: File name to which the nucleotide sequences of the TTA
containing genes are written in the fasta format.

* -faafn|-protfas: File name to which the protein sequences of the TTA containing
genes are written in the fasta format.

* Any remaining arguments after the above options are considered to be
input filenames to be processed.

### Notes

If the outfile is not specified tabular output is written to STDOUT (terminal).

If the nucleotides output file is not specified the nucleotide sequences of the
TTA containing genes found are not saved anywhere.

If the proteins output file is not specified the protein sequences of the
TTA containing genes found are not saved anywhere.

