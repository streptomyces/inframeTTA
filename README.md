# inframeTTA.pl

### Description

A simple Perl script to find genes containing inframe TTA codons in
them. TTA recognizing tRNA plays a regulatory role in the life cycle
of (at least some) _Streptomyces_ species. This role was first
discovered in _Streptomyces coelicolor_. It is often of interest to
look for TTA containing genes in newly sequenced _Streptomyces_
genomes, and possibly in other high GC bacterial genomes.

The input to this script is a Genbank (or EMBL) file. It looks for TTA
codons in all features with the primary tag of CDS.

### Help

    perl inframeTTA.pl -help

will show some help. It calls the `perldoc` command.


### Examples

~~~ {.sh}

 perl inframeTTA.pl -help

 perl inframeTTA.pl -out sco_plasmids_tta.txt -fnafn sco_plasmids_tta.fna \
 -faafn sco_plasmids_tta.faa sco_plasmids.gbk 

 perl inframeTTA.pl -out sco_plasmids_tta.txt -ntfas sco_plasmids_tta.fna \
 -protfas sco_plasmids_tta.faa sco_plasmids.gbk 

 perl inframeTTA.pl -ntfas sco_plasmids_tta.fna \
 -protfas sco_plasmids_tta.faa sco_plasmids.gbk 

 perl inframeTTA.pl -outfile sco_plasmids_tta.txt \
 -ntfas sco_plasmids_tta.fna sco_plasmids.gbk 

~~~

### Options

* `-outfile`: File name to which tabular (tab delimited) output is written.


* `-fnafn` or `-ntfas`: File name to which the nucleotide sequences of the TTA
containing genes are written in the fasta format.

* `-faafn` or `-protfas`: File name to which the protein sequences of the TTA containing
genes are written in the fasta format.

* Any remaining arguments after the above options are considered to be
input filenames to be processed.

### Notes

If the outfile is not specified tabular output is written to STDOUT (terminal).

If the nucleotides output file is not specified the nucleotide sequences of the
TTA containing genes found are not saved anywhere.

If the proteins output file is not specified the protein sequences of the
TTA containing genes found are not saved anywhere.

### Dependencies

This script depends on [BioPerl](https://bioperl.org) being installed
and available. It only needs Bio::SeqIO (and whatever it depends upon)
to work.

### Author

govind.chandra@jic.ac.uk

