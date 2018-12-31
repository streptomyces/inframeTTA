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

# inframeCodons.pl

### Examples

     perl inframeCodons.pl -help

     perl inframeCodons.pl -codons tta,ttt -out sco_plasmids_codons.txt -fnafn sco_
     -faafn sco_plasmids_codons.faa sco_plasmids.gbk 

     perl inframeCodons.pl -codons tta,ttt -out sco_plasmids_codons.txt -ntfas sco_
     -protfas sco_plasmids_codons.faa sco_plasmids.gbk 

     perl inframeCodons.pl -codons tta,ttt -ntfas sco_plasmids_codons.fna \
     -protfas sco_plasmids_codons.faa sco_plasmids.gbk 

     perl inframeCodons.pl -codons tta,ttt,ctt -outfile sco_plasmids_codons.txt \
     -ntfas sco_plasmids_codons.fna sco_plasmids.gbk

### Options

* -outfile: File name to which tabular (tab delimited) output is written.

* -fnafn|-ntfas: File name to which the nucleotide sequences of the genes
  containing the codons are written in the fasta format.

* -faafn|-protfas: File name to which the protein sequences of the genes
  containing the codons are written in the fasta format.

* -codons: Codons to look for and report. Comma separated. No spaces in
  the argument. For example `-codons tta,ttt,ctt`. If no codons are specified
  then TTA is searched for.

* Any remaining arguments after the above options are considered to be
  input filenames to be processed.

### Notes

The number of columns in the tabular output depends upon the number of
codons searched.

If the outfile is not specified tabular output is written to STDOUT
(terminal).

If the nucleotides output file is not specified the nucleotide sequences
of the TTA containing genes found are not saved anywhere.

If the proteins output file is not specified the protein sequences of the
TTA containing genes found are not saved anywhere.

# pretty_translate.pl

### Examples

     perl pretty_translate.pl sco_plasmids_codons.fna
 
     perl pretty_translate.pl -outfile sco_plasmids_codons_pretty.txt \
     sco_plasmids_codons.fna
 
     perl pretty_translate.pl -frame 2 sco_plasmids_codons.fna
 
     perl pretty_translate.pl -strand -1 -frame 2 sco_plasmids_codons.fna

# Options

* -strand: 1 or -1. -1 will revcom before translating. Defaults to 1

* -frame: 0 or 1 or 2. Reading frame to translate. Defaults to 0.

* -outfile: If specified, output is written to this file. Otherwise to
  STDOUT.

* Remaining non-option arguments are considered to be fasta file names to
  be translated.



### Dependencies

These scripts depend on [BioPerl](https://bioperl.org) being installed
and available. They only need Bio::SeqIO (and whatever it depends upon)
to work.

### Author

govind.chandra@jic.ac.uk

