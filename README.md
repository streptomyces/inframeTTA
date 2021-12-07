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

### Description

Look for codons in CDSs in genbank files.

### Examples

    perl inframeCodons.pl -help

    perl inframeCodons.pl -codons tta,ttt -out sco_plasmids.txt \
    -fnafn sco_plasmids.fna -faafn sco_plasmids.faa sco_plasmids.gbk 

    perl inframeCodons.pl -codons tta,ttt -out sco_plasmids.txt \
    -ntfas sco_plasmids.fna -protfas sco_plasmids.faa sco_plasmids.gbk 

    perl inframeCodons.pl -codons tta,ttt -ntfas sco_plasmids.fna \
    -protfas sco_plasmids.faa sco_plasmids.gbk 

    perl inframeCodons.pl -codons tta,ttt,ctt -outfile sco_plasmids.txt \
    -ntfas sco_plasmids.fna sco_plasmids.gbk 

    perl inframeCodons.pl -codons tta,ttt -out outdir/sco_plasmids.txt \
    -fnafn outdir/sco_plasmids.fna -faafn outdir/sco_plasmids.faa \
    sco_plasmids.gbk

### Notes

The number of columns in the tabular output written to STDOUT depends
upon the number of codons searched. This information is also in the per
codon tabular output files so this output can be safely ignored. It is
just to give you something to watch while the script is working.

If the nucleotides output file is not specified the nucleotide sequences
of the codon containing genes found are not saved anywhere.

If the proteins output file is not specified the protein sequences of
the codon containing genes found are not saved anywhere.

#### Output file names

Tabular, protein and, nucleotide fasta output file names are generated
using the output filenames provided and the list of codons to search.
For example, if the options used are

    -codons ttt,tta,ctt -fnafn directory/out.fna -outfile directory/out.txt

then the three file names generated are

    directory/out_TTT.fna
    directory/out_TTA.fna
    directory/out_CTT.fna
    directory/out_TTT.txt
    directory/out_TTA.txt
    directory/out_CTT.txt

If the directory path in the output file name does not already exist
then an attempt will be made to make it. Failure of this attempt will
stop the script.

### Options

* -outfile: File name to which tabular (tab delimited) output is
written. See *output file names* above.

* -fnafn|-ntfas: File name to which the nucleotide sequences of the
genes containing the codons are written in the fasta format.
See *output file names* above.

* -faafn|-prottfas: File name to which the protein sequences of the
genes containing the codons are written in the fasta format.
See *output file* names above.

* -codons: Codons to look for and report. Comma separated. No spaces in
the argument. e.g. `-codons tta,ttt,ctt`.
If no codons are specified then TTA is searched for. All specified
codons are converted to upper case without any warning. All
occurrences of *U* in codons are converted to *T* without warning.
Codons containing anything other than *A,C,G,T,U* are dropped as are
codons not exactly three characters long.

* Any remaining arguments after the above options are considered to be
input filenames to be processed.

# pretty_translate.pl

### Description

Translates the nucleotide sequences in the given fasta files and the
prints out both the nucleotide and amino acid sequences so that it is easy
to see which amino acid is coming from which codon. See below for an
example output.

    >SCP1.11c tta - ttt 15,288,318 hypothetical protein
    TTGATTCCCTCGGGCTTTTCCCTCAAGGGCGGTGTCAGTCCGATACGTCAGCATGACCGT
    L  I  P  S  G  F  S  L  K  G  G  V  S  P  I  R  Q  H  D  R

    CGGGAACGACATGCCGCCCTGGAGGTCGCCATGCGCAAGCCACACGATCCTCACGAGCAG
    R  E  R  H  A  A  L  E  V  A  M  R  K  P  H  D  P  H  E  Q

### Examples

     perl pretty_translate.pl sco_plasmids_TTA.fna
 
     perl pretty_translate.pl -outfile sco_plasmids_TTA_pretty.txt \
     sco_plasmids_TTA.fna
 
     perl pretty_translate.pl -frame 2 sco_plasmids_TTA.fna
 
     perl pretty_translate.pl -strand -1 -frame 2 sco_plasmids_TTA.fna

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

Govind Chandra (govind.chandra@jic.ac.uk)

