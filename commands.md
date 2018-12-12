# Natasha Soares

#### Tue 04 Dec 2018

Hi Natasha,

inframeTTA.pl is attached.

Run it as

  perl inframeTTA.pl sco_plasmids.gbk

I have attached sco_plasmids.gbk for testing.

You should get a tabular output on stdout as well as the files out.fna
and out.faa. These are fasta files containing the nucleotide and
protein sequences of the TTA containing genes.

In the tabular output

Column 1: is the gene identifier

Column 2: is the position of the TTA codon in the nucleotide sequence
of the gene. Note that this position is zero based.

Column 3: is the value of the "product" tag in the annotation for the gene.


As you can probably guess from the script itself, you may specify more
than one genbank files on the command line.

Let me know whether the script works for you or not. If it does not,
please paste any error messages you get in an email to me.


I have deliberately kept the script very simple i.e. haven't used
command line options processing or complex error reporting. Of course,
the complexity coming from using BioPerl (Bio::SeqIO) cannot be
avoided. The idea is to test a very minimal script before going any
further.

Best wishes.

Govind




~~~ {.sh}

zip -j -r ~/mnt/wstemp/natasha/ttaScriptAndTestInput code/inframeTTA.pl sco_plasmids.gbk

    git remote add origin git@github.com:streptomyces/inframeTTA.git
    git push -u origin master

~~~

