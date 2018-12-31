# use 5.14.0;
use strict;
use Bio::SeqIO;
use Carp;
use File::Basename;
use File::Spec;
use File::Path qw(make_path remove_tree);

# {{{ Getopt::Long
use Getopt::Long;
my $outfile;
my $fnafn;
my $faafn;
my $help;
my $glcodons = qq(TTA);
GetOptions (
"outfile:s" => \$outfile,
"fnafn|ntfas:s" => \$fnafn,
"faafn|protfas:s" => \$faafn,
"codons:s" => \$glcodons,
"help" => \$help
);
# }}}

if($help) {
exec("perldoc $0");
exit;
}

# {{{ POD

=head1 Name

inframeCodons.pl

=head2 Description

Look for codons in CDSs in genbank files.

=head2 Examples

 perl inframeCodons.pl -help

 perl inframeCodons.pl -codons tta,ttt -out sco_plasmids_codons.txt \
 -fnafn sco_plasmids_codons.fna -faafn sco_plasmids_codons.faa sco_plasmids.gbk 

 perl inframeCodons.pl -codons tta,ttt -out sco_plasmids_codons.txt \
 -ntfas sco_plasmids_codons.fna -protfas sco_plasmids_codons.faa sco_plasmids.gbk 

 perl inframeCodons.pl -codons tta,ttt -ntfas sco_plasmids_codons.fna \
 -protfas sco_plasmids_codons.faa sco_plasmids.gbk 

 perl inframeCodons.pl -codons tta,ttt,ctt -outfile sco_plasmids_codons.txt \
 -ntfas sco_plasmids_codons.fna sco_plasmids.gbk 


=head2 Notes

The number of columns in the tabular output depends upon the number of codons
searched.

If the outfile is not specified tabular output is written to STDOUT (terminal).

If the nucleotides output file is not specified the nucleotide sequences of the
TTA containing genes found are not saved anywhere.

If the proteins output file is not specified the protein sequences of the
TTA containing genes found are not saved anywhere.

=head3 Output file names

Protein and nucleotide fasta output file names are generated using the output
filenames provided and the list of codons to search. For example, if the
options used are

 -codons ttt,tta,ctt -fnafn directory/out.fna

then the three file names generated are

 directory/out_TTT.fna
 directory/out_TTA.fna
 directory/out_CTT.fna

If the directory path in the output file name does not already exist then an
attempt will be made to make it. Failure of this attempt will stop the script.

=head2 Options

=over 2

=item

-outfile: File name to which tabular (tab delimited) output is written.

=item

-fnafn|-ntfas: File name to which the nucleotide sequences of the genes
containing the codons are written in the fasta format.

=item

-faafn|-prottfas: File name to which the protein sequences of the genes
containing the codons are written in the fasta format.

=item

-codons: Codons to look for and report. Comma separated. No spaces
in the argument. See example below.

 -codons tta,ttt,ctt

If no codons are specified then TTA is searched for.


=item 

Any remaining arguments after the above options are considered to be input
filenames to be processed.

=back


=cut

# }}}

# Command line arguments not eaten up by Getopt::Long above are
# considered to be input file names.
my @infiles = @ARGV;
unless(@infiles) {
  croak("No input files provided. Nothing to do. Exiting");
}


# {{{ The output filehandle. Default to STDOUT.
my $ofh;
if($outfile) {
  my ($noex, $dir, $ext)= fileparse($outfile, qr/\.[^.]*/);
  if($dir) {
    unless(-d $dir) {
      make_path($dir); # croaks or carps anyway. Hence no checks.
    }
  }
  open($ofh, ">", $outfile);
}
else {
  open($ofh, ">&STDOUT");
}
select($ofh);
# }}}


# {{{ Build @codons
my @temp = split(/(,|\s)+/, $glcodons);
my @codons;
for my $cod (@temp) {
  if($cod =~ m/[^acgtu]/i or length($cod) != 3) {
  } else {
    $cod = uc($cod);
    $cod =~ s/U/T/g;
    push(@codons, $cod);
  }
}
# }}}

# {{{ Hashes of filehandles and Bio::SeqIO objects for each codon
# for nucleotide and protein fasta output.
# fasta output for nucleotide sequences
my %ntfh;
my %ntout;
if($fnafn) {
  my ($noex, $dir, $ext)= fileparse($fnafn, qr/\.[^.]*/);
  unless(-d $dir) {
    make_path($dir); # croaks or carps anyway. Hence no checks.
  }
  for my $cod (@codons) {
    my $codfnafn = File::Spec->catfile($dir, $noex . "_" . $cod . $ext);
    open($ntfh{$cod}, ">", $codfnafn); 
    $ntout{$cod} = Bio::SeqIO->new(-fh => $ntfh{$cod}, -format => "fasta");
  }
}

# fasta output for protein sequences
my %protfh;
my %protout;
if($faafn) {
  my ($noex, $dir, $ext)= fileparse($faafn, qr/\.[^.]*/);
  unless(-d $dir) {
    make_path($dir); # croaks or carps anyway. Hence no checks.
  }
  for my $cod (@codons) {
    my $codfaafn = File::Spec->catfile($dir, $noex . "_" . $cod . $ext);
    open($protfh{$cod}, ">", $codfaafn); 
    $protout{$cod} = Bio::SeqIO->new(-fh => $protfh{$cod}, -format => "fasta");
  }
}
# }}}

# {{{ Loop through all CDSs inside each sequence object inside each genbank
# file. A single genbank file can have multiple sequences in it.

for my $infile (@infiles) {
  my $seqio=Bio::SeqIO->new(-file => $infile);
  while(my $seqobj=$seqio->next_seq()) {
    my $cdsCnt = 0;
    for my $feature ($seqobj->all_SeqFeatures()) {
      if($feature->primary_tag() eq 'CDS') {
        $cdsCnt += 1;
        my $codon_start = 1;
        if($feature->has_tag('codon_start')) {
          ($codon_start) = $feature->get_tag_values('codon_start');
        }
        my $offset = 0;
        if($codon_start > 1) { $offset += $codon_start - 1; }
        my $subobj = $feature->spliced_seq(-nosort => 1);
        my %codpos = inframeCodons($subobj, $offset);
        my @keycod = keys(%codpos);
        if(@keycod) {
          my @outlist;
          my ($id, $product) = idAndProd($feature, $cdsCnt);
          push(@outlist, $id);
          for my $codon (@codons) {
            if(exists($codpos{$codon})) {
              push(@outlist, $codon, join(",", @{$codpos{$codon}}) );
            }
            else {
              push(@outlist, $codon, "-");
            }
          }
          push(@outlist, $product);
          print(join("\t", @outlist), "\n");
          my @outdesc = @outlist[1..$#outlist];

          if(%ntout) {
            my $featobj=$feature->spliced_seq(-nosort => 1);
            $featobj->display_id($id);
            for my $cod (@keycod) {
              $featobj->description("$cod at: " . join(" ", @{$codpos{$cod}}) . ". " . $product);
              $ntout{$cod}->write_seq($featobj);
            }

          }
          if(%protout) {
            my $aaobj = featTranslate($feature);
            $aaobj->display_id($id);
            for my $cod (@keycod) {
              $aaobj->description("$cod at: " . join(" ", @{$codpos{$cod}}) . ". " . $product);
              $protout{$cod}->write_seq($aaobj);
            }
          }
        
        }
      }
    }
  }
}
# }}}

exit;


# {{{ Multiple END blocks run in reverse order of definition.
END {
  if($ofh) {
    close($ofh);
  }
  for my $cod (keys %ntfh) {
    close($ntfh{$cod});
  }
  for my $cod (keys %protfh) {
    close($protfh{$cod});
  }
}

# }}}


### subroutines below ###
=head1 Subroutines

=cut

# {{{ sub idAndProd

=head2 sub idAndProd

=head3 Arguments

=over 2

=item 1.
A sequence feature object.

=item 2.
A number. This is used to generate an identifier if the feature does not
contain any tag which might be used as an identifier.

=back

=head3 Returns

Returns two strings.

=cut

sub idAndProd {
  my $feature = shift(@_);
  my $cdsCnt = shift(@_);
  my %anno;
  for my $tag (qw(gene note product locus_tag gene_synonym)) {
    if($feature->has_tag($tag)) {
      my @taglist = $feature->get_tag_values($tag);
      $anno{$tag} = \@taglist;
    }
  }
  my $id;
  if(exists($anno{locus_tag})) {
    $id = $anno{locus_tag}->[0];
  }
  elsif(exists($anno{gene_synonym})) {
    $id = $anno{gene_synonym}->[0];
  }
  elsif(exists($anno{gene})) {
    $id = $anno{gene}->[0];
  }
  else {
    $id = "CDS_" . sprintf("%04d", $cdsCnt);
  }
  my $pproduct;
  if(exists($anno{product})) { $pproduct = join(" ", @{$anno{product}}); }
  else { $pproduct = "No annotated product"; }
  return($id, $pproduct);
}
# }}}


# {{{ sub inframeCodons

=head2 sub inframeCodons

=head3 Arguments

=over 2

=item 1.
A nucleotide sequence string or a Bio::Seq object.

=item 2.
An offset (integer) to start looking at triplets from. Zero to start
looking from the first nucleotide.

=back

=head3 Returns

Returns a list of (zero based) nucleotide positions.

=cut

sub inframeCodons {
  my $inseq=shift(@_);
  my $offset = shift(@_);
  my %ret;
  my $seq;
  if(ref($inseq)) { $seq = lc($inseq->seq()); }
  else { $seq = lc($inseq); }
  my $pos = $offset;
  while(my $triplet = substr($seq, $pos, 3)) {
    for my $codon (@codons) {
    if($triplet =~ m/$codon/i) {
      push(@{$ret{$codon}}, $pos);
    }
    }
    $pos += 3;
  }
  return(%ret);
}
# }}}



# {{{ sub featTranslate

=head2 sub featTranslate

=head3 Arguments

=over 2

=item 1.
A sequence feature object.

=back

=head3 Returns

Returns a Bio::Seq (protein) object.

=cut

sub featTranslate {
  my $feature=shift(@_);
  my $codon_start=1;
  if($feature->has_tag('codon_start')) {
    ($codon_start) = $feature->get_tag_values('codon_start');
  }
  my $aaobj;
  my $frame = 0;
  if($codon_start > 1) { $frame += ($codon_start-1);}
  my $featobj = $feature->spliced_seq(-nosort => 1);
  $aaobj = $featobj->translate(-frame => $frame, -complete => 0);
  return($aaobj);
}
# }}}


=head1 Author

Govind Chandra E<lt>govind.chandra@jic.ac.ukE<gt>

=cut

