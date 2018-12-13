# use 5.14.0;
use strict;
use Bio::SeqIO;
use Carp;

# {{{ Getopt::Long
use Getopt::Long;
my $outfile;
my $fnafn;
my $faafn;
my $help;
GetOptions (
"outfile:s" => \$outfile,
"fnafn|ntfas:s" => \$fnafn,
"faafn|protfas:s" => \$faafn,
"help" => \$help
);
# }}}

if($help) {
exec("perldoc $0");
exit;
}

# {{{ POD

=head1 Name

inframeTTA.pl

=head2 Examples

 perl inframeTTA.pl -help

 perl inframeTTA.pl -out sco_plasmids_tta.txt -fnafn sco_plasmids_tta.fna \
 -faafn sco_plasmids_tta.faa sco_plasmids.gbk 

 perl inframeTTA.pl -out sco_plasmids_tta.txt -ntfas sco_plasmids_tta.fna \
 -protfas sco_plasmids_tta.faa sco_plasmids.gbk 

 perl inframeTTA.pl -ntfas sco_plasmids_tta.fna \
 -protfas sco_plasmids_tta.faa sco_plasmids.gbk 

 perl inframeTTA.pl -outfile sco_plasmids_tta.txt \
 -ntfas sco_plasmids_tta.fna sco_plasmids.gbk 

=head2 Options

=over 2

=item

-outfile: File name to which tabular (tab delimited) output is written.

=item

-fnafn|-ntfas: File name to which the nucleotide sequences of the TTA
containing genes are written in the fasta format.

=item

-faafn|-protfas: File name to which the protein sequences of the TTA containing
genes are written in the fasta format.

=item 

Any remaining arguments after the above options are considered to be input
filenames to be processed.

=back

=head2 Notes

If the outfile is not specified tabular output is written to STDOUT (terminal).

If the nucleotides output file is not specified the nucleotide sequences of the
TTA containing genes found are not saved anywhere.

If the proteins output file is not specified the protein sequences of the
TTA containing genes found are not saved anywhere.

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
  open($ofh, ">", $outfile);
}
else {
  open($ofh, ">&STDOUT");
}
select($ofh);
# }}}


# fasta output for nucleotide sequences
my $ntout;
if($fnafn) {
  $ntout = Bio::SeqIO->new(-file => ">$fnafn", -format => "fasta");
}

# fasta output for protein sequences
my $aaout;
if($faafn) {
  $aaout = Bio::SeqIO->new(-file => ">$faafn", -format => "fasta");
}


# Loop through all CDSs inside each sequence object inside each genbank file.
# A single genbank file can have multiple sequences in it.

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
        my @ttapos = inframeTTA($subobj, $offset);
        if(@ttapos) {
          my ($id, $product) = idAndProd($feature, $cdsCnt);
          print(join("\t", $id, join(" ", @ttapos), $product), "\n");
          if($ntout) {
            my $featobj=$feature->spliced_seq(-nosort => 1);
            $featobj->display_id($id);
            $featobj->description($product);
            $ntout->write_seq($featobj);
          }
          if($aaout) {
            my $aaobj = featTranslate($feature);
            $aaobj->display_id($id);
            $aaobj->description($product);
            $aaout->write_seq($aaobj);
          }
        }
      }
    }
  }
}

# Multiple END blocks run in reverse order of definition.
END {
  if($ofh) {
    close($ofh);
  }
}
exit;


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


# {{{ sub inframeTTA

=head2 sub inframeTTA

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

sub inframeTTA {
  my $inseq=shift(@_);
  my $offset = shift(@_);
  my @retlist;
  my $seq;
  if(ref($inseq)) { $seq = lc($inseq->seq()); }
  else { $seq = lc($inseq); }
  my $pos = $offset;
  while(my $triplet = substr($seq, $pos, 3)) {
    if($triplet =~ m/[tu][tu]a/i) {
      push(@retlist, $pos);
    }
    $pos += 3;
  }
  return(@retlist);
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

