# use 5.14.0;
use strict;
use Bio::SeqIO;

my @infiles = @ARGV;

my $fnafn = "out.fna";
my $faafn = "out.faa";

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
            $ntout->write_seq($featobj);
          }
          if($aaout) {
            my $aaobj = featTranslate($feature);
            $aaobj->display_id($id);
            $aaout->write_seq($aaobj);
          }
        }
      }
    }
  }
}

exit;


### subroutines below ###

# {{{ sub idAndProd
# Returns two strings.
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
  else { $pproduct = "None"; }
  return($id, $pproduct);
}
# }}}


# {{{ sub inframeTTA
# Returns a list of (zero based) nucleotide positions.
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
# Returns a Bio::Seq
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
