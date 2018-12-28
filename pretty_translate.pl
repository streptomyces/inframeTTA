# use 5.14.0
use Getopt::Long;
use Bio::SeqIO;
use strict;
use Carp;

my $strand = 1;
my $frame = 0;
my $help;
my $outfile;

GetOptions(
"strand:i" => \$strand,
"frame:i" => \$frame,
"outfile:s" => \$outfile,
"help" => \$help
);


if($help) {
exec("perldoc $0");
exit;
}

# {{{ POD

=head1 Name

pretty_translate.pl

=head2 Examples

 perl pretty_translate.pl sco_plasmids_codons.fna
 
 perl pretty_translate.pl -outfile sco_plasmids_codons_pretty.txt \
 sco_plasmids_codons.fna
 
 perl pretty_translate.pl -frame 2 sco_plasmids_codons.fna
 
 perl pretty_translate.pl -strand -1 -frame 2 sco_plasmids_codons.fna

=head2 Options

=over 2

=item

-strand: 1 or -1. -1 will revcom before translating. Defaults to 1

=item

-frame: 0 or 1 or 2. Reading frame to translate. Defaults to 0.

=item

-outfile: If specified, output is written to this file. Otherwise to STDOUT.

=item

Remaining non-option arguments are considered to be fasta file names to be
translated.

=back


=cut

# }}}


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



my @infiles = @ARGV;

for my $infas (@infiles) {
  my $seqio=Bio::SeqIO->new('-file' => $infas);

  while(my $seqobj=$seqio->next_seq()) {
    my $seqid = $seqobj->display_name();
    my $seqdesc = $seqobj->description();
    my $nt_seq;
    my $aa_seq;
    my $ppos=1;
    if($strand==1) {
      $nt_seq=$seqobj->seq();
      $aa_seq=$seqobj->translate(undef, undef, $frame, undef, undef, undef, undef)->seq();
    }
    elsif($strand==-1) {
      $aa_seq=$seqobj->revcom()->translate(undef, undef, $frame, undef, undef, undef, undef)->seq();
      $nt_seq=$seqobj->revcom()->seq();
    }
    &pretty($seqid, $seqdesc, $nt_seq, $aa_seq, 60, $frame);
  }
}

exit;

### subs begin ###

# {{{ sub pretty
sub pretty {
  my $seqid = shift(@_);
  my $seqdesc = shift(@_);
  my $ntseq = shift(@_);
  my $aaseq = shift(@_);
  my $width = shift(@_);
  my $frame = shift(@_);
  print(">$seqid $seqdesc\n");
  my $pwidth = $width/3;
  my @ntlines;
  my $pos = 0;
  while(my $chunk = substr($ntseq, $pos, $width)) {
    push(@ntlines, $chunk);
    $pos += $width;
  }
  my @trlines;
  $pos = 0;
  while(my $chunk = substr($aaseq, $pos, $pwidth)) {
    push(@trlines, $chunk);
    $pos += $pwidth;
  }
  my $four_lc = 0;
  foreach my $ntline (@ntlines) {
    print("$ntline\n");
    my $aastr = $trlines[$four_lc];
    print(" " x $frame);
    print(join("  ", split(//, $aastr)), "\n");
    print("\n");
    $four_lc += 1;
  }
}
# }}}
