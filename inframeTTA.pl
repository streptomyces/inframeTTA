#!/usr/bin/perl
use 5.14.0;
use utf8;
use Carp;
use lib qw(/home/sco /home/sco/mnt/smoke/perllib);
use File::Basename;
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
    tablistV tablistVE linelistV linelistVE tablistH linelistH
    tablistER tablistVER linelistER linelistVER tabhashER tabhashVER csvsplit);
use File::Spec;
use Bio::SeqIO;

# {{{ Getopt::Long
use Getopt::Long;
my $outdir;
my $indir;
my $fofn;
my $outex; # extension for the output filename when it is derived on infilename.
my $conffile = qq(local.conf);
my $errfile;
my $runfile;
my $outfile;
my $testCnt = 0;
our $verbose;
my $skip = 0;
my $help;
GetOptions (
"outfile:s" => \$outfile,
"dirout:s" => \$outdir,
"indir:s" => \$indir,
"fofn:s" => \$fofn,
"extension:s" => \$outex,
"conffile:s" => \$conffile,
"errfile:s" => \$errfile,
"runfile:s" => \$runfile,
"testcnt:i" => \$testCnt,
"skip:i" => \$skip,
"verbose" => \$verbose,
"help" => \$help
);
# }}}


if($help) {
exec("perldoc $0");
exit;
}

# {{{ open the errfile
if($errfile) {
open(ERRH, ">", $errfile);
print(ERRH "$0", "\n");
close(STDERR);
open(STDERR, ">&ERRH"); 
}
# }}}

# {{{ Populate %conf if a configuration file 
my %conf;
if(-s $conffile ) {
  open(my $cnfh, "<", $conffile);
  my $keyCnt = 0;
  while(my $line = readline($cnfh)) {
    chomp($line);
    if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
    my @ll=split(/\s+/, $line, 2);
    $conf{$ll[0]} = $ll[1];
    $keyCnt += 1;
  }
  close($cnfh);
  linelistE("$keyCnt keys placed in conf.");
}
elsif($conffile ne "local.conf") {
linelistE("Specified configuration file $conffile not found.");
}
# }}}

# {{{ outdir and outfile business.
my $ofh;
my $idofn = 0;    # Flag for input filename derived output filenames. 
if($outfile) {
  my $ofn;
  if($outdir) {
    unless(-d $outdir) {
      unless(mkdir($outdir)) {
        croak("Failed to make $outdir. Exiting.");
      }
    }
    $ofn = File::Spec->catfile($outdir, $outfile);
  }
  else {
    $ofn = $outfile;
  }
  open($ofh, ">", $ofn);
}
elsif($outdir) {
linelistE("Output filenames will be derived from input");
linelistE("filenames and placed in $outdir");
    unless(-d $outdir) {
      unless(mkdir($outdir)) {
        croak("Failed to make $outdir. Exiting.");
      }
    }
$idofn = 1;
}
else {
  open($ofh, ">&STDOUT");
}
select($ofh);
# }}}

# {{{ populate @infiles
my @infiles;
if(-e $fofn and -s $fofn) {
open(FH, "<", $fofn);
while(my $line = readline(FH)) {
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my $fn;
if($indir) {
$fn = File::Spec->catfile($indir, $line);
}
else {
$fn = $line;
}

push(@infiles, $fn);
}
close(FH);
}
else {
@infiles = @ARGV;
}

# }}}


for my $infile (@infiles) {

my $seqio=Bio::SeqIO->new(-file => $infile);
# my $seqout=Bio::SeqIO->new(-fh => $ofh, -format => 'fasta');


while(my $seqobj=$seqio->next_seq()) {
my $cdsCnt = 0;
  for my $feature ($seqobj->all_SeqFeatures()) {
    if($feature->primary_tag() eq 'CDS') {
      $cdsCnt += 1;
      my $codon_start = 1;
      if($feature->has_tag('codon_start')) {
        ($codon_start) = $feature->get_tag_values('codon_start');
      }

      my %anno;
      for my $tag (qw(gene note product locus_tag gene_synonym)) {
        if($feature->has_tag($tag)) {
          my @taglist = $feature->get_tag_values($tag);
          $anno{$tag} = \@taglist;
        }
      }
      my $subobj = $feature->spliced_seq(-nosort => 1);
      my $offset = 0;
      if($codon_start > 1) { $offset += $codon_start - 1; }
      my @ttapos = inframeTTA($subobj, $offset);
      if(@ttapos) {
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
        if(not exists($anno{product})) { push(@{$anno{product}}, "None"); }
        print(join("\t", $id, join(" ", @ttapos), @{$anno{product}}), "\n");
      }
    }
  }
}
}


exit;

# Multiple END blocks run in reverse order of definition.
END {
close($ofh);
close(STDERR);
close(ERRH);
}

sub inframeTTA {
  my $inseq=shift(@_);
  my $offset = shift(@_);
  my @retlist;
  my $seq;
  if(ref($inseq)) { $seq = lc($inseq->seq()); }
  else { $seq = lc($inseq); }
  my $pos = $offset;
  while(my $triplet = substr($seq, $pos, 3)) {
    if($triplet=~m/[tu][tu]a/i) {
      push(@retlist, $pos);
    }
    $pos+=3;
  }
  return(@retlist);
}



# # {{{ sub feat_translate {
# sub feat_translate {
#   my $feature=shift(@_);
#   my $codon_start=1;
#   if($feature->has_tag('codon_start')) {
#     ($codon_start) = $feature->get_tag_values('codon_start');
#   }
#   my $aaobj;
#   eval {
#     my $offset=1;
#     if($codon_start > 1) { $offset = $codon_start;}
#     my $featobj=$feature->spliced_seq(-nosort => 1);
#     $aaobj = $featobj->translate(-offset => $offset, -complete => 1);
#   };
#   return($aaobj);
# }
# # }}}
