#!/usr/bin/env perl

use strict;
use warnings;



my @dirs = qw (./txmon_new ./txmon_old);
my $prefix = "TXMon_txsec";
my $suffix = "root";


(my $prog = $0) =~ s|.+/||g;
my $help = "Usage: $prog [-flags]
Options:
-f filename       => use the Run Numbers in 'filename'
-h                => this Help screen
-i                => do Not convert eps to gifs
-w                => do Not generate html folders
-r                => do Not copy root files
";



my ($filename, $noconvert, $nofoldergenerator, $nocopy);
while (@ARGV && $ARGV[0] =~ /^-/)
  {
    my $arg = shift @ARGV;

    if ($arg =~ /^-{1,2}h/i)
    {
        die $help;
    }

    if ($arg =~ /^-f/i)
      {
        $filename = shift;
        next;
      }
    if ($arg =~ /^-i/i)
      {
        $noconvert = "true";
        next;
      }
    if ($arg =~ /^-w/i)
      {
        $nofoldergenerator = "true";
        next;
      }
    if ($arg =~ /^-r/i)
      {
        $nocopy = "true";
        next;
      }

    warn "I don't understand '$arg'.\n";
  }




my @runs;
foreach (@ARGV)
  {
    if (/^\d{6,}$/)
      {
        push @runs, $_;
      } # if digits
  } # foreach argv




# when converting this pl file into a py file
# please rewrite this section, entirely.
# this segment assumes that each run number is 6
# digits long. This will not be so for the first
# 99,999 runs. One should use copyRoot.py for
# obtaining the run numbers. A run file of 001
# has been added to the testdata file for trouble
# -shooting
# Happy Hunting, ~richard~
if ($filename)
  {
    open (SOURCE, $filename) or die;
    while (<SOURCE>)
      {
        chomp;
        s/\#.+//;
        s/\s+//g;
        if (/^\d{6,}$/)
	  {
            push @runs, $_;
	  }
      } # while source
  } # if filename

die "Usage: $0 123456 [123457] ...\n" unless (@runs);



@runs = reverse sort {$a <=> $b} @runs;
print "runs @runs\n";
my @rootfiles;
foreach my $run (@runs)
  {
    foreach my $dir (@dirs)
      {
        my @files = glob ("$dir/$prefix$run*$suffix");
        if (@files)
	  {
            push @rootfiles, @files;
            last;
	  } # if found
      } # foreach dir
  } # foreach run



print "files:\n",join("\n", @rootfiles),"\n";
if (@rootfiles)
  {
    my $convert;
    if (!$noconvert)
      {
        $convert = "; ./makeImgs.tcsh";
      }
    my $foldergenerator;
    if (!$nofoldergenerator)
      {
        $foldergenerator = "; ./htmldirectorygen.py ./L*/*.html";
      }
    my $rootcopier;
    if (!$nocopy)
      {
        $rootcopier = "; ./copyRoot.py $filename";
      }



    my $command = "(./AllTXMonFitter.exe @rootfiles $convert $foldergenerator $rootcopier) >& all.out &";
    print $command,"\n";
    system $command;

  }
