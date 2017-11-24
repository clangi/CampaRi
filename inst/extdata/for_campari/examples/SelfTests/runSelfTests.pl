#!/usr/bin/perl -w
# $Id: runSelfTests.pl,v 1.1 2008/05/06 18:09:24 matt Exp $
# Utility to automatically run and analyze a series of consistency checks
# All the tests here are determinisic, meaning that we expect answers to
# be identical.

use strict;
use Getopt::Std;
use File::Copy;
$|++; # turn on autoflush

sub usage {
    my $reason = shift;
    my $data = <<DATA;

Usage: 
  runSelfTests.pl [-h][-v][-f def_file][-c][-d][-r][-e] [tests]

 Perform a series of consistency checks, and compare against "known
 good" results.

 -h Prints this help message.
 -v Verbose output
 -f Use given input file (default: TestCases.pl)
 -c Create and populate work directories
 -d Clears work directories
 -r Execute the tests
 -e Analyze tests results

 tests Names of the tests to perform.  If not specified, all tests
    listed in the input file are performed

DATA

    print $reason ."\n" if $reason;
    print $data;
    exit;
}


# recurses through all the test directories and sets up the work directory.
# To set up the work directory, creates (if it does not exist) it, and
# then copies with replace all the regular files in the test directory.
sub do_create_test {
   no warnings 'once';
   my @tests = @_ or die "Tests not passed";


   foreach my $test_d (@tests) { # loop through all test directories

      my @files;

      foreach my $f (glob("$test_d/*")) { # @files becomes a list of all regular files in $test_d
         push(@files, $f) if (-f $f);
      }

      my $work_d = $test_d.$main::work_dir;
      if (not -d $work_d) {
         printf "Making dir $work_d\n" if $main::verbose;
         mkdir $work_d or die "Cannot mkdir $work_d: $!";
      }

      # copy with replacement (using hash $main::Replacement) all @files
      foreach my $f (@files) {
         my $wf = $f;
         $wf =~ s/$test_d/$work_d/;  # new filename has dir replaced by work dir
         copy_replace($f, $wf, \%main::Replacements);
      }
   }
}


sub do_delete_test {
   my @tests = @_ or die "Tests not passed";

# safety warning, confirm
# from http://www.devdaily.com/perl/edu/articles/pl010005/pl010005.shtml
   printf "Delete all test results? (y/n) ";
   my $key;
   chomp ($key = <STDIN>);

   if ($key !~ /^y/) {
      printf "Quitting...\n";
      exit(1);
   }

   foreach my $test_d (@tests) { # loop through all test directories
      my $dir = $test_d.$main::work_dir;
      next if not -d $dir;
      printf "Deleting $dir/*, $dir\n" if $main::verbose;
      unlink <$dir/*>;
      rmdir $dir or die "Cannot rmdir $dir: $!\n";
   }
}

sub do_run_test {
   my @tests = @_ or die "Tests not passed";

   my $cmd = "sh run_test.sh";  # this will be executed in each directory
   foreach my $test_d (@tests) { # loop through all test directories
      print "Running $test_d test...";
      my $dir = $test_d.$main::work_dir;
      do_execute($dir, $cmd);
   }
}

      
# The files which we diff are defined in the file DiffFiles.pl, in each 
# test directory.  
sub do_evaluate_test {
   no warnings 'once';
   my @tests = @_ or die "Tests not passed";
   my $diff_def = "DiffFiles.pl"; # file with filenames to diff
   my $diff_fn = "diffs.txt";     # where output stored

   # this will be executed in each directory for each file to be compared against standard
   #my $cmd_prefix = "diff --side-by-side --suppress-common-lines ";  
   my $cmd_prefix = "diff -c --report-identical-files";  

   foreach my $test_d (@tests) { # loop through all test directories

      my $dfn = "$test_d/$diff_def";
      if (not -f $dfn) {
         die "Cannot find the input definition file, $dfn";
      }
      require $dfn;  # this is where @main::DiffFiles is defined

      my $work_d = $test_d.$main::work_dir;
   
      print "Evaluating test $test_d...\n";
      # delete the diff file if it exists
      $dfn = "$work_d/$diff_fn";
      unlink $dfn if (-e $dfn);

      my $diffs = 0;
      foreach my $f (@main::DiffFiles) {
         
         my $cmd = "$cmd_prefix $f ../$main::standard_dir/$f >> $diff_fn";
         my $status = do_execute($work_d, $cmd, 1);  # 0 on no diffs, 1 on diffs
         $diffs += $status;
         print "   Diffs found for $f in $test_d\n" if ($status == 1);
      }
      if ($diffs) {
         print "Please examine file $work_d/$diff_fn\n\n";
      } else {
         print "OK.\n\n";
      }
   }
}

# Executes given command in given directory, then comes back to current dir
# If integer error_flag is passed, then error is reported iff return value of
# command equals that value, otherwise anything except 0 is an error.
# returns the return value of the command
sub do_execute {
   my $dir = shift or die "dir not passed";
   my $cmd = shift or die "cmd not passed";
   my $err;
   $err = shift or $err = 0;

   printf "Executing $cmd in $dir\n" if $main::verbose;
   my $curdir = `pwd`;
   chomp $curdir;
   chdir $dir or die "Cannot chdir to $dir: $!\n";
   my $status = system $cmd;
   $status /= 256;
   die "$cmd exited funny: $?" unless ($status == 0 or $status == $err);
   chdir $curdir or die "Cannot chdir to $curdir: $!\n";
   return $status;
}

sub copy_replace {
   my $in_file = shift;
   my $out_file = shift;
   my $rr = shift; # this is a hash reference.
   my %reps = %$rr;

   printf "Copy with replace: $in_file -> $out_file\n" if $main::verbose;
   open(IN, "<$in_file") or die "Cannot open $in_file: $!\n";
   open(OUT, ">$out_file") or die "Cannot open $out_file: $!\n";

# copy the file line by line, attempting to do all replacements on all lines.
# if the lines are huge, this could be slow!
   while (my $line = <IN>) {
      foreach my $key (keys %reps) {
         $line =~ s/$key/$reps{$key}/g;
      }
      print OUT $line;
   }
   close IN;
   close OUT;
}
   

no warnings 'once';
$main::work_dir = "/work";
$main::standard_dir = "/standard";

# parameter flag stuff...
use vars qw($opt_h $opt_v $opt_f $opt_c $opt_d $opt_r $opt_e);
usage("Error processing parameter flags") if not
    getopts("hvf:cdre"); 
usage("Help requested (-h)") if $opt_h;

$main::verbose = $opt_v;
$main::def_file = $opt_f ? $opt_f : "TestCases.pl";

usage("Nothing to do!  Please specify one or more of the flags, -c -d -r -e") 
      unless ($opt_d or $opt_c or $opt_r or $opt_e);

# Read the Test Cases file.  see, Perl Cookbook p. 399
if (not -f $main::def_file) {
   die "Cannot find the input definition file, $main::def_file";
}
printf "Reading $main::def_file\n" if $main::verbose;
require $main::def_file;

# names of test directories can come from command line, or from the definition file
my @tests;
if ($ARGV[0]) { # tests passed on command line
   @tests = @ARGV;
} else {     # tests taken from def_file.  
   @tests = @main::TestCases;  
}

do_delete_test(@tests) if $opt_d;
do_create_test(@tests) if $opt_c;
do_run_test(@tests) if $opt_r;
do_evaluate_test(@tests) if $opt_e;
