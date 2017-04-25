#!/usr/bin/perl -w
# $ Id: $

# Definition file for runSelfTests.pl

# TEST CASES
# A list of all the test cases to be performed automatically
@main::TestCases = (
    "WaterBox",
    "GLY15_EV",
    "ProteinMoveset"
);

# REPLACEMENT ARRAYS
# key => value pairs of replacements
%main::Replacements = (
   "_ABSINTH_HOME_" => "/packages/absinth",
);


1;  # this is important
