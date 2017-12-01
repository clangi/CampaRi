#!/usr/bin/perl -w
# $ Id: $

# Defines, in each test directory, the files which will be compared against the standard ones.

# COMPARE FILES
# Files which are diffed 
@main::DiffFiles = (
"ACCEPTANCE.dat",
"ENERGY.dat",
"TEST_END.int",
"TEST_END.pdb",
"TEST_START.int",
"TEST_START.pdb",
);

1;  # this is important
