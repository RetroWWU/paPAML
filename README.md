# paPAML

## Description

paPAML simplifies, amplifies, and accelerates selection analyses via parallel processing, including detection of negatively selected sites. paPAML compiles results of site, branch, and branch-site models and detects site-specific negative selection with the output of a codon list labelling significance values. The tool simplifies selection analyses for casual and inexperienced users and accelerates computing speeds up to the number of allocated computer threads

# Installation

## FreeBSD (version 13.*)

The operating system, where paPAML was developed and tested first is/was FreeBSD, so there will be a first description how to install paPAML here.  The installation on FreeBSD is simple.  You don't need any conda environment, all the needed packages are available in the system repositories.  As root (or pre sudo) install the packages with

\# pkg install paml hyphy p5-Proc-ProcessTable p5-Statistics-Distributions

As a "normal" user download paPAML.pl and place it where you like, the best maybe in $HOME/bin.  Additionally it would be best to adjust the PATH to the directory where you put paPAML.pl into.  That's it!

You can also download paPAML.pl direct from the git repository.  First you need the git software

\# pkg install git

and then you just type

\# 

## Linux

# Usage

Call paPAML.pl to get help

`\# paPAML.pl
U`SAGE
    paPAML.pl -p runs [-f controlfiles] [-t tests] [-s significance] [-d] {codemlparams}
    paPAML.pl -i
    paPAML.pl -c

VERSION 1.20

WHERE
    runs         - the number of parallel runs
    controlfiles - a list of control files.  It is assumed they are named
                   with a suffix "".ctl"!  If not given all files with
                   that suffix are taken to be calculated
    tests        - the used tests (1, 2, 3 or h) to run the data against.
                   They can be written like "1" or "12" . The order does
                   not matter.
                   (default: 123h)
    significance - the maximum p value to print marked trees.  Used for`


