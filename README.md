# paPAML

## Description

paPAML is a tool that combines parallel execution of codeml and hyphy for different trees

# Installation

## FreeBSD (version 13.*)

I know, not so many people use the brillant BSD operating system.  But it is definitly worth to take a look! The Instalaltion on FreeBSD is simple.  You don not need any conda environment, all the needed packages are available in the system repositories.  As root (or pre sudo) install the packages with

\# pkg add paml hyphy p5-Proc-ProcessTable p5-Statistics-Distributions

Then download paPAML.pl - and place it somewhere you like.  Additionally it would be best to adjust the PATH to the directory where you put paPAML.pl into.  That's it!  

Call paPAML.pl to get help

\# paPAML.pl

USAGE

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
                   
    significance - the maximum p value to print marked trees.  Used for

## Linux

