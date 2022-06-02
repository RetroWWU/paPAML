# Description

paPAML simplifies, amplifies, and accelerates selection analyses via parallel processing, including detection of negatively selected sites. paPAML compiles results of site, branch, and branch-site models and detects site-specific negative selection with the output of a codon list labelling significance values. The tool simplifies selection analyses for casual and inexperienced users and accelerates computing speeds up to the number of allocated computer threads

# Installation

## FreeBSD (Version 13.*)

The operating system, where paPAML was first developed and tested under FreeBSD, so there will be a first description how to install paPAML here.  The installation on FreeBSD is very simple - you don't need any conda environment, all the needed packages are available in the system repositories.  As root (or using sudo) install the packages with

    # pkg install wget paml hyphy p5-Proc-ProcessTable p5-Statistics-Distributions p5-File-Which

Now you can (1) download as a "normal" user paPAML.pl direct from the github repository and place it where you like, the best maybe is $HOME/bin.  Additionally it would be good to adjust the PATH environment to the directory where you put paPAML.pl into.  Or you (2) do following, what may be even more easy

    # wget https://raw.githubusercontent.com/RetroWWU/paPAML/main/paPAML.pl
    # chmod u+x paPAML.pl

and place it in $HOME/bin - if you want

    # mkdir $HOME/bin
    # mv paPAML.pl $HOME/bin

That's it!

## Linux (Ubuntu 22.04)

In pricipal the installation is the same on Linux. First install as root following packages

    # apt update
    # apt install hyphy-common hyphy-pt paml libfile-which-perl libproc-processtable-perl libstatistics-distributions-perl

Now you can (1) download as a "normal" user paPAML.pl direct from the github repository and place it where you like, the best maybe is $HOME/bin.  Additionally it would be good to adjust the PATH environment to the directory where you put paPAML.pl into.  Or you (2) do following, what may be even more easy

    # wget https://raw.githubusercontent.com/RetroWWU/paPAML/main/paPAML.pl
    # chmod u+x paPAML.pl

and place it in $HOME/bin - if you want

    # mkdir $HOME/bin
    # mv paPAML.pl $HOME/bin

But unfortunately are the binaries (/usr/bin/codeml and /usr/bin/hyphy) in Ubuntu only wrapper scripts with some standard definitions - so to make paPAML correctly you have to adjust the PATH

    # export PATH=/usr/lib/hyphy/bin:/usr/lib/paml/bin:$PATH

You may place it in your shell resource file .bashrc to make it work permanently

That's it, too!

# Usage

to get a small help about parameters and usage just type the command itself

    # paPAML.pl
    USAGE
        paPAML.pl -p runs [-f controlfiles] [-t tests] [-s significance] [-d] {codemlparams}
        paPAML.pl -i
        paPAML.pl -c

    VERSION 1.20

    WHERE
        runs         - the number of parallel runs
        ...

For more details please have a look into the paper :-)
