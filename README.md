![paPAML Logo](https://github.com/RetroWWU/paPAML/blob/main/paPAML.jpg)

# Description

paPAML simplifies, amplifies, and accelerates selection analyses via parallel processing, including detection of negatively selected sites. paPAML compiles results of site, branch, and branch-site models and detects site-specific negative selection with the output of a codon list labelling significance values. The tool simplifies selection analyses for casual and inexperienced users and accelerates computing speeds up to the number of allocated computer threads

# Installation

## FreeBSD (Version 13.*)

The operating system for which paPAML was first developed and tested for was FreeBSD. The installation on this OS is very simple - you don't need any conda environment (but if you need it: see later, there is a "quick" howto), all the needed packages are available in the system repositories. As root (or using sudo), install the following packages with the command:

    # pkg install wget paml hyphy p5-Proc-ProcessTable p5-Statistics-Distributions p5-File-Which

Now you can either (1) download paPAML.pl as a "normal" user directly from this github repository and place it wherever you like (best option might be your home directory $HOME/bin). If possible adjust the PATH environment to the directory where you put paPAML.pl into (though this is optional, it let’s you access paPAML easier). Or (2) use the following, potentially easier step:

    # wget https://raw.githubusercontent.com/RetroWWU/paPAML/main/paPAML.pl
    # chmod u+x paPAML.pl

and place the extracted directory in your home directory $HOME/bin or another target directory.

    # mkdir $HOME/bin
    # mv paPAML.pl $HOME/bin

That's it!

## Linux (Ubuntu 22.04)

In pricipal the installation is the same on Linux. First install the following packages as root (or using sudo):

    # apt update
    # apt install hyphy-common hyphy-pt paml libfile-which-perl libproc-processtable-perl libstatistics-distributions-perl

Now you can again either (1) download paPAML.pl as a "normal" user directly from this github repository and place it wherever you like (best option might be your home directory $HOME/bin). Additionally it would be good to adjust the PATH environment to the directory where you put paPAML.pl into. Or (2) use the following potentially easier step:

    # wget https://raw.githubusercontent.com/RetroWWU/paPAML/main/paPAML.pl
    # chmod u+x paPAML.pl

and place the extracted directory in your home directory $HOME/bin or another target directory. 

    # mkdir $HOME/bin
    # mv paPAML.pl $HOME/bin

Unfortunately binaries in Ubuntu (like /usr/bin/codeml and /usr/bin/hyphy) are only wrapper scripts with some standard definitions, which will cause problems. To circumvent this you have to adjust the PATH of the following directories before starting the program:

    # export PATH=/usr/lib/hyphy/bin:/usr/lib/paml/bin:$PATH

Copy-paste the export command into your shell resource file .bashrc to make it permanent.

That's it, too!

## macOS X

The installation under Conda in OSX seems not to be as easy as possible.  So here another way to get it working

* Install PAML

http://abacus.gene.ucl.ac.uk/software/paml.html

* Install HYPHY

https://www.hyphy.org/

* Install cpanm

    # curl -L https://cpanmin.us | perl - --sudo App::cpanminus
    # cpanm install Proc::ProcessTable
    # cpanm install File::Which
    # cpanm install Statistics::Distributions
    # cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)

* download paPAML.pl

    # wget https://raw.githubusercontent.com/RetroWWU/paPAML/main/paPAML.pl
    # chmod u+x paPAML.pl
    # mv paPAML.pl /usr/local/bin

* Add path in .bashrc

export PATH=$PATH:/usr/local/bin 

## Conda

The following coding is a quick solution to the conda environment, but it provides a usable conda container environment with name papaml and the main program paPAML.pl which is stored in $HOME/bin folder and is executable.

(1) Install the conda / anaconda environment and enter "yes" to write the changes to the shell

(2) Start or restart bash and execute the following commands (note: your prompt should like like: "(base)...")

    conda create -n papaml
    conda activate papaml
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install -y perl hyphy paml perl-file-which perl-statistics-distributions perl-app-cpanminus gcc make
    cpanm install Proc::ProcessTable

(3) Download the paPAML.pl program and put it as an executalbe in $HOME/bin

    wget https://raw.githubusercontent.com/RetroWWU/paPAML/main/paPAML.pl .
    mkdir $HOME/bin
    mv paPAML.pl $HOME/bin
    chmod u+x $HOME/bin/paPAML.pl

(4) Finally you can adjust the PATH in your $HOME/.bashrc that enables you to call paPAML.pl everywhere without a path:

    echo "export PATH=$HOME/bin:\$PATH" >> $HOME/.bashrc
    
(5) Potentially you will need to restart bash for the changes to take effect.

and it is done, too :-)

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
