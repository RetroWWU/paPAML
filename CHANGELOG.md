**2026-05-15 - v2.16**

Remove bayes values in omega graphs - meaning to remove the misleading info and programming error "ARRAY(...)" from the files *.B und *.BS in *result.omega

**2026-04-21 - v2.15a**

The input data of sequences is checked against (1) If the seqeunces length is dividable by 3 and (2) not containing stop codons.  By default the stop codons are TAA, TAG and TGA.  This can be changed using "-seqtype 2" which sets them to AGA, AGG, TAA, TAG, TGA.

The logging is not written into a file paPAML-yyyy-mm-dd-hh-mm.log anymore, by default it is written to the console.  The parameter "-l" may be used to redirect to a logging file.

The deletion of symlinks to codeml and hyphy are (always) deleted on program end.

**2025-12-04 - v2.13**

The run writes a file called pid with the process id of paPAML.pl.  So it is easier to terminate a run - even if it would be running in the background like:

    # nohup paPAML.pl param1 param2 > run.log &

Terminate it (when entered the directory of choice) with

    # kill $(cat pid); rm pid

**2025-11-21 - v2.12**

enhance output by adding number of trees

    # +--------------------------------------------------------------------------------------+
    # |  Results Test 2 - branch-site specific                                               |
    # |  p-value_significance_limit: 0.050000 / corrected_for_multiple_testing: 0.001087     |
    # |  Tested number of trees: 46                                                          |
    # +--------------------------------------------------------------------------------------+

by making "Tested branches" for sigma and signam / #trees (here 0.05 and 0.001...)

    # ==> Tested branches with p < 0.050000
    # ...
    # ==> Tested branches with p < 0.001087
    # ...

Correct the codon counting and add the positive and negative counts

    # +--------------------------------------------------------------------------------------+
    # |  Results Test 4 - HyPhy FEL                                                          |
    # +--------------------------------------------------------------------------------------+
    # ...
    # 5       -       0.0000
    # 11      -       0.0001
    # 15      -       0.0401
    # ...
    # ==> Positive_selection: 4
    # ==> Negative_selection: 98

**2024-07-10 - v2.11**

add maxtime parameter -m

**2024-04-04 - v2.10**

add site values / [s] in branch-site specific svg graphs

**2024-03-07 - v2.9**

adjust path setting for seqfile and treefile

**2024-02-14 - v2.8**

Write logging into a file named paPAML-yyyy-mm-dd-hh-mm-ss.log

**2023-12-13**

(1) add tree to svg graphs and (2) make -icode work with different amino acids

**2023-01-20**

add *.result_aa.fa file for amino acids

**2023-01-16**

Remove a typo in code and skip invalid test 3 scenario

**2022-09-17**

Enhance output format and overall runtime message

**2022-06-07**

Correct tree with weights

**2022-06-01**

Add fasta sequences

**2022-05-30**

Enable termination of subprocesses on terminte

**2022-05-18**

Disable termination of subprograms by interrupt

**2022-05-09**

Remove double paramters in ctl file

**2022-05-05**

Change Model -> Test in result

**2022-04-11**

Redesign a bit and remove a small bug

**2022-04-01**

Remove a bug in getCodonInfo

**2022-03-23**

Continue revert sequence orientation and values

**2022-03-22**

Revert sequence orientation and values

**2022-03-17**

Add 3 lines in sequence output

**2022-02-13**

Add hyphy calculation

**2022-01-18**

Allow tree files with ambigous names

**2021-12-14**

handle codeml errors

**2021-12-13**

Tree name and tree in result

**2021-12-09**

Add additional parameters
