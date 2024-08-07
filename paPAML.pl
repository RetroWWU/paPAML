#!/usr/bin/env perl

#
# ==============================================================================
# [2024-07-10] v2.11: add parameter -m (maxtime)
# [2024-04-04] v2.10: add [s] in svg graphs
# [2024-03-07] v2.9: adjust path setting for seqfile and treefile
# [2024-02-13] v2.8: add logging file
# [2023-01-11] v2.7: add -o all
# [2023-12-12] v2.6: Implement -icode dependend settings
# [2023-12-06] v2.5: SVG graphs with tree
# [2023-06-21] v2.4: Print corrected p-value in test header only
# [2023-03-28] v2.3: Add omega calculation
# [2023-01-20] v2.2: add *.result_aa.fa file for amino acids
# [2023-01-16] v2.1: correct a type and test 3 results
# [2022-09-15] v2.0: Extend runtime info
# [2022-06-07] v1.23: Correct tree with weights
# [2022-06-01] v1.22: Add fasta sequences
# [2022-05-30] v1.21: Enable termination of subprocesses on terminte
# [2022-05-18] v1.20: Disable termination of subprograms by interrupt
# [2022-05-09] v1.19: Remove double paramters in ctl file
# [2022-05-05] v1.18: Change Model -> Test in result
# [2022-04-11] v1.17: Redesign a bit and remove a small bug
# [2022-04-01] v1.15: Remove a bug in getCodonInfo
# [2022-03-23] v1.14: Continue revert sequence orientation and values
# [2022-03-22] v1.13: Revert sequence orientation and values
# [2022-03-17] v1.12: Add 3 lines in sequence output
# [2022-02-13] v1.11: Add hyphy calculation
# [2022-01-18] v1.10: Allow tree files with ambigous names
# [2021-12-14] v1.9: handle codeml errors
# [2021-12-13] v1.8: Tree name and tree in result
# [2021-12-09] Add additional parameters
# ==============================================================================
#

use strict;

use Cwd qw(cwd realpath);
use File::Basename;
use File::Path  qw(rmtree);
use File::Which qw(which);
use Proc::ProcessTable;
use Statistics::Distributions;

my $RUNTIMEFILE = "runtime";

# The start time - used in combination with maxtime
my $start = time;

# The logfile handle;
my $logfile;

# The tag/extension for the links of running codeml, hyphy by special name
my $tag       = time;
my $codemlpgm = "codeml-$tag";
my $hyphypgm  = "hyphy-$tag";

# Parameters
my $significance = 0.05;
my $para         = 0;
my $tests        = "123h";
my $debug        = 0;
my $omega;
my $maxtime;

my ($info, $clean);

my @ctlnames;

# The starttime of this run and the last time check and the interval.
# The lasttimecheck is put in "future" to have at least some time left
# for a prognose
my ($starttime, $lasttimecheck, $lasttimeinterval) = (time, time + 600, 600);

# The actual (total) number of runs
my $runindex = 0;

# The total runs for a tree / ctl file
my ($ctlcount, $ctlindex);

# Translate codons in amino acids
my %aacode = (
	AAA => "K",
	AAC => "N",
	AAG => "K",
	AAT => "N",
	ACA => "T",
	ACC => "T",
	ACG => "T",
	ACT => "T",
	AGA => "R",
	AGC => "S",
	AGG => "R",
	AGT => "S",
	ATA => "I",
	ATC => "I",
	ATG => "M",
	ATT => "I",
	CAA => "Q",
	CAC => "H",
	CAG => "Q",
	CAT => "H",
	CCA => "P",
	CCC => "P",
	CCG => "P",
	CCT => "P",
	CGA => "R",
	CGC => "R",
	CGG => "R",
	CGT => "R",
	CTA => "L",
	CTC => "L",
	CTG => "L",
	CTT => "L",
	GAA => "E",
	GAC => "D",
	GAG => "E",
	GAT => "D",
	GCA => "A",
	GCC => "A",
	GCG => "A",
	GCT => "A",
	GGA => "G",
	GGC => "G",
	GGG => "G",
	GGT => "G",
	GTA => "V",
	GTC => "V",
	GTG => "V",
	GTT => "V",
	TAA => "*",
	TAC => "Y",
	TAG => "*",
	TAT => "Y",
	TCA => "S",
	TCC => "S",
	TCG => "S",
	TCT => "S",
	TGA => "*",
	TGC => "C",
	TGG => "W",
	TGT => "C",
	TTA => "L",
	TTC => "F",
	TTG => "L",
	TTT => "F",
);

# Mapping icode => codon => amini acids
my %aacode2 = (
	1 => {
		ATA => "M",
		TGA => "W"
	},
	2 => {
		ATA => "M",
		CTC => "*",
		CTG => "*",
		CTT => "*",
		TGA => "W"
	},
	3 => {
		TGA => "W"
	},
	4 => {
		AGA => "S",
		AGG => "S",
		ATA => "M",
		TGA => "W"
	},
	5 => {
		TAA => "Q",
		TAG => "Q"
	},
	6 => {
		AAA => "N",
		AGA => "S",
		AGG => "S",
		TGA => "W"
	},
	7 => {
		TGA => "C"
	},
	8 => {
		CTG => "S"
	},
	9 => {
		AGA => "G",
		AGG => "G",
		ATA => "M",
		TGA => "W"
	},
	10 => {
		TAG => "Q"
	}
);

# Mapping from codeml -icode to hyphy --code
my %icodemap = (
	0  => "Universal",
	1  => "Vertebrate-mtDNA",
	2  => "Yeast-mtDNA",
	3  => "Mold-Protozoan-mtDNA",
	4  => "Invertebrate-mtDNA",
	5  => "Ciliate-Nuclear",
	6  => "Echinoderm-mtDNA",
	7  => "Euplotid-Nuclear",
	8  => "Alt-Yeast-Nuclear",
	9  => "Ascidian-mtDNA",
	10 => "Blepharisma-Nuclear",
);

# The default codeml parameters.  They can be changed over command line
my %params = (
	"CodonFreq"    => "2",
	"Malpha"       => "0",
	"Mgene"        => "0",
	"RateAncestor" => "0",
	"Small_Diff"   => "0.5e-6",
	"aaDist"       => "0",
	"alpha"        => "0",
	"cleandata"    => "1",
	"clock"        => "0",
	"estFreq"      => "0",
	"fix_alpha"    => "1",
	"fix_blength"  => "-1",
	"fix_rho"      => "1",
	"getSE"        => "0",
	"icode"        => "0",
	"method"       => "0",
	"ndata"        => "1",
	"omega"        => "1",
	"outfile"      => "mlc",
	"rho"          => "0",
	"runmode"      => "0",
	"seqtype"      => "1"
);

#
# ------------------------------------------------------------------------
# Display usage
# ------------------------------------------------------------------------
#
sub usage {
	my $p = join("\n", map {$_ = sprintf("    %-14s %s", "-$_", $params{$_})} (sort keys %params));

	print <<EOF;
USAGE
    paPAML.pl -p runs [-f controlfiles] [-t tests] [-s significance] [-o omega]
              [-d] [-m maxtime] {codemlparams}
    paPAML.pl -i [-f controlfiles]
    paPAML.pl -c

VERSION 2.11

WHERE
    runs         - the number of parallel runs
    controlfiles - a list of control files.  It is assumed they are named
                   with a suffix ".ctl"!  If not given all files with
                   that suffix are taken to be calculated
    tests        - the used tests (1, 2, 3 or h) to run the data against.
                   They can be written like "1" or "12". The order does
                   not matter.
                   (default: $tests)
    significance - the maximum p value to print marked trees.  Used for
                   printing bayes values and in hyphy call
                   (default: $significance)
    omega        - a single species or a comma separated species list like
                   "hom,sap,xyz" where the omega values are calculated or
                   "all" for all species
    maxtime      - the maximum runtime in minutes, hours or days.  If not
                   provided it is unlimited
                   (examples: "12m" or "34h" or "2d")
    -d           - the generated result directories are kept and not deleted
    -i           - info about your runs
    -c           - clean all temporary folders
    codemlparams - additional parameters for the ctl file, if not provided
                   the default parameters will be used.
                   (example: -Mgene 9 -rho 34)

DEFAULT CODEML-PARAMETERS
$p

DESCRIPTION
    This program takes the specified or all control files (*.ctl) of
    the actual folder and calls paml/codeml and hyphy (in paralell as
    background processes) for them.  There are several tests that are
    done:

    * test 1:
      Site specific (one run without marker)
      model = "0", nssites = "1 2 7 8", fixomega = "0", omega = "1"

    * test 2:
      Branch site model with selection (several runs with markers)
      model = "2", nssites = "2", fixomega = "0", omega = "1"
      Branch site model without selection (several runs with markers)
      model = "2", nssites = "2", fixomega = "1", omega = "1"

    * test 3:
      Branch model with selection (several runs with markers)
      model = "2", nssites = "0", fixomega = "0", omega = "1"
      Branch model without selection (several runs with markers)
      model = "2", nssites = "0", fixomega = "1", omega = "1"

    * test h:
      Run hyphy test

    For every test and all the marked trees there will be 1 run and 1
    (temporary) subdirectory with codeml or hyphy results created.
    When a run is finished the subfolder will contain a file called
    DONE or ERROR - marking finished jobs or those with an error.

    Assume your control file is called abc.ctl, these subfolders for
    all neccessary codeml runs are called like abc-00-00000,
    abc-01-00002 etc.  For hyphy the folder is abc-hyphy.  So the
    naming of the temporary codeml folders is in general

        ctlfilename "-" branchtype "-" treenumber

    and for hyphy

        ctlfilename "-" "hyphy"

    The branchtypes are "10", "20", "21", 30", "31" and "hyphy".  "10"
    is for the site specific branch. "20" and "21" for branch site
    model without or with fixomega.  "30" and "31" for branch model
    without or with fixomega.  The suffix "hyphy" is used for the
    hyphy run.  The treenumber is the "tree number".  Note: tree
    numbers start with zero. The folders abc-20-00012 and abc-30-00012
    correspond together, the first for branch site mode and the second
    for the branch model, both without fixomega.

    The -o (omega) specification has a special meaning.  It uses the
    model = "2", nssites = "0", fixomega = "0" setting for codeml and
    generates an additionl output file called *.result.omega with the
    dN/dN values of a specified species related branches.  The output is

      Tree dn/ds_background dn/ds_foreground
      ((1 #1,2),3) 0.40097 0.94020
      ((1,2 #2),3) 0.12 1.2
      ...

    The log output will have the form

      [T yyyy-mm-dd hh:mm:ss] Same message

    where T is an indicator.  "E" is Error, "I" is additional info and
    ">" marking a special event.  In the end you will get an info
    about the total runetime taken like

      [> 2022-09-17 11:12:11] The total runtime was 37.6 minutes

    To establish this the program calculates the runtime it took.  On
    an interruption the runtime is put in the file called "runtime".
    On continuation this file is added to the runtime the program
    took.  This will happen even several times.  If the run finnished
    successfully this file is removed.

    As already mentioned, a finished run in a subfolder will have a
    file called DONE - this is the indicator that the calculation for
    that specific run/tree is done.  On the other hand on an error
    there exists a file called ERROR.

    The program stays running, as long as your codeml or hyphy jobs
    are working.  If you press CTRL-C the program terminates and all
    running codeml and hyphy runs are canceled.

    When all runs finnished succesfully (or if in the meantime all
    needed jobs of a control file are done) there will be a *.result a
    *.result.fa and a *.result_aa.fa file for every *.ctl file found
    and the generated subfolders will be removed.  You can skip the
    deletion of the generated subfolders, when you use the -d (debug)
    parameter.

    Restarting paPAML.pl again after a termination or with error runs
    will first remove the subfolders without a DONE file and rerun
    those runs again.  Meaning: you can restart as many times as long
    as not all runs are done.

    If a *.result file exists for a *.ctl file, the run(s) will be
    skipped if it is started again.
    
    Additionally a logfile called paPAML-yyyy-mm-dd-hh-mm-ss.log is
    created where the logging will be written
EOF
	exit(0);
}

if (!@ARGV) {
	usage();
}

# ------------------------------------------------------------------------
# Open and close the logfile
# ------------------------------------------------------------------------

{
	my ($sec, $min, $hour, $mday, $mon, $year) = localtime(time);
	my $filename = sprintf("paPAML-%4d-%02d-%02d-%02d-%02d-%02d.log", $year + 1900, $mon + 1, $mday, $hour, $min, $sec);
	open($logfile, ">", $filename);
	print "Run is started!  See $filename for logging information...\n";
	message("I", "The command is: " . basename($0) . " " . join(" ", @ARGV));
}

END {
	close($logfile);
}

#
# ------------------------------------------------------------------------
# Returns the date as a string
# ------------------------------------------------------------------------
#
sub getDate {
	my ($sec, $min, $hour, $mday, $mon, $year) = localtime(time);
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec);
}

#
# ------------------------------------------------------------------------
# Prints a message
# ------------------------------------------------------------------------
#
sub message {
	my ($type, $message) = @_;
	printf $logfile ("[%1s %s] %s\n", $type ? $type : "-", getDate(), $message);
}

#
# ------------------------------------------------------------------------
# Clean all generated files for control files not needed any more
# ------------------------------------------------------------------------
#
sub cleanRuns {
	for my $ctlname (@_) {
		message(">", "Clean $ctlname...");
		my @dirs = grep {-d $_} <$ctlname-*>;
		for my $dir (@dirs) {
			message("", "Clean run $dir...");
			rmtree($dir);
		}
	}
}

#
# ------------------------------------------------------------------------
# Clean undone runs
# ------------------------------------------------------------------------
#
sub cleanUndones {
	for my $ctlname (@_) {
		my @dirs = grep {-d $_} <$ctlname-*>;
		for my $dir (@dirs) {
			if (!-f "$dir/DONE") {
				message(">", "Clean error run $dir...");
				rmtree($dir);
			}
		}
	}
}

#
# ------------------------------------------------------------------------
# Clean program symlinks
# ------------------------------------------------------------------------
#
sub cleanSymlinks {
	for my $pgm ($codemlpgm, $hyphypgm) {
		if (-l $pgm) {
			message(">", "Clean symlink $pgm");
			unlink($pgm);
		}
	}
}

#
# ------------------------------------------------------------------------
# Get parameters from command line
# ------------------------------------------------------------------------
#
sub getParams {
	for (my $i = 0 ; $i < @ARGV ; $i++) {
		my $p = $ARGV[$i];
		if ($p eq "-i") {
			$info = 1;
		}
		elsif ($p eq "-c") {
			$clean = 1;
		}
		elsif ($p eq "-p") {
			$para = $ARGV[++$i];
			if (!($para =~ m/^\d+$/)) {
				message("E", "Parameter -p needs an integer value!");
				exit(1);
			}
		}
		elsif ($p eq "-f") {
			my $s      = $ARGV[++$i];
			my @a      = split(/,/, $s);
			my @errors = ();
			for my $x (@a) {
				my $file = "$x.ctl";
				push(@errors, $file) if (!-f $file);
			}
			if (@errors) {
				message("E", "Control file(s) " . join(",", @errors) . " do not exist!");
				exit(1);
			}
			@ctlnames = @a;
		}
		elsif ($p eq "-t") {
			my $s = $ARGV[++$i];
			$tests = "";
			$tests .= "1" if ($s =~ m/1/);
			$tests .= "2" if ($s =~ m/2/);
			$tests .= "3" if ($s =~ m/3/);
			$tests .= "h" if ($s =~ m/h/);
			if (!($tests =~ m/^[123h]+$/)) {
				message("E", "Parameter -t needs value(s) of 1, 2, 3 and/or h!");
				exit(1);
			}
			if ($tests ne $s) {
				message("I", "Parameter tests (-t) is changed to $tests!");
			}
		}
		elsif ($p eq "-s") {
			$significance = $ARGV[++$i];
			if (!($significance =~ m/^\d*\.?\d+$/)) {
				message("E", "Parameter -s needs a float value!");
				exit(1);
			}
		}
		elsif ($p eq "-o") {
			$omega = $ARGV[++$i];
			$tests .= "3" if (!($tests =~ m/3/));
		}
		elsif ($p eq "-d") {
			$debug = 1;
		}
		elsif ($p eq "-m") {
			$maxtime = $ARGV[++$i];
			if ($maxtime =~ m/^(\d+)([mhd])$/) {
				my ($value, $dim) = ($1, $2);
				$maxtime = $value * 60;
				$maxtime *= 60 if ($dim ne "m");
				$maxtime *= 24 if ($dim eq "d");
			}
			else {
				message("E", "Parameter -m is incorrect!");
				exit(1);
			}
		}
		elsif ($p eq "-runmode") {
			my $runmode = $ARGV[++$i];
			if ($runmode != 0) {
				message("W", "Runmode <> 0 may produce no or invalid results!");
			}
		}
		else {
			if ($p =~ m/^-/) {
				my $value = $ARGV[++$i];
				$params{substr($p, 1)} = $value;
			}
		}
	}

	my $count = 0;
	$count++ if ($clean);
	$count++ if ($info);
	$count++ if ($para);
	if ($count != 1) {
		message("E", "The parameters are incorrect!");
		exit(1);
	}
}

#
# ------------------------------------------------------------------------
# Write into a file
# ------------------------------------------------------------------------
#
sub writeFile {
	my ($filename, $content) = @_;

	open(F, ">", $filename);
	print F $content;
	close(F);
}

#
# ------------------------------------------------------------------------
# Read a file and return the lines
# ------------------------------------------------------------------------
#
sub readFile {
	my ($filename) = @_;

	my @lines;
	open(F, "<", $filename);
	while (my $line = <F>) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		push(@lines, $line);
	}
	close(F);

	return @lines;
}

#
# ------------------------------------------------------------------------
# Get the sub process ids
# ------------------------------------------------------------------------
#
sub getSubpids {
	my @pids;

	my $proc = Proc::ProcessTable->new();
	foreach my $p (@{$proc->table}) {
		if ($p->uid == $<) {
			if ($p->cmndline =~ m/^..\/($codemlpgm|$hyphypgm)/) {
				push(@pids, $p->pid);
			}
		}
	}

	return @pids;
}

#
# ------------------------------------------------------------------------
# Kills all subprocesses
# ------------------------------------------------------------------------
#
sub terminate {
	message(">", "Terminate...");

	my $runtime = -f $RUNTIMEFILE ? join("", readFile($RUNTIMEFILE)) : 0;
	writeFile($RUNTIMEFILE, $runtime + (time - $starttime));

	my @pids = getSubpids();
	for my $pid (@pids) {
		message("", "process $pid...");
		kill(15, $pid) if (kill(0, $pid));
	}

	sleep(5);

	@pids = getSubpids();
	for my $pid (@pids) {
		message("", "process $pid...");
		kill(9, $pid) if (kill(0, $pid));
	}

	cleanSymlinks();

	exit(1);
}

#
# ------------------------------------------------------------------------------
# Checks maximum runtime
# ------------------------------------------------------------------------------
#
sub checkMaxtime {
	if ($maxtime && $start + $maxtime < time) {
		message("I", "Maximum runtime is over...");
		terminate();
	}
}

#
# ------------------------------------------------------------------------------
# Returns the success state
# ------------------------------------------------------------------------------
#
sub getSuccess {
	return sprintf("%.1f%%", $ctlindex / $ctlcount * 100.0);
}

#
# ------------------------------------------------------------------------------
# Returns the estimated runtime
# ------------------------------------------------------------------------------
#
sub getFinnish {
	my $t = time - $starttime;
	return sprintf("%d minutes", (($t * $ctlcount / $ctlindex) - $t) / 60 + 1);
}

#
# ------------------------------------------------------------------------------
# Wait for the next free slot
# ------------------------------------------------------------------------------
#
sub dowait {
	while (1) {
		checkMaxtime();
		if ($lasttimecheck + $lasttimeinterval < time) {
			message("I", sprintf("Estimated time to finnish all runs: %s", getFinnish()));
			$lasttimecheck = time;
		}
		my @pids = getSubpids();
		return if (@pids < $para);
		sleep(5);
	}
}

#
# ------------------------------------------------------------------------------
# Print the tree file of a directory
# ------------------------------------------------------------------------------
#
sub printTree {
	my ($dir, $treeno) = @_;

	my @fs = <$dir/*.ctl>;
	if (!@fs) {
		print RESULT "[E] No control file in $dir!\n";
		return;
	}

	my $filename;
	for my $line (readFile($fs[0])) {
		if ($line =~ m/treefile\s*=\s*([^\s]+)/) {
			$filename = $1;
		}
	}

	printf RESULT ("%s\t", "Tree_" . ($treeno + 1));
	my @a = readFile("$dir/$filename");
	for (my $i = 1 ; $i < @a ; $i++) {
		print RESULT $a[$i], "\n";
	}
}

#
# ------------------------------------------------------------------------------
# Look for errors in the end
# ------------------------------------------------------------------------------
#
sub printErrors {
	for my $ctlname (@ctlnames) {
		my @dirs = grep {-d $_} <$ctlname-*>;

		my $error = 0;
		for my $dir (@dirs) {
			if (-f "$dir/ERROR") {
				message("E", "There are errors in run $dir!");
				$error = 1;
			}
			my @a = qx(grep "is missing in the tree" $dir/codeml.log) if (-f "$dir/codeml.log");
			for my $s (@a) {
				chomp($s);
				message("E", "Error from codeml run $dir: $s");
				$error = 1;
			}
		}
	}
}

#
# ------------------------------------------------------------------------------
# Returns the text/comment inside braces
# ------------------------------------------------------------------------------
#
sub enbrace {
	return join("\n", map {$_ = sprintf("|  %-82s  |", $_); $_} @_);
}

#
# ------------------------------------------------------------------------------
# Returns omega values dn, ds ad the index of the first matching line
# ------------------------------------------------------------------------------
#
sub getOmega {
	my ($test, $lines) = @_;

	my ($dn, $ds, $index) = (0, 0, -1);

	if ($test eq "2") {
		my (@p, @b, @f);
		my $i = 0;
		while ($i < @$lines) {
			if ($lines->[$i++] =~ m/^site\s+class\s+/) {
				$index = $i - 1;
				last;
			}
		}
		while ($i < @$lines) {
			my $line = $lines->[$i++];
			if ($line =~ m/^proportion\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)/) {
				@p = ($1, $2, $3, $4);
			}
			elsif ($line =~ m/^background\s+w\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)/) {
				@b = ($1, $2, $3, $4);
			}
			elsif ($line =~ m/^foreground\s+w\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)/) {
				@f = ($1, $2, $3, $4);
			}
			else {
				last;
			}
		}
		for (my $i = 0 ; $i < 4 ; $i++) {
			$dn += $p[$i] * $b[$i];
			$ds += $p[$i] * $f[$i];
		}
	}
	elsif ($test eq "3") {
		for (my $i = 0 ; $i < @$lines ; $i++) {
			if ($lines->[$i] =~ m/^w \(dN\/dS\)[^:]+:\s+([\d\.]+)\s+([\d\.]+)/) {
				($dn, $ds) = ($1, $2);
				$index = $i;
				last;
			}
		}
	}

	return ($dn, $ds, $index);
}

#
# ------------------------------------------------------------------------------
# Prints codeml data
# ------------------------------------------------------------------------------
#
sub generateCodeml {
	my ($ctlname) = @_;

	my ($bayes12, $bayes78, $codondata) = ({}, {}, {});

	print RESULT <<EOS;
+--------------------------------------------------------------------------------------+
|  Codeml test results                                                                 |
+--------------------------------------------------------------------------------------+
EOS

	# Calculate model 1
	if ($tests =~ m/1/) {
		print RESULT <<EOS;

+--------------------------------------------------------------------------------------+
|  Results Test 1 - site models                                                        |
+--------------------------------------------------------------------------------------+
|  EXAMPLE                                                                             |
|                                                                                      |
|  p-value_significance_limit: 0.05                                                    |
|                                                                                      |
|  Model_1_M0              23                      -1549.837332     0.6728             |
|  Model_2_M               25                      -1545.728933     1.0040             |
|  Model_7_M0              23                      -1549.873442     0.7000             |
|  Model_8_M               25                      -1545.726929     1.0033             |
|                                                                                      |
|  Site models (1,2,7,8)   np (degrees of freedom)  lnl-values      w (omega)          |
|                                                                   -values (dn/ds)    |
|                                                                                      |
|  Four site models, model 1 and 7 allow for neutral and negative selection            |
|  (the null models), while model 2 and 8 also allow for positive selection            |
|  and are compared with their counterparts (2 with 1 and 8 with 7).                   |
|                                                                                      |
|  Test_1_2                8.216798                2                0.016434           |
|  Test_7_8                8.293026                2                0.015819           |
|  CHI-square              difference lnl x2       Difference np    P-value            |
|                                                                                      |
|  The Model 1-2 comparison is usually more stringent than the model 7-8 comparison    |
|  and therefore not as sensitive.                                                     |
|                                                                                      |
|  Bayes_1_2                                BEB                                        |
|  76                       P               0.959*           4.837      +/- 2.203      |
|  Bayes_7_8                                                                           |
|  13                       H               0.972*           4.025      +/- 1.882      |
|  76                       P               0.978*           4.053      +/- 1.874      |
|  Codon number             AA (reference)  P-values sign.   w-values   SD of w-values |
|                                            *P > 95%                                  |
|                                           **P > 99%                                  |
|  Sites under significant positive selection found with Bayes Empirical Bayes (BEB)   |
|  analysis, both for model 1-2 and model 7-8 comparison.                              |
|                                                                                      |
|  Conclusions: In this example, the site model comparisons show a significant         |
|  difference from the null models (P-values =  0.016434 and 0.015819, respectively)   |
|  and found two sites under significant positive selection across the phylogeny       |
|  (13 H and 76 P). 13 H was only discovered by the less stringent and more sensitive  |
|  model 7-8 comparison.                                                               |
+--------------------------------------------------------------------------------------+

EOS

		my @dirs  = <$ctlname-10-*>;
		my @lines = readFile("$dirs[0]/mlc");

		my @a = ();
		for my $line (@lines) {
			if ($line =~ m/lnL\(.*?np:\s*(\d+)\):\s*([0-9\.\-]+)/) {
				push(@a, [$1, $2]);
			}
		}

		my $index = 0;
		for (my $i = 0 ; $i < @lines ; $i++) {
			if ($lines[$i] =~ m/branch\s+t\s+N\s+S/) {
				my $line = $lines[$i + 2];
				push(@{$a[$index++]}, (split(/\s+/, $line))[4]) if ($index < @a);
			}
		}

		my @ms = ("1_M0", "2_M", "7_M0", "8_M");
		$index = 0;
		for my $x (@a) {
			printf RESULT ("%s\t%s\n", "Model_" . $ms[$index++], join("\t", @$x));
		}

		my $dltr = abs(2 * ($a[0]->[1] - $a[1]->[1]));
		my $dn   = abs($a[1]->[0] - $a[0]->[0]);
		printf RESULT ("%s\t%f\t%d\t%f\n", "Test_1_2", $dltr, $dn, Statistics::Distributions::chisqrprob($dn, $dltr));

		$dltr = abs(2 * ($a[2]->[1] - $a[3]->[1]));
		$dn   = abs($a[3]->[0] - $a[2]->[0]);
		printf RESULT ("%s\t%f\t%d\t%f\n", "Test_7_8", $dltr, $dn, Statistics::Distributions::chisqrprob($dn, $dltr));

		$index = 0;
		for (my $i = 0 ; $i < @lines ; $i++) {
			if ($lines[$i] =~ m/Bayes Empirical Bayes/) {
				$i += 6;
				print RESULT ($index == 0 ? "Bayes_2\n" : "Bayes_7\n");
				while ($lines[$i] =~ m/\d+\s+[A-Z\*]\s+/) {
					my @b = split(/\s+/, $lines[$i++]);
					if ($b[2] =~ m/\*/) {
						printf RESULT ("%d\t%s\t%s\t%.3f\t%s\t%.3f\n", $b[0], $b[1], $b[2], $b[3], $b[4], $b[5]);
						$b[2] =~ s/\*+$//g;
						($index == 0 ? $bayes12->{$b[0] - 1} : $bayes78->{$b[0] - 1}) = 1 - $b[2];
					}
				}
				$index++;
			}
		}
	}

	# Calculate model 2 and 3
	for my $test ("2", "3") {
		next if (!($tests =~ m/$test/));

		# Get folders for "without or with" fixomega
		my @dirs0 = <$ctlname-${test}0-*>;
		my @dirs1 = <$ctlname-${test}1-*>;

		# Write header if no data
		if (!@dirs0) {
			print RESULT ($test == 2 ? "\n# Test 2 - branch-site specific" : "\n# Test 3 - branch specific"), "\n\n";
		}

		# If there is a result and header is missing
		my $headerprinted = 0;

		# Loop over all directories (trees)
		for (my $treeno = 0 ; $treeno < @dirs0 ; $treeno++) {
			my @lines0 = readFile("$dirs0[$treeno]/mlc");
			my @lines1 = readFile("$dirs1[$treeno]/mlc");

			my (@b0,  @b1);
			my ($np0, $lnl0);
			my ($np1, $lnl1);

			# Get Bayes values
			for (my $k = 0 ; $k < @lines0 ; $k++) {
				if ($lines0[$k] =~ m/Bayes Empirical Bayes/) {
					$k += 2;
					while ($lines0[$k] =~ m/\d+\s+[A-Z\*]\s+/) {
						my @a = split(/\s+/, $lines0[$k++]);
						push(@b0, [$a[0], $a[1], $a[2]]) if ($a[2] =~ m/\*/);
					}
					last;
				}
			}
			for (my $k = 0 ; $k < @lines1 ; $k++) {
				if ($lines1[$k] =~ m/Bayes Empirical Bayes/) {
					$k += 2;
					while ($lines1[$k] =~ m/\d+\s+[A-Z\*]\s+/) {
						my @a = split(/\s+/, $lines1[$k++]);
						push(@b1, [$a[0], $a[1], $a[2]]) if ($a[2] =~ m/\*/);
					}
					last;
				}
			}

			# Get np and lnl values
			for my $line (@lines0) {
				if ($line =~ m/lnL\(.*?np:\s*(\d+)\):\s*([0-9\.\-]+)/) {
					($np0, $lnl0) = ($1, $2);
					last;
				}
			}
			for my $line (@lines1) {
				if ($line =~ m/lnL\(.*?np:\s*(\d+)\):\s*([0-9\.\-]+)/) {
					($np1, $lnl1) = ($1, $2);
					last;
				}
			}

			next if (($np0 == $np1) || ($lnl0 == $lnl1));

			my $dltr = abs(2 * ($lnl1 - $lnl0));
			my $dn   = abs($np1 - $np0);
			my $p    = Statistics::Distributions::chisqrprob($dn, $dltr);

			# Write header if first tree
			if (!$headerprinted) {
				my $s = enbrace(
					sprintf(
						"p-value_significance_limit: %f / corrected_for_multiple_testing: %f",
						$significance, $significance / @dirs0
					)
				);
				if ($test == 2) {
					print RESULT <<EOS;

+--------------------------------------------------------------------------------------+
|  Results Test 2 - branch-site specific                                               |
$s
+--------------------------------------------------------------------------------------+
|  EXAMPLE                                                                             |
|                                                                                      |
|  p-value_significance_limit: 0.05                                                    |
|  p-value_significance_limit (corrected_for_multiple_testing): 0.00455                |
|  correction: 0.05/11 (11 = number of foreground branches/tested trees                |
|                                                                                      |
|  Tree_3 ((Hsa_Human,Hla_Gibbon) #1,((Cgu/Can_colobus, Pne_langur),Mmu_rhesus),       |
|  (Ssc_squirrelM,Cja_marmorset));                                                     |
|                                                                                      |
|  analyzed tree with foreground branch labeled #1                                     |
|                                                                                      |
|  Tree_3_M                16                      -898.514392                         |
|  omega free (0 < w)      np (degrees of freedom) lnl value                           |
|                                                                                      |
|  Bayes_3_M                                                                           |
|   79                     L                        0.963*                             |
|  122                     R                        0.956*                             |
|                                                                                      |
|  codon number            amino acid (reference)   P-value and significance           |
|                                                   of selection on site               |
|                                                    *P > 95%                          |
|                                                   **P > 99%                          |
|                                                                                      |
|  Sites under significant positive sel. found with Bayes Empirical Bayes analysis     |
|                                                                                      |
|  Tree_3_M0               15                      -901.979225                         |
|  Omega fixed (w=1)       np (degrees of freedom) lnl value                           |
|                                                                                      |
|  Test_3                  6.929666                 1                 0.000706         |
|  CHI-square              difference lnl x2        difference np     P-value          |
|                                                                                      |
|  Tree_3_site_classes                                                                 |
|  Site class              0          1          2a          2b                        |
|  Proportion              0.39899    0.47265    0.05875     0.06960                   |
|  Background_w            0.00000    1.00000    0.00000     1.00000                   |
|  Foreground_w            0.00000    1.00000  999.00000   999.00000                   |
|                                                                                      |
|  Different site classes for background (B) and foreground (F) branch(es) with        |
|  corresponding w-values                                                              |
|                                                                                      |
|  Site class nomenclature     0: 0 < w < 1     for B and F                            |
|                              1: w = 1         for B and F                            |
|                             2a: 0 < w < 1     for B and w ≧ 1 for F                  |
|                             2b: w = 1         for B and w ≧ 1 for F                  |
|                                                                                      |
|  Proportion of overall sites of specific site classes (combined equal to 1)          |
|                                                                                      |
|  Tree_3_branch_omega     0.542250                 128.694300                         |
|  w-values of different   overall w-value of B     Overall w-value of F               |
|  branch types overall                                                                |
|                                                                                      |
|  Conclusions: Significant positive selection (P = 0.043569) and specific sites       |
|  under positive selection (79 L and 122 R) found in the foreground branch (exact     |
|  overall w-value difficult to indicate due to the upper limit of 999 being reached   |
|  in the site classes).                                                               |
+--------------------------------------------------------------------------------------+

EOS
				}
				else {
					print RESULT <<EOS;

+--------------------------------------------------------------------------------------+
|  Results Test 3 - branch model - branch specific                                     |
$s
+--------------------------------------------------------------------------------------+
|  EXAMPLE                                                                             |
|                                                                                      |
|  p-value_significance_limit: 0.05                                                    |
|  p-value_significance_limit (corrected_for_multiple_testing): 0.0125                 |
|  correction: 0.05/4 (4 = number of foreground branches/tested trees)                 |
|                                                                                      |
|  Tree_3 (pon_abe,(pan_tro,hom_sap #1));                                              |
|                                                                                      |
|  analyzed tree with foreground branch labeled #1                                     |
|                                                                                      |
|  Tree_3_M                7                      -1882.816024                         |
|  Tree_3_M0               6                      -1884.852625                         |
|  omega free (0 < w)      np (degrees of freedom) lnl value                           |
|  omega fixed to (w = 1)  np (degrees of freedom) lnl value                           |
|                                                                                      |
|  Test_3          4.073202            1              0.043569                         |
|  CHI-square      difference lnl x2   difference np  P-value                          |
|                                                                                      |
|  Tree_3_branch_omega     0.008030               1.015790                             |
|  w-values of different   w-value                w-value                              |
|  branch types            background             foreground branch                    |
|                                                                                      |
|  Conclusions: In this example, the w-value of #1 is significantly different          |
|  (P = 0.043569) than the background w. Nevertheless, after applying the correction   |
|  for multiple testing (new P-boundary = 0.0125), the w-value  (P = 0.043569) is no   |
|  longer significant.                                                                 |
+--------------------------------------------------------------------------------------+

EOS
				}
				$headerprinted = 1;
			}

			if ($p < $significance) {
				printTree($dirs0[$treeno], $treeno);

				printf RESULT ("%s\t%d\t%s\n", "Tree_" . ($treeno + 1) . "_M", $np0, $lnl0);
				if (@b0) {
					printf RESULT ("Bayes_%d_O\n", $treeno + 1);
					print RESULT join("\n", map {sprintf(("%d\t%s\t%s", $_->[0], $_->[1], $_->[2]))} @b0), "\n";
					map {$codondata->{$_->[0] - 1}->{"2"} .= sprintf(",Tree_%d:%0.3f", ($treeno + 1), 1 - $_->[2])} @b0;
				}

				printf RESULT ("%s\t%d\t%s\n", "Tree_" . ($treeno + 1) . "_M0", $np1, $lnl1);
				if (@b1) {
					printf RESULT ("Bayes_%d_M\n", $treeno + 1);
					print RESULT join("\n", map {sprintf(("%d\t%s\t%s", $_->[0], $_->[1], $_->[2]))} @b1), "\n";
					map {$codondata->{$_->[0] - 1}->{"2"} .= sprintf(",Tree_%d:%0.3f", ($treeno + 1), 1 - $_->[2])} @b1;
				}

				printf RESULT ("%s\t%f\t%d\t%f\n", "Test_" . ($treeno + 1), $dltr, $dn, $p);

				# Calculate background/foreground values

				my ($dn, $ds, $index) = getOmega($test, \@lines0);
				if ($test eq "2") {
					printf RESULT ("Tree_%d_site_classes\n", $treeno + 1);
					for (my $i = 0 ; $i < 4 ; $i++) {
						my $line = $lines0[$index + $i];
						$line =~ s/site class/site_class/;
						$line =~ s/d w/d_w/;
						$line =~ s/[ ]+/\t/g;
						print RESULT $line, "\n";
					}
				}
				printf RESULT ("Tree_%d_branch_omega\t%f\t%f\n", $treeno + 1, $dn, $ds);
			}
		}
	}

	return ($bayes12, $bayes78, $codondata);
}

#
# ------------------------------------------------------------------------------
# Prints hyphy data
# ------------------------------------------------------------------------------
#
sub generateHyphy {
	my ($ctlname, $codondata) = @_;

	print RESULT <<EOS;

+--------------------------------------------------------------------------------------+
|  Results Test 4 - HyPhy FEL                                                          |
+--------------------------------------------------------------------------------------+
|  EXAMPLE                                                                             |
|                                                                                      |
|  P-value_significance_limit: 0.05                                                    |
|                                                                                      |
|  codon number (reference)   "+"/"-"                      P-value                     |
|  66                            -                         0.0068                      |
|  76                            +                         0.028                       |
|                                positive selection(+)                                 |
|                                negative selection(-)                                 |
|                                                                                      |
|  Conclusions: In this example, the HyPhy FEL algorithm detected one site with        |
|  significant negative selection (P = 0.0068) and one site with significant           |
|  positive selection (P = 0.0288) across the phylogeny.                               |
+--------------------------------------------------------------------------------------+

EOS

	my $lineno = 0;
	my @lines  = readFile("$ctlname-hyphy/out");
	while ($lineno < @lines && !($lines[$lineno] =~ m/For partition 1 these sites/))           {$lineno++}
	while ($lineno < @lines && !($lines[$lineno] =~ m/(\e\[[0-9;]*m(?:\e\[K)?|[\x80-\xff]+)/)) {$lineno++}
	while ($lineno < @lines && $lines[$lineno] =~ m/(\e\[[0-9;]*m(?:\e\[K)?|[\x80-\xff]+)/) {
		my $line = $lines[$lineno++];
		$line =~ s/\s+//g;
		my @a = split(/\|/, $line);
		last if (@a < 5);
		my $pos = $a[6] =~ m/Pos/;
		$a[6] =~ m/p=(\d\.\d+)/;
		my $value = $1;
		$codondata->{int($a[1]) - 1}->{$pos ? "h+" : "h-"} = $value;
		printf RESULT ("%d\t%s\t%0.4f\n", $a[1], $pos ? "+" : "-", $value);
	}
}

#
# ------------------------------------------------------------------------------
# Prints the final sequence table and the fasta file
# ------------------------------------------------------------------------------
#
sub generateSequence {
	my ($ctlname, $bayes12, $bayes78, $codondata) = @_;

	# Get sequenz filename name from control file
	my @lines = grep {$_ =~ /seqfile/} readFile("$ctlname.ctl");
	my ($tmp, $seqfile) = split(/=/, $lines[0]);
	$seqfile =~ s/\s+//g;

	# Extract sequence
	my $content = join("\n", readFile($seqfile));
	my $seq;
	if ($content =~ m/>.*?\n(([ACGT\-]+\n)+)/is) {
		$seq = $1;
	}
	$seq =~ s/[\-\s]//g;
	$seq = uc($seq);

	print RESULT <<EOS;

+--------------------------------------------------------------------------------------+
|  Sequence overview of site specific results                                          |
+--------------------------------------------------------------------------------------+
|  Codon   Codon  Amino  T1_Bayes_1_2  T1_Bayes_7_8  T2_Bayes  Hyphy     Hyphy         |
|  number         acid                                         negative  positive      |
+--------------------------------------------------------------------------------------+

EOS

	my ($b12,  $b78,  $b,  $hn,  $hp);
	my ($b12_, $b78_, $b_, $hn_, $hp_);

	my $count = length($seq) / 3;
	my $icode = $params{icode};

	# Loop over codons
	for (my $i = 0 ; $i < $count ; $i++) {
		my $codon = substr($seq, $i * 3, 3);
		my ($cu, $cl) = (uc($codon), lc($codon));
		my $u  = uc($codon);
		my $aa = $aacode2{$icode}{$u} ? $aacode2{$icode}{$u} : $aacode{$u};
		my ($aau, $aal) = ($aa, lc($aa));
		my $cd    = $codondata->{$i};
		my $trees = substr($cd->{"2"}, 1);

		my $x_b12 = exists $bayes12->{$i};
		my $x_b78 = exists $bayes78->{$i};
		my $x_hm  = exists $cd->{"h-"};
		my $x_hp  = exists $cd->{"h+"};

		printf RESULT (
			"%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			$i + 1,
			$codon,
			$aacode2{$icode}{$codon} ? $aacode2{$icode}{$codon}         : $aacode{$codon},
			$x_b12                   ? sprintf("%0.4f", $bayes12->{$i}) : ".",
			$x_b78                   ? sprintf("%0.4f", $bayes78->{$i}) : ".",
			$trees                   ? $trees                           : ".",
			$x_hm                    ? sprintf("%0.4f", $cd->{"h-"})    : ".",
			$x_hp                    ? sprintf("%0.4f", $cd->{"h+"})    : "."
		);

		$b12 .= $x_b12 ? $cu : $cl;
		$b78 .= $x_b78 ? $cu : $cl;
		$b   .= $trees ? $cu : $cl;
		$hn  .= $x_hm  ? $cu : $cl;
		$hp  .= $x_hp  ? $cu : $cl;

		$b12_ .= $x_b12 ? $aau : $aal;
		$b78_ .= $x_b78 ? $aau : $aal;
		$b_   .= $trees ? $aau : $aal;
		$hn_  .= $x_hm  ? $aau : $aal;
		$hp_  .= $x_hp  ? $aau : $aal;
	}

	open(FASTA, ">", "$ctlname.result.fa");
	print FASTA qq(>T1_Bayes_1_2\n$b12\n);
	print FASTA qq(>T1_Bayes_7_8\n$b78\n);
	print FASTA qq(>T2_Bayes\n$b\n);
	print FASTA qq(>Hyphy_negative\n$hn\n);
	print FASTA qq(>Hyphy_positive\n$hp\n);
	close(FASTA);

	open(FASTA, ">", "$ctlname.result_aa.fa");
	print FASTA qq(>T1_Bayes_1_2\n$b12_\n);
	print FASTA qq(>T1_Bayes_7_8\n$b78_\n);
	print FASTA qq(>T2_Bayes\n$b_\n);
	print FASTA qq(>Hyphy_negative\n$hn_\n);
	print FASTA qq(>Hyphy_positive\n$hp_\n);
	close(FASTA);
}

#
# ------------------------------------------------------------------------------
# Return the tree if species matches
# ------------------------------------------------------------------------------
#
sub isOmega {
	my ($tree, $species) = @_;

	# Species and #1 direct together
	if ($tree =~ m/[^A-Za-z0-9\_\-]$species #1/) {
		return 1;
	}

	# Species inside balanced braces and a following #1
	elsif ($tree =~ m/^(.*?\)) #1/) {
		my $s = $1;
		my $n = 0;
		for (my $i = length($s) - 1 ; $i >= 0 ; $i--) {
			my $c = substr($s, $i, 1);
			if ($c eq ")") {
				$n++;
			}
			elsif ($c eq "(") {
				$n--;
			}
			if ($n == 0) {
				my $t = substr($s, $i);
				return $t =~ m/[^A-Za-z0-9\_\-]$species[^A-Za-z0-9\_\-]/ ? 1 : 0;
			}
		}
	}

	return 0;
}

# ------------------------------------------------------------------------------
# Constants for svg
# ------------------------------------------------------------------------------

my $FONT     = qq(font-family="monospace" font-size="12" font-weight="normal");
my $LINETYPE = qq(fill-opacity:1.0; stroke-linecap:square; stroke-opacity:1.0; stroke-width:1);
my $RECTTYPE = qq(fill-opacity:1.0; stroke-opacity:1.0; stroke-width:1);

# The x offset of "layers" in the tree graph
my $XDELTA = 40;

# Graph settings
my $LEFTOFFSET    = 30;
my $TOPPADDING    = 100;
my $LEFTPADDING   = 10;
my $RIGHTPADDING  = 40;
my $BOTTOMPADDING = 10;
my $CELLWIDTH     = 52;
my $HEIGHT        = 200;

my %COLORS = (
	white   => "rgb(255,255,255)",
	black   => "rgb(32,32,32)",
	gray    => "rgb(128,128,128)",
	icolor  => "rgb(64,64,64)",
	scolor  => "rgb(32,32,192)",
	hscolor => "rgb(192,32,32)"
);

#
# ------------------------------------------------------------------------------
# Return a float as readable string
# ------------------------------------------------------------------------------
#
sub getFloat {
	my $s;
	if ($_[0] < 0.0001) {
		$s = sprintf("%0.10e", $_[0]);
		$s =~ s/0+e/e/g;
	}
	else {
		$s = sprintf("%0.10f", $_[0]);
		$s =~ s/0+$//g;
	}
	return $s;
}

#
# ------------------------------------------------------------------------------
# Return the whole svg data
# ------------------------------------------------------------------------------
#
sub svgCreate {
	my ($width, $height, $elems) = @_;
	my $xmlns =
		qq(xmlns="http://www.w3.org/2000/svg")
	  . qq( xmlns:svg="http://www.w3.org/2000/svg")
	  . qq( xmlns:xlink="http://www.w3.org/1999/xlink");
	return qq(<svg height="$height" width="$width" $xmlns>\n) . join("\n", @$elems) . qq(\n</svg>);
}

#
# ------------------------------------------------------------------------------
# Return an svg text element
# ------------------------------------------------------------------------------
#
sub svgText {
	my ($x, $y, $text, $color) = @_;
	return qq(<text fill="$color" $FONT x="$x" y="$y">$text</text>);
}

#
# ------------------------------------------------------------------------------
# Return a vertical svg text element
# ------------------------------------------------------------------------------
#
sub svgTextUp {
	my ($x, $y, $text, $color) = @_;
	$x += 11;
	return qq(<text fill="$color" $FONT transform="translate($x,$y) rotate(-90)">$text</text>);
}

#
# ------------------------------------------------------------------------------
# Return an svg line eleemtn
# ------------------------------------------------------------------------------
#
sub svgLine {
	my ($x1, $y1, $x2, $y2, $color) = @_;
	return qq(<line style="fill:$color; stroke:$color; $LINETYPE" x1="$x1" y1="$y1" x2="$x2" y2="$y2"/>);
}

#
# ------------------------------------------------------------------------------
# Return an svg rectangle element
# ------------------------------------------------------------------------------
#
sub svgRectangle {
	my ($x, $y, $width, $height, $color) = @_;
	return qq(<rect style="fill:$color; stroke:$color; $RECTTYPE" width="$width" height="$height" x="$x" y="$y"/>);
}

#
# ------------------------------------------------------------------------------
# Return an svg circle element
# ------------------------------------------------------------------------------
#
sub svgCircle {
	my ($x, $y, $radius, $color) = @_;
	return qq(<circle style="fill:$color; stroke:$color; $RECTTYPE" cx="$x" cy="$y" r="$radius"/>);
}

#
# ------------------------------------------------------------------------------
# Return an y coordiante
# ------------------------------------------------------------------------------
#
sub svgY {
	my ($value, $height, $unit, $toppaddng) = @_;
	return int($height - $unit * $value + $toppaddng);
}

#
# ------------------------------------------------------------------------------
# Fill the y coordinate info into nodes, who do not have (that are no leaves)
# ------------------------------------------------------------------------------
#
sub svgFillY {
	my ($node) = @_;

	for my $n (@{$node->{subnodes}}) {
		svgFillY($n);
	}

	if (!$node->{y}) {
		my $first = $node->{subnodes}->[0];
		my $last  = $node->{subnodes}->[-1];
		if ($first->{y} && $last->{y}) {
			$node->{y1} = $first->{y};
			$node->{y2} = $last->{y};
			$node->{y}  = int(($first->{y} + $last->{y}) / 2);
		}
	}
}

#
# ------------------------------------------------------------------------------
# Fill the name/text elements in elems
# ------------------------------------------------------------------------------
#
sub svgFillNames {
	my ($node, $elems, $maxdepth) = @_;

	my $x = $LEFTPADDING + ($maxdepth + 1) * $XDELTA + 4;

	if ($node->{name}) {
		my $s = $node->{name} . ($node->{marked} ? " " . $node->{marked} : "");
		push(@$elems, svgText($x, $node->{y} + 4, $s, $COLORS{black}));
	}
	else {
		map {svgFillNames($_, $elems, $maxdepth)} @{$node->{subnodes}};
	}
}

#
# ------------------------------------------------------------------------------
# Fill the line elements in elems
# ------------------------------------------------------------------------------
#
sub svgFillLines {
	my ($node, $elems, $maxdepth) = @_;

	return if ($node->{name});

	push(@$elems, svgCircle($node->{x}, $node->{y}, 2, $COLORS{black}));

	if ($node->{marked}) {
		push(@$elems, svgText($node->{x} + 4, $node->{y} + 4, $node->{marked}, $COLORS{black}));
	}
	for my $n (@{$node->{subnodes}}) {
		svgFillLines($n, $elems, $maxdepth);
		push(@$elems, svgLine($n->{x}, $n->{y1}, $n->{x}, $n->{y2}, $COLORS{black}));
		if ($n->{name}) {
			my $x2 = $node->{x} + ($maxdepth - $n->{depth} + 1) * $XDELTA;
			push(@$elems, svgLine($node->{x}, $n->{y}, $x2, $n->{y}, $COLORS{black}));
		}
		else {
			my $x2 = $node->{x} + $XDELTA;
			push(@$elems, svgLine($node->{x}, $n->{y}, $x2, $n->{y}, $COLORS{black}));
		}
	}
}

#
# ------------------------------------------------------------------------------
# Fill the tree into elems
# ------------------------------------------------------------------------------
#
sub svgFillTree {
	my ($otree, $elems, $maxdepth, $y) = @_;

	my $depth   = 0;
	my $yoffset = 16;
	my $node    = {x => $LEFTPADDING + 2};

	while ($otree =~ m/(\(|,|\)(\s+#\d+)?|[a-z][^,:\(\) ]*(\s+#\d+)?|:[^\)\(, ]+)/gi) {
		my $x = $1;
		if ($x eq "(") {
			$depth++;
			my $subnode = {parent => $node, x => $LEFTPADDING + $depth * $XDELTA};
			push(@{$node->{subnodes}}, $subnode);
			$node     = $subnode;
			$maxdepth = $depth if ($depth > $maxdepth);
		}
		elsif ($x =~ m/^\)\s*(#\d+)?/) {
			my $marked = $1;
			$node->{marked} = $marked         if ($marked);
			$node           = $node->{parent} if ($node->{parent});
			$depth--;
		}
		elsif ($x =~ m/([a-z][^,:\(\) ]*)(\s+#\d+)?/i) {
			my ($name, $marked) = ($1, $2);
			my $subnode = {name => $name, x => $LEFTPADDING + $depth * $XDELTA, y => $y, depth => $depth};
			$subnode->{marked} = $marked if ($marked);
			push(@{$node->{subnodes}}, $subnode);
			$y += $yoffset;
		}
	}

	svgFillY($node);
	svgFillNames($node, $elems, $maxdepth);
	svgFillLines($node, $elems, $maxdepth);

	return ($y, $maxdepth * $XDELTA + 200);
}

#
# ------------------------------------------------------------------------------
# Create the omega / svg file
# ------------------------------------------------------------------------------
#
sub generateOmegaGraph {
	my ($file, $data, $title, $s, $hs, $otree, $test) = @_;

	my $y;
	my @elems;
	my @bayes;

	# The maximum omega value
	my $max = 2.0;

	for my $d (@$data) {
		$TOPPADDING = 200 if ($d->[3] >= $max);
	}

	my $unit = $HEIGHT / $max;

	# Width of the graph - not of the whole image
	my $width = $LEFTOFFSET + (@$data * $CELLWIDTH);
	$width = 650 if ($width < 650);

	push(@elems, svgText($LEFTPADDING + $LEFTOFFSET, 15, $title, $COLORS{black}));
	push(@elems, svgLine($LEFTPADDING, $TOPPADDING,           $LEFTPADDING, $TOPPADDING + $HEIGHT, $COLORS{black}));
	push(@elems, svgLine($LEFTPADDING, $TOPPADDING + $HEIGHT, $width,       $TOPPADDING + $HEIGHT, $COLORS{black}));
	push(@elems, svgTextUp($LEFTPADDING - 7, $TOPPADDING - 5, "Omega", $COLORS{black}));
	push(@elems, svgText($LEFTPADDING + $width - 4, $TOPPADDING + $HEIGHT + 4, "Time", $COLORS{black}));

	# Draw omega coordinates
	my $n = "1" . "0" x (length(sprintf("%d", $max)) - 1);
	for (my $o = 0 ; $o < $max ; $o += $n) {
		my $y = svgY($o, $HEIGHT, $unit, $TOPPADDING);
		push(@elems, svgText($LEFTPADDING + 5, $y - 5, sprintf("%d", $o), $COLORS{black}));
		push(@elems, svgLine($LEFTPADDING, $y, $width, $y, $COLORS{gray}));
	}

	# Draw values
	my $x = $LEFTOFFSET;
	for my $d (reverse @$data) {
		my $color;
		if ($d->[4] == 1 || $d->[3] == 999) {
			$color = $COLORS{gray};
		}
		elsif ($d->[4] >= $s) {
			$color = $COLORS{icolor};
		}
		else {
			$color = ($d->[4] <= $hs ? $COLORS{hscolor} : $COLORS{scolor});
		}

		my $y = svgY($d->[3] < $max ? $d->[3] : $max, $HEIGHT, $unit, $TOPPADDING);
		push(@elems, svgTextUp($x - 2,      $y - 5,  "[#]",   $color));
		push(@elems, svgTextUp($x - 2,      $y - 30, $d->[0], $color));
		push(@elems, svgTextUp($x - 2 + 12, $y - 5,  "[w]",   $color));
		push(@elems, svgTextUp($x - 2 + 12, $y - 30, $d->[3], $color));
		push(@elems, svgTextUp($x - 2 + 24, $y - 5,  "[p]",   $color));
		push(@elems, svgTextUp($x - 2 + 24, $y - 30, $d->[4], $color));
		if ($test == 2) {
			push(@elems, svgTextUp($x - 2 + 36, $y - 5, "[s]", $color));
			my @bs  = ();
			my @bs2 = ();
			if (@{$d->[5]} <= 3) {
				map {$_->[2] =~ s/[^\*]*//g; push(@bs, join("", @$_))} @{$d->[5]};
			}
			else {
				map {$_->[2] =~ s/[^\*]*//g; push(@bs, join("", @$_))} @{$d->[5]}[0 .. 2];
				push(@bs, "...");
				map {$_->[2] =~ s/[^\*]*//g; push(@bs2, join("", @$_))} @{$d->[5]};
				push(@bayes, sprintf("[#] %d [s] %s", $d->[0], join("/", @bs2)));
			}
			push(@elems, svgTextUp($x - 2 + 36, $y - 30, @bs ? join("/", @bs) : "-", $color));
		}
		elsif ($test == 3) {
			my $fb;
			$fb = "F&gt;B" if ($d->[3] > $d->[2]);
			$fb = "F=B"    if ($d->[3] == $d->[2]);
			$fb = "F&lt;B" if ($d->[3] < $d->[2]);
			push(@elems, svgTextUp($x - 2 + 36, $y - 5,  "[d]", $color));
			push(@elems, svgTextUp($x - 2 + 36, $y - 30, $fb,   $color));
		}
		push(@elems, svgRectangle($x, $y - 2, $CELLWIDTH - 4, 4, $color));
		$x += $CELLWIDTH;
	}

	$y = $TOPPADDING + $HEIGHT + 20;
	$s = sprintf("Significance: %f\n", $significance);
	$s .= sprintf(" / High Significance (p-value limit / number of tests): %f\n\n", $hs);
	push(@elems, svgText($LEFTPADDING, $y, $s, $COLORS{black}));

	$y += 4;
	for my $b (@bayes) {
		$y += 16;
		push(@elems, svgText($LEFTPADDING, $y, $b, $COLORS{black}));
	}

	$y += 4;
	my @a = (
		"Red: Significant with multiple testing correction",
		"Blue: Significant with original p-value limit (0.05)",
		"Black: Not significant",
		"Grey: Not significant and nonsensical data points (w=999 or p=1)",
		"#: Number of foreground branch starting with target species",
		"w: Omega value of foreground branch",
		"p: p-value of branch model comparison (foreground branch (F) against background branches (B))"
	);
	push(@a, "d: Direction of selection change (F&lt;B: towards neg. selection, F&gt;B: towards pos. selection)")
	  if ($test == 3);
	push(@a, "s: Relevant site values, separated by /")
	  if ($test == 2);

	for my $s (@a) {
		$y += 16;
		push(@elems, svgText($LEFTPADDING, $y, $s, $COLORS{gray}));
	}

	my ($h, $w) = svgFillTree($otree, \@elems, 0, $y + 20);
	$width = $w if ($w > $width);

	open(F, ">", $file) || die;
	binmode F;
	print F svgCreate($width + $LEFTPADDING + $RIGHTPADDING + 60, $h + $BOTTOMPADDING, \@elems);
	close F;
}

#
# ------------------------------------------------------------------------------
# Prints the species / omega based data and create a graph
# ------------------------------------------------------------------------------
#
sub generateOmega {
	my ($ctlname, $test) = @_;

	my $dir = "$ctlname.result.omega";
	mkdir($dir) if (!-d $dir);

	# Get folders without fixomega
	my @dirs0 = <$ctlname-${test}0-*>;
	my @dirs1 = <$ctlname-${test}1-*>;

	# Read omega tree
	my $otree;
	my @a = readFile("$ctlname.ctl");
	for (my $i = 0 ; !$otree && $i < @a ; $i++) {
		if ($a[$i] =~ m/treefile\s*=\s*([^\s]+)/) {
			$otree = join("", readFile($1));
		}
	}

	# Get all species from the omega parameter
	my @specs;
	if ($omega eq "all") {
		while ($otree =~ m/([a-z].*?)[^a-z0-9\_\-]/gi) {
			push(@specs, $1);
		}
	}
	else {
		@specs = split(",", $omega);
	}

	# Extract marked tree and values.  Do it before the loop over species,
	# to avoid uneccessary runtime

	my @mtrees;
	my @values = ();

	for (my $treeno = 0 ; $treeno < @dirs0 ; $treeno++) {
		my $mtree;
		my @a = readFile("$dirs0[$treeno]/codeml.ctl");
		for (my $i = 0 ; $i < @a ; $i++) {
			if ($a[$i] =~ m/treefile\s*=\s*([^\s]+)/) {
				$mtree = (readFile("$dirs0[$treeno]/$1"))[1];
				last;
			}
		}

		my @lines0 = readFile("$dirs0[$treeno]/mlc");
		my @lines1 = readFile("$dirs1[$treeno]/mlc");

		my ($dn, $ds) = getOmega($test, \@lines0);

		my @b0;

		# Get Bayes values
		for (my $k = 0 ; $k < @lines0 ; $k++) {
			if ($lines0[$k] =~ m/Bayes Empirical Bayes/) {
				$k += 2;
				while ($lines0[$k] =~ m/\d+\s+[A-Z\*]\s+/) {
					my @a = split(/\s+/, $lines0[$k++]);
					push(@b0, [$a[0], $a[1], $a[2]]) if ($a[2] =~ m/\*/);
				}
				last;
			}
		}

		my ($np0, $lnl0);
		my ($np1, $lnl1);

		# Get np and lnl values
		for my $line (@lines0) {
			if ($line =~ m/lnL\(.*?np:\s*(\d+)\):\s*([\d\.\-]+)/) {
				($np0, $lnl0) = ($1, $2);
				last;
			}
		}
		for my $line (@lines1) {
			if ($line =~ m/lnL\(.*?np:\s*(\d+)\):\s*([\d\.\-]+)/) {
				($np1, $lnl1) = ($1, $2);
				last;
			}
		}

		my $p;
		if (($np0 != $np1) && ($lnl0 != $lnl1)) {
			my $dltr = abs(2 * ($lnl1 - $lnl0));
			$p = Statistics::Distributions::chisqrprob(abs($np1 - $np0), $dltr);
		}
		else {
			$p = 1;
		}

		push(@values, [getFloat($dn), getFloat($ds), getFloat($p), \@b0]);
		push(@mtrees, $mtree);
	}

	# High significance
	my $hs = $significance / @dirs0;

	# Error species
	my %errspecs;

	for my $species (@specs) {
		my $diff = 0;
		my $ot   = $otree;
		my @data = ();

		for (my $i = 0 ; $i < @mtrees ; $i++) {
			my $mtree = $mtrees[$i];

			# Check and get omega values
			if (isOmega($mtree, $species)) {
				my $v = $values[$i];
				push(@data, [@data + 1, $mtree, $v->[0], $v->[1], $v->[2], $v->[3]]);

				# Change next " #1" in otree by values
				my $s = " #" . @data;
				substr($ot, index($mtree, " #1") + $diff, 0) = $s;
				$diff += length($s);
			}
		}

		if (!@data) {
			$errspecs{$species} = 1;
			next;
		}

		my $omegafile = "$dir/${species}." . ($test == 2 ? "BS" : "B");
		open(F, ">", $omegafile);
		printf F ("# Significance: %f\n", $significance);
		printf F (
			"# High Significance (corrected for multiple testing via p-value limit divided by number of tests): %f\n\n",
			$hs
		);
		print F "# Number\tTree\tdN/dS_background\tdN/dS_foreground\tP-Value\n\n";

		for my $d (@data) {
			print F join("\t", @$d), "\n";
		}
		close(F);

		writeFile("$omegafile.tree", $ot);

		my $title = "Forground Omega Graph - " . ($test == 2 ? "branch-site specific" : "branch specific");
		generateOmegaGraph("$omegafile.svg", \@data, $title, $significance, $hs, $ot, $test);
	}

	if (%errspecs) {
		message(
			"I",
			sprintf(
				"Following omega %s species have no data or are invalid: %s",
				($test == 2 ? "branch-site specific" : "branch specific"),
				join(",", sort keys %errspecs)
			)
		);
	}
}

#
# ------------------------------------------------------------------------------
# Look for finished control files and generate results
# ------------------------------------------------------------------------------
#
sub generate {
  LOOP:
	for my $ctlname (@ctlnames) {

		# Check for results file - don't regenerate
		my $resultfile = "$ctlname.result";
		next if (-f $resultfile);

		# Nothing done so far
		my @dirs = grep {-d $_} <$ctlname-*>;
		next if (!@dirs);

		# Check for errors or error run
		for my $dir (@dirs) {
			my @a = qx(grep "is missing in the tree" $dir/codeml.log) if (-f "$dir/codeml.log");
			next LOOP                                                 if (-f "$dir/ERROR" || @a);
		}

		# Check if all runs are ready
		my $count = 0;
		map {$count++ if (-f "$_/DONE")} @dirs;
		next if ($count < @dirs);

		message(">", "Generate result for $ctlname...");

		open(RESULT, ">", $resultfile);

		my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime();
		my $s = enbrace(
			sprintf("%04d-%02d-%02d %02d:%02d", $year + 1900, $mon + 1, $mday, $hour, $min),
			sprintf("Results for %s.ctl with tests %s and significance %f", $ctlname, $tests, $significance)
		);
		print RESULT <<EOS;
+--------------------------------------------------------------------------------------+
$s
+--------------------------------------------------------------------------------------+

EOS

		my ($bayes12, $bayes78, $codondata) = generateCodeml($ctlname);
		generateHyphy($ctlname, $codondata);
		generateSequence($ctlname, $bayes12, $bayes78, $codondata);

		close(RESULT);

		if ($omega) {
			generateOmega($ctlname, "2");
			generateOmega($ctlname, "3");
		}

		if (!$debug) {
			cleanRuns($ctlname);
		}
	}
}

#
# ------------------------------------------------------------------------------
# Mark the first/next element in the tree with "#1" and return the new
# tree.  Return empty, if the tree is finished, meaning the #1 is at
# the last node/element.
# ------------------------------------------------------------------------------
#
sub mark {
	my ($tree) = @_;

	my $index = index($tree, " #1");
	if ($index < 0) {
		$tree =~ s/([a-z][a-z0-9\_\-]+|\))([\d\.\:\[\]]+|)?/$1$2 #1/i;
	}
	else {
		$tree =~ m/^(.*) #1(.*)$/;
		my ($head, $tail) = ($1, $2);
		$tail =~ s/([a-z][a-z0-9\_\-]+|\))([\d\.\:\[\]]+|)?/$1$2 #1/i;
		$tree = "$head$tail";
		$tree = "" if ($tree =~ m/#1;$/ || !($tree =~ m/#1/));
	}

	return $tree;
}

#
# ------------------------------------------------------------------------------
# Extend the control file by parameters
# ------------------------------------------------------------------------------
#
sub extendCtl {
	my ($ctl, $model, $nssites, $fixomega) = @_;
	my @a = grep {!($_ =~ m/(\s+|)(model|NSsites|fix_omega)\s*=/)} split(/\n/, $ctl);
	return "model = $model\nNSsites = $nssites\nfix_omega = $fixomega\n" . join("\n", @a);
}

#
# ------------------------------------------------------------------------------
# Prepare and start a single run
# ------------------------------------------------------------------------------
#
sub runCodeml {
	my ($ctlname, $seqfile, $treefile, $ctl, $tree, $treetype, $treeno) = @_;

	my $subdir = sprintf("%s-%02d-%05d", $ctlname, $treetype, $treeno);
	return if (-f "$subdir/DONE");

	generate() if ((++$runindex % 10) == 0);

	dowait();

	$ctlindex++;
	message(">", sprintf("Start %s-%02d-%05d [started or finished %s]", $ctlname, $treetype, $treeno, getSuccess()));

	mkdir($subdir);
	chdir($subdir);

	my $sf = basename($seqfile);
	my $tf = basename($treefile);

	symlink($seqfile, $sf);
	writeFile("codeml.ctl", $ctl);
	writeFile($tf,          " 1\n$tree\n");

	my $command = "../$codemlpgm 2>&1 >codeml.log";
	writeFile("RUN", "($command; if [ \$? -eq 0 ]; then touch DONE; else touch ERROR; fi)&");
	system("sh RUN paPAML $subdir");

	chdir("..");
}

#
# ------------------------------------------------------------------------------
# Prepare and start a hyphy run
# ------------------------------------------------------------------------------
#
sub runHyphy {
	my ($ctlname, $alignfile, $treefile) = @_;

	my $subdir = "$ctlname-hyphy";
	return if (-f "$subdir/DONE");

	generate() if ((++$runindex % 10) == 0);

	dowait();

	$ctlindex++;
	message(">", sprintf("Start %s [started or finished %s]", $subdir, getSuccess()));

	mkdir($subdir);
	chdir($subdir);

	my $af = basename($alignfile);
	my $tf = basename($treefile);

	symlink($alignfile, $af);
	symlink($treefile,  $tf);

	my $code = $icodemap{$params{icode}};
	if (!$code) {
		message("E", "Codon usage table is not implemented in codeml and/or HyPhy FEL");
		return;
	}
	$code = "--code $code" if ($code);

	my $command = "../$hyphypgm fel $code --pvalue $significance --alignment $af --tree $tf >out 2>/dev/null";
	writeFile("RUN", "($command; if [ \$? -eq 0 ]; then touch DONE; else touch ERROR; fi)&");
	system("sh RUN paPAML $subdir");

	chdir("..");
}

#
# ------------------------------------------------------------------------------
# Returns the total number of runs depending on the tree
# ------------------------------------------------------------------------------
#
sub getRuns {
	my ($tree) = @_;

	my $count = 0;

	my ($n, $t) = (0, $tree);
	for (my $i = 0 ; $t = mark($t) ; $i++) {
		$n++;
	}
	$count++ if ($tests =~ m/1/);
	$count += 2 * $n if ($tests =~ m/2/);
	$count += 2 * $n if ($tests =~ m/3/);
	$count++ if ($tests =~ m/h/);

	return $count;
}

#
# ------------------------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------------------------
#
sub prepare {
	my ($ctlname) = @_;

	($ctlcount, $ctlindex) = (0, 0);

	my ($seqfile, $treefile, $tree, $ctl);

	# Read control file and skip some parameters
	open(F, "<", "$ctlname.ctl");
	while (my $line = <F>) {
		$line =~ s/\s*\*.*$//;
		next if ($line =~ m/^\s*$/);

		my ($key, $value) = split(/=/, $line);
		$key   =~ s/\s+//g;
		$value =~ s/^\s+//;
		$value =~ s/\s+$//;

		next if ($key =~ m/model|NSsites|fix_omega|omega/);
		next if ($params{$key});

		if ($key eq "seqfile") {
			$seqfile = $value;
			$value   = basename($value);
		}
		elsif ($key eq "treefile") {
			$treefile = $value;
			$value    = basename($value);
		}

		$ctl .= "$key = $value\n";
	}
	close(F);

	for my $key (sort keys %params) {
		$ctl .= "$key = $params{$key}\n";
	}

	if (!$seqfile || !-f $seqfile) {
		message("E", "Missing seqfile $seqfile!");
		return;
	}
	if (!$treefile || !-f $treefile) {
		message("E", "Missing treefile $treefile!");
		return;
	}

	$tree = join("", readFile($treefile));
	$tree =~ s/\s+//g;
	return if !($tree);

	$ctlcount = getRuns($tree);
	$ctlindex = grep {-f "$_/DONE"} <$ctlname-*>;

	if (!$ctlcount) {
		message("E", "It seems that there are no runs to be done!");
		return;
	}

	return (realpath($seqfile), realpath($treefile), $tree, $ctl);
}

#
# ------------------------------------------------------------------------------
# Loop over all *.ctl files and start runs
# ------------------------------------------------------------------------------
#
sub loop {
	cleanUndones(@ctlnames);

	for my $ctlname (@ctlnames) {
		my $resultfile = "$ctlname.result";
		if (-f $resultfile) {
			message("I", "The run for $ctlname is already done!");
			next;
		}

		message(">", "Take run $ctlname...");

		my @a = prepare($ctlname);
		next if (!@a);

		my ($seqfile, $treefile, $tree, $ctl) = @a;

		message("I", "There are $ctlindex of $ctlcount runs already done...");

		# Site specific
		if ($tests =~ m/1/) {
			my $c = extendCtl($ctl, "0", "1 2 7 8", "0");
			my $t = $tree;
			runCodeml($ctlname, $seqfile, $treefile, $c, $t, "10", 0);
		}

		# Branch site model with and without selection
		if ($tests =~ m/2/) {
			my $c = extendCtl($ctl, "2", "2", "0");
			my $t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "20", $treeno);
			}

			$c = extendCtl($ctl, "2", "2", "1");
			$t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "21", $treeno);
			}
		}

		# Branch model with and without selection
		if ($tests =~ m/3/) {
			my $c = extendCtl($ctl, "2", "0", "0");
			my $t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "30", $treeno);
			}

			$c = extendCtl($ctl, "2", "0", "1");
			$t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "31", $treeno);
			}
		}

		# Calculate hyphy results
		if ($tests =~ m/h/) {
			runHyphy($ctlname, $seqfile, $treefile);
		}
	}

	message(">", "Waiting for remaining runs...");
	sleep(2);
	while (1) {
		checkMaxtime();
		my @pids = getSubpids();
		last if (!@pids);
		sleep(10);
	}

	generate();

	printErrors();

	message(">", "Finished!");
}

#
# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
#
sub main {
	$SIG{INT}  = \&terminate;
	$SIG{TERM} = \&terminate;

	getParams();

	@ctlnames = map {$_ =~ s/\.ctl$//; $_} <*.ctl> if (!@ctlnames);

	if ($info) {

		#print "==> Running paPAML processes (pid: treefolder)\n";
		#my $proc = Proc::ProcessTable->new();
		#foreach my $p (@{$proc->table}) {
		#	if ($p->uid == $< && $p->cmndline =~ m/sh RUN paPAML/) {
		#		my @a = split(/\s+/, $p->cmndline);
		#		printf("%d: %s\n", $p->pid, $a[3]);
		#	}
		#}

		foreach my $ctlname (@ctlnames) {
			if (-f "$ctlname.result") {
				message("I", "Run $ctlname is finished!");
				next;
			}
			my @a = grep {$_ =~ m/treefile\s*=/} readFile("$ctlname.ctl");
			if (@a && $a[0] =~ m/=\s*([^\s]+)/) {
				my $n = getRuns(join("", readFile($1)));
				my $m = grep {-f "$_/DONE"} <$ctlname-*>;
				message("I", sprintf("Run $ctlname finished: %.1f%% (%d/%d)", $m / $n * 100.0, $m, $n));
			}
			else {
				message("E", "No treefile in run $ctlname found!");
			}
		}
		exit(0);
	}
	elsif ($clean) {
		cleanRuns(@ctlnames);
		exit(0);
	}

	my $codeml = which("codeml");
	if (!$codeml) {
		message("E", "There is no program called codeml in path!");
		exit(1);
	}

	my $hyphy = which("hyphy");
	if (!$hyphy) {
		message("E", "There is no program called hyphy in path!");
		exit(1);
	}

	if (!@ctlnames) {
		message("E", "There is no control file!");
		exit(1);
	}

	symlink($codeml, $codemlpgm);
	symlink($hyphy,  $hyphypgm);

	loop();

	cleanSymlinks();

	my $runtime = time - $starttime;
	if (-f $RUNTIMEFILE) {
		$runtime += join("", readFile($RUNTIMEFILE));
		unlink($RUNTIMEFILE);
	}

	message(">", sprintf("The total runtime was %.1f minutes", $runtime / 60));

	# Set exit code if at least one run and one (sub-)tree failed
	if (!$debug) {
		for my $ctlname (@ctlnames) {
			my @dirs = grep {-d $_} <$ctlname-*>;
			exit(1) if (@dirs);
		}
	}
}

main();
