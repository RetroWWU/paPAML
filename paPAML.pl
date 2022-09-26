#!/usr/bin/env perl

#
# ==============================================================================
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

use File::Basename;
use File::Path  qw(rmtree);
use File::Which qw(which);
use Proc::ProcessTable;
use Statistics::Distributions;

my $RUNTIMEFILE = "runtime";

# The tag/extension for the links of running codeml, hyphy by special name
my $tag       = time;
my $codemlpgm = "codeml-$tag";
my $hyphypgm  = "hyphy-$tag";

# Parameters
my $significance = 0.05;
my $para         = 0;
my $tests        = "123h";
my $debug        = 0;

my ($info, $clean);

my @ctlnames;

# The starttime of this run and the last time check and the interval.
# The lasttimecheck is put in "future" to have at least some time left
# for a prognose
my ($starttime, $lasttimecheck, $lasttimeinterval) = (time, time + 60, 60);

# The actual (total) number of runs
my $runindex = 0;

# The total runs for a tree / ctl file
my ($ctlcount, $ctlindex);

# Translate codons in amino acids
my %aacode = (
	TTT => "F",
	TTC => "F",
	TTA => "L",
	TTG => "L",
	TCT => "S",
	TCC => "S",
	TCA => "S",
	TCG => "S",
	TAT => "Y",
	TAC => "Y",
	TAA => "*",
	TAG => "*",
	TGT => "C",
	TGC => "C",
	TGA => "*",
	TGG => "W",
	CTT => "L",
	CTC => "L",
	CTA => "L",
	CTG => "L",
	CCT => "P",
	CCC => "P",
	CCA => "P",
	CCG => "P",
	CAT => "H",
	CAC => "H",
	CAA => "Q",
	CAG => "Q",
	CGT => "R",
	CGC => "R",
	CGA => "R",
	CGG => "R",
	ATT => "I",
	ATC => "I",
	ATA => "I",
	ATG => "M",
	ACT => "T",
	ACC => "T",
	ACA => "T",
	ACG => "T",
	AAT => "N",
	AAC => "N",
	AAA => "K",
	AAG => "K",
	AGT => "S",
	AGC => "S",
	AGA => "R",
	AGG => "R",
	GTT => "V",
	GTC => "V",
	GTA => "V",
	GTG => "V",
	GCT => "A",
	GCC => "A",
	GCA => "A",
	GCG => "A",
	GAT => "D",
	GAC => "D",
	GAA => "E",
	GAG => "E",
	GGT => "G",
	GGC => "G",
	GGA => "G",
	GGG => "G",
);

# The default parameters.  They can be changed over command line
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
	"fix_alpha"    => "1",
	"fix_blength"  => "-1",
	"fix_rho"      => "1",
	"getSE"        => "0",
	"icode"        => "0",
	"method"       => "0",
	"ndata"        => "1",
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
    paPAML.pl -p runs [-f controlfiles] [-t tests] [-s significance] [-d] {codemlparams}
    paPAML.pl -i [-f controlfiles]
    paPAML.pl -c

VERSION 2.0

WHERE
    runs         - the number of parallel runs
    controlfiles - a list of control files.  It is assumed they are named
                   with a suffix "".ctl"!  If not given all files with
                   that suffix are taken to be calculated
    tests        - the used tests (1, 2, 3 or h) to run the data against.
                   They can be written like "1" or "12" . The order does
                   not matter.
                   (default: $tests)
    significance - the maximum p value to print marked trees.  Used for
                   printing bayes values and in hyphy call
                   (default: $significance)
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
    needed jobs of a control file are done) there will be a *.result
    and a *.result.fa file for every *.ctl file found and the
    generated subfolders will be removed.  You can skip the deletion
    of the generated subfolders, when you use the -d (debug)
    parameter.

    Restarting paPAML.pl again after a termination or with error runs
    will first remove the subfolders without a DONE file and rerun
    those runs again.  Meaning: you can restart as many times as long
    as not all runs are done.

    If a *.result file exists for a *.ctl file, the run(s) will be
    skipped if it is started again.
EOF
	exit(0);
}

if (!@ARGV) {
	usage();
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
	printf("[%1s %s] %s\n", $type ? $type : "-", getDate(), $message);
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
			$para = $ARGV[$i++ + 1];
			if (!($para =~ m/^\d+$/)) {
				message("E", "Parameter -p needs an integer value!");
				exit(1);
			}
		}
		elsif ($p eq "-f") {
			my $s      = $ARGV[$i++ + 1];
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
			my $s = $ARGV[$i++ + 1];
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
				print "--> Parameter tests (-t) is changed to $tests!\n";
			}
		}
		elsif ($p eq "-s") {
			$significance = $ARGV[$i++ + 1];
			if (!($significance =~ m/^\d*\.?\d+$/)) {
				message("E", "Parameter -s needs a float value!");
				exit(1);
			}
		}
		elsif ($p eq "-d") {
			$debug = 1;
		}
		else {
			if ($p =~ m/^-/) {
				my $value = $ARGV[$i++ + 1];
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
			if ($p->cmndline =~ m/$codemlpgm|$hyphypgm/) {
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

	sleep(2);

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
# Prints codeml data
# ------------------------------------------------------------------------------
#
sub generateCodeml {
	my ($ctlname) = @_;

	my ($bayes12, $bayes78, $codondata) = ({}, {}, {});

	print RESULT qq(
# -----------------------------------------------------------------------------
# Codeml test results
# -----------------------------------------------------------------------------

# Reference
#
# (tree);
# Model_n    np lnl dNdS
# Tree_n_O   np lnl
# Tree_n_M   np lnl
# Test_1_2   abs(2*(lnl?Model_2-lnl?Model_1)) abs(np?Model_2-np?Model_1) pvalue
# Test_7_8   abs(2*(lnl?Model_8-lnl?Model_2)) abs(np?Model_8-np?Model_7) pvalue
# Test_n     abs(2*(lnl?Tree_n_O-lnl?Tree_n_M)) abs(np?Tree_n_O-np?Tree_n_M) pvalue
# Bayes_n
# nr aminoacid probability [postmean += stddev]
#
# n - number
# ? - value referenced to the following parameter

);

	# Calculate model 1
	if ($tests =~ m/1/) {
		print RESULT qq(# Test 1 - site specific\n\n);

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

		my @ms = (1, 2, 7, 8);
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
	for my $t ("2", "3") {
		next if (!($tests =~ m/$t/));

		print RESULT ($t == 2 ? "\n# Test 2 - branch-site specific" : "\n# Test 3 - branch specific"), "\n\n";

		# Get folders for "without or with" fixomega
		my @dirs0 = <$ctlname-${t}0-*>;
		my @dirs1 = <$ctlname-${t}1-*>;

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

			my $dltr = abs(2 * ($lnl1 - $lnl0));
			my $dn   = abs($np1 - $np0);
			my $p    = Statistics::Distributions::chisqrprob($dn, $dltr);

			if ($p < $significance) {
				printTree($dirs0[$treeno], $treeno);

				printf RESULT ("%s\t%d\t%s\n", "Tree_" . ($treeno + 1) . "_O", $np0, $lnl0);
				if (@b0) {
					printf RESULT ("Bayes_%d_O\n", $treeno + 1);
					print RESULT join("\n", map {sprintf(("%d\t%s\t%s", $_->[0], $_->[1], $_->[2]))} @b0), "\n";
					map {$codondata->{$_->[0] - 1}->{"2"} .= sprintf(",Tree_%d:%0.3f", ($treeno + 1), 1 - $_->[2])} @b0;
				}

				printf RESULT ("%s\t%d\t%s\n", "Tree_" . ($treeno + 1) . "_M", $np1, $lnl1);
				if (@b1) {
					printf RESULT ("Bayes_%d_M\n", $treeno + 1);
					print RESULT join("\n", map {sprintf(("%d\t%s\t%s", $_->[0], $_->[1], $_->[2]))} @b1), "\n";
					map {$codondata->{$_->[0] - 1}->{"2"} .= sprintf(",Tree_%d:%0.3f", ($treeno + 1), 1 - $_->[2])} @b1;
				}

				printf RESULT ("%s\t%f\t%d\t%f\n", "Test_" . ($treeno + 1), $dltr, $dn, $p);
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

	print RESULT qq(
# -----------------------------------------------------------------------------
# Test 4 - Hyphy FEL
# -----------------------------------------------------------------------------

# Reference
#
# codon-number "+"|"-" pvalue

);

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

	print RESULT qq(
# -----------------------------------------------------------------------------
# Sequence overview of site specific results
# -----------------------------------------------------------------------------

# Reference
#
# Codon_number\tCodon\tAmino_acid\tT1_Bayes_1_2\tT1_Bayes_7_8\tT2_Bayes\tHyphy_negative\tHyphy_positive

);

	my ($b12, $b78, $b, $hn, $hp);
	my $count = length($seq) / 3;

	# Loop over codons
	for (my $i = 0 ; $i < $count ; $i++) {
		my $codon = substr($seq, $i * 3, 3);
		my ($cu, $cl) = (uc($codon), lc($codon));
		my $cd    = $codondata->{$i};
		my $trees = substr($cd->{"2"}, 1);

		printf RESULT (
			"%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			$i + 1,
			$codon,
			$aacode{$codon},
			exists $bayes12->{$i} ? sprintf("%0.4f", $bayes12->{$i}) : ".",
			exists $bayes78->{$i} ? sprintf("%0.4f", $bayes78->{$i}) : ".",
			$trees                ? $trees                           : ".",
			exists $cd->{"h-"}    ? sprintf("%0.4f", $cd->{"h-"})    : ".",
			exists $cd->{"h+"}    ? sprintf("%0.4f", $cd->{"h+"})    : "."
		);

		$b12 .= exists $bayes12->{$i} ? $cu : $cl;
		$b78 .= exists $bayes78->{$i} ? $cu : $cl;
		$b   .= $trees                ? $cu : $cl;
		$hn  .= exists $cd->{"h-"}    ? $cu : $cl;
		$hp  .= exists $cd->{"h+"}    ? $cu : $cl;
	}

	open(FASTA, ">", "$ctlname.result.fa");
	print FASTA qq(>T1_Bayes_1_2\n$b12\n);
	print FASTA qq(>T1_Bayes_7_8\n$b78\n);
	print FASTA qq(>T2_Bayes\n$b\n);
	print FASTA qq(>Hyphy_negative\n$hn\n);
	print FASTA qq(>Hyphy_positive\n$hp\n);
	close(FASTA);
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
		printf RESULT (
			"# %04d-%02d-%02d %02d:%02d: Results for %s.ctl with tests %s and significance %f\n",
			$year + 1900,
			$mon + 1, $mday, $hour, $min, $ctlname, $tests, $significance
		);

		my ($bayes12, $bayes78, $codondata) = generateCodeml($ctlname);
		generateHyphy($ctlname, $codondata);
		generateSequence($ctlname, $bayes12, $bayes78, $codondata);

		close(RESULT);

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
		$tree =~ s/(\w+|\))([\d\.\:\[\]]+)?/$1$2 #1/;
	}
	else {
		$tree =~ m/^(.*) #1(.*)$/;
		my ($head, $tail) = ($1, $2);
		return "" if (!($tail =~ m/[\w\:\.]/));
		$tail =~ s/(\w+|\))([\d\.\:\[\]]+)?/$1$2 #1/;
		$tree = "$head$tail";
	}

	return $tree;
}

#
# ------------------------------------------------------------------------------
# Extend the control file by parameters
# ------------------------------------------------------------------------------
#
sub extendCtl {
	my ($ctl, $model, $nssites, $fixomega, $omega) = @_;
	my @a = grep {!($_ =~ m/(\s+|)(model|NSsites|fix_omega|omega)\s*=/)} split(/\n/, $ctl);
	return "model = $model\nNSsites = $nssites\nfix_omega = $fixomega\nomega = $omega\n" . join("\n", @a);
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

	symlink("../$seqfile", $seqfile);
	writeFile("codeml.ctl", $ctl);
	writeFile($treefile,    " 1\n$tree\n");

	writeFile("RUN", "(../$codemlpgm 2>&1 >codeml.log; if [ \$? -eq 0 ]; then touch DONE; else touch ERROR; fi)&");
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

	symlink("../$alignfile", $alignfile);
	symlink("../$treefile",  $treefile);

	my $command = "../$hyphypgm fel --pvalue $significance --alignment $alignfile --tree $treefile >out 2>/dev/null";
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
		}
		elsif ($key eq "treefile") {
			$treefile = $value;
		}

		$ctl .= "$key = $value\n";
	}
	close(F);

	for my $key (sort keys %params) {
		$ctl .= "$key = $params{$key}\n";
	}

	if (!$seqfile || !-f $seqfile) {
		messaeg("E", "Missing seqfile $seqfile!");
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

	return ($seqfile, $treefile, $tree, $ctl);
}

#
# ------------------------------------------------------------------------------
# Loop over all *.ctl files and start runs
# ------------------------------------------------------------------------------
#
sub loop {
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
			my $c = extendCtl($ctl, "0", "1 2 7 8", "0", "1");
			my $t = $tree;
			runCodeml($ctlname, $seqfile, $treefile, $c, $t, "10", 0);
		}

		# Branch site model with and without selection
		if ($tests =~ m/2/) {
			my $c = extendCtl($ctl, "2", "2", "0", "1");
			my $t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "20", $treeno);
			}

			$c = extendCtl($ctl, "2", "2", "1", "1");
			$t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "21", $treeno);
			}
		}

		# Branch model with and without selection
		if ($tests =~ m/3/) {
			my $c = extendCtl($ctl, "2", "0", "0", "1");
			my $t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "30", $treeno);
			}

			$c = extendCtl($ctl, "2", "0", "1", "1");
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

	message(">", "Waiting to finish remaining runs...");
	while (1) {
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

		message(">", "Runs");
		foreach my $ctlname (@ctlnames) {
			if (-f "$ctlname.result") {
				message("I", "Run $ctlname is finished!");
				next;
			}
			my @a = grep {$_ =~ m/treefile\s*=/} readFile("$ctlname.ctl");
			if (@a && $a[0] =~ m/=\s*([^\s]+)/) {
				my $runs  = getRuns(join("", readFile($1)));
				my @dones = grep {-f "$_/DONE"} <$ctlname-*>;
				message("I", sprintf("Run $ctlname finished: %.1f%%", @dones / $runs * 100.0));
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

	cleanUndones(@ctlnames);

	loop();

	cleanSymlinks();

	my $runtime = time - $starttime;
	if (-f $RUNTIMEFILE) {
		$runtime += join("", readFile($RUNTIMEFILE));
		unlink($RUNTIMEFILE);
	}

	message(">", sprintf("The total runtime was %.1f minutes", $runtime / 60));
}

main();
