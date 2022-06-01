#!/usr/bin/env perl

#
# ==============================================================================
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
use File::Path qw(rmtree);
use File::Which qw(which);
use Proc::ProcessTable;
use Statistics::Distributions;

# The tag/extension for the links of running codeml, hyphy by special name
my $tag = time;
my $codemlpgm = "codeml-$tag";
my $hyphypgm = "hyphy-$tag";

# Parameters
my $significance = 0.05;
my $para         = 0;
my $tests        = "123h";
my $debug        = 0;

my ($info, $clean);

my @ctlnames;

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
    paPAML.pl -i
    paPAML.pl -c

VERSION 1.22

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
    This program takes specified or all control files (*.ctl) of the
    actual folder and calls paml/codeml and hyphy (in paralell as
    background processes) for them.  There are several tests that are
    done, either all or selected ones:

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

    As already mentioned, a finished run in a subfolder will have a
    file called DONE - this is the indicator that the calculation for
    that specific run/tree is done.  On the other hand on an error
    there exists a file called ERROR.

    The program stays running, as long as your codeml or hyphy jobs
    are working.  If you press CTRL-C the program finnishes and all
    running codeml and hyphy runs are canceled.

    When all runs finnished succesfully (or if in the meantime all
    needed jobs of a control file are done) there will be a *.result
    file for every *.ctl file found and the generated subfolders will
    be removed (calling without the "-d" parameter!).

    Restarting paPAML.pl again with error runs will first remove the
    folders with ERROR files and rerun those runs again.

    An interrupted run can/will start from the point (subfolder) where
    it was stopped (it has no DONE file!) - as long as the
    intermediate data was not deleted.

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
# Clean all generated files for control files not needed any more
# ------------------------------------------------------------------------
#
sub cleanRuns {
	for my $ctlname (@_) {
		print "==> Clean $ctlname...\n";
		my @dirs = grep {-d $_} <$ctlname-*>;
		for my $dir (@dirs) {
			print "Clean run $dir...\n";
			rmtree($dir);
		}
	}
}

#
# ------------------------------------------------------------------------
# Clean error runs
# ------------------------------------------------------------------------
#
sub cleanErrors {
	for my $ctlname (@_) {
		my @dirs = grep {-d $_} <$ctlname-*>;
		for my $dir (@dirs) {
			if (-f "$dir/ERROR") {
				print "==> Clean error run $dir...\n";
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
			print "==> Clean symlink $pgm\n";
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
				print "[E] parameter -p needs an integer value!\n";
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
				print "[E] control file(s) " . join(",", @errors) . " do not exist!\n";
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
				print "[E] parameter -t needs value(s) of 1, 2, 3 and/or h!\n";
				exit(1);
			}
			if ($tests ne $s) {
				print "--> Parameter tests (-t) is changed to $tests!\n";
			}
		}
		elsif ($p eq "-s") {
			$significance = $ARGV[$i++ + 1];
			if (!($significance =~ m/^\d*\.?\d+$/)) {
				print "[E] parameter -s needs an float value!\n";
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
		print "[E] The parameters are incorrect!\n";
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
# Kills all codeml processes
# ------------------------------------------------------------------------
#
sub terminate {
	print "==> Terminate...\n";

	my @pids = getSubpids();
	for my $pid (@pids) {
		print "--> process $pid...\n";
		kill(15, $pid) if (kill(0, $pid));
	}

	sleep(2);

	@pids = getSubpids();
	for my $pid (@pids) {
		print "--> process $pid...\n";
		kill(9, $pid) if (kill(0, $pid));
	}

	cleanSymlinks();

	exit(1);
}

#
# ------------------------------------------------------------------------------
# Wait for the next free slot
# ------------------------------------------------------------------------------
#
sub dowait {
	while (1) {
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
				print "[E] There are errors in run $dir!\n";
				$error = 1;
			}
			my @a = qx(grep "is missing in the tree" $dir/codeml.log) if (-f "$dir/codeml.log");
			for my $s (@a) {
				chomp($s);
				print "[E] Error from codeml run $dir: $s\n";
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
# Prints the final sequnce
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

	# Loop over codons
	my $count = length($seq) / 3;
	for (my $i = 0 ; $i < $count ; $i++) {
		my $codon = substr($seq, $i * 3, 3);
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
	}

	open(FASTA, ">", "$ctlname.result.fa");
	printf FASTA ">T1_Bayes_1_2\n";
	for (my $i = 0 ; $i < $count ; $i++) {
		my $codon = substr($seq, $i * 3, 3);
		printf FASTA ("%s",	exists $bayes12->{$i} ? uc($codon) : lc($codon));
	}
	printf FASTA "\n";
	printf FASTA ">T1_Bayes_7_8\n";
	for (my $i = 0 ; $i < $count ; $i++) {
		my $codon = substr($seq, $i * 3, 3);
		printf FASTA ("%s",	exists $bayes78->{$i} ? uc($codon) : lc($codon));
	}
	printf FASTA "\n";
	printf FASTA ">T1_Bayes_1_2\n";
	for (my $i = 0 ; $i < $count ; $i++) {
		my $codon = substr($seq, $i * 3, 3);
		my $trees = substr($codondata->{$i}->{"2"}, 1);
		printf FASTA ("%s",	$trees ? uc($codon) : lc($codon));
	}
	printf FASTA "\n";
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

		print "==> Generate result for $ctlname...\n";

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

	if (!($tree =~ m/#1/)) {
		$tree =~ s/\((\w+)/\($1 #1/;
	}
	else {
		$tree =~ m/^(.*) #1(.*)$/;
		my ($head, $tail) = ($1, $2);
		return "" if (!($tail =~ m/\w/));
		if (!($tail =~ s/^\)/\) #1/)) {
			if (!($tail =~ s/,(\(*\w+)/,$1 #1/)) {
				$tail =~ s/\((\w+,\w+)\)/\($1\) #1/;
			}
		}
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
my $runcount = 0;

sub runCodeml {
	my ($ctlname, $seqfile, $treefile, $ctl, $tree, $treetype, $treeno) = @_;

	my $subdir = sprintf("%s-%02d-%05d", $ctlname, $treetype, $treeno);
	return if (-f "$subdir/DONE");

	generate() if ((++$runcount % 10) == 0);

	dowait();

	printf("==> Start %s-%02d-%05d\n", $ctlname, $treetype, $treeno);

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

	generate() if ((++$runcount % 10) == 0);

	dowait();

	printf("==> Start %s\n", $subdir);

	mkdir($subdir);
	chdir($subdir);

	symlink("../$alignfile", $alignfile);
	symlink("../$treefile",  $treefile);

	my $command = "../$hyphypgm fel --pvalue $significance --alignment $alignfile --tree $treefile >out";
	writeFile("RUN", "($command; if [ \$? -eq 0 ]; then touch DONE; else touch ERROR; fi)&");
	system("sh RUN paPAML $subdir");

	chdir("..");
}

#
# ------------------------------------------------------------------------------
# Loop over all *.ctl files and start runs
# ------------------------------------------------------------------------------
#
sub loop {
	for my $ctlname (@ctlnames) {
		my $resultfile = "$ctlname.result";
		next if (-f $resultfile);

		print "==> Run control file $ctlname...\n";

		my ($seqfile, $treefile);

		# Read control file and skip some parameters
		my $ctl;
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
			print "[E] Missing seqfile $seqfile!\n";
			next;
		}
		if (!$treefile || !-f $treefile) {
			print "[E] Missing treefile $treefile!\n";
			next;
		}

		# Read the tree
		open(F, "<", $treefile);
		my $tree = join("", <F>);
		close(F);
		$tree =~ s/\s+//g;

		next if (!$tree);

		if ($tests =~ m/1/) {

			# Site specific
			my $c = extendCtl($ctl, "0", "1 2 7 8", "0", "1");
			my $t = $tree;
			runCodeml($ctlname, $seqfile, $treefile, $c, $t, "10", 0);
		}

		if ($tests =~ m/2/) {

			# Branch site model with selection
			my $c = extendCtl($ctl, "2", "2", "0", "1");
			my $t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "20", $treeno);
			}

			# Branch site model without selection
			$c = extendCtl($ctl, "2", "2", "1", "1");
			$t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "21", $treeno);
			}
		}

		if ($tests =~ m/3/) {

			# Branch model with selection
			my $c = extendCtl($ctl, "2", "0", "0", "1");
			my $t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "30", $treeno);
			}

			# Branch model without selection
			$c = extendCtl($ctl, "2", "0", "1", "1");
			$t = $tree;
			for (my $treeno = 0 ; $t = mark($t) ; $treeno++) {
				runCodeml($ctlname, $seqfile, $treefile, $c, $t, "31", $treeno);
			}
		}

		if ($tests =~ m/h/) {

			# Calculate hyphy results
			runHyphy($ctlname, $seqfile, $treefile);
		}
	}

	print "==> Waiting to finish remaining runs...\n";
	while (1) {
		my @pids = getSubpids();
		last if (!@pids);
		sleep(10);
	}

	generate();

	printErrors();

	print "==> Finished!\n";
}

#
# ------------------------------------------------------------------------------
# Loop over all input files
# ------------------------------------------------------------------------------
#
sub main {
	$SIG{INT}  = \&terminate;
	$SIG{TERM} = \&terminate;

	getParams();

	@ctlnames = map {$_ =~ s/\.ctl$//; $_} <*.ctl> if (!@ctlnames);

	if ($info) {
		print "==> Running paPAML processes (pid: tree)\n";
		my $proc = Proc::ProcessTable->new();
		foreach my $p (@{$proc->table}) {
			if ($p->uid == $< && $p->cmndline =~ m/sh RUN paPAML/) {
				my @a = split(/\s+/, $p->cmndline);
				printf("%d: %s\n", $p->pid, $a[3]);
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
		print "[E] There is no program called codeml!\n";
		exit(1);
	}

	my $hyphy = which("hyphy");
	if (!$hyphy) {
		print "[E] There is no program called hyphy!\n";
		exit(1);
	}

	if (!@ctlnames) {
		print "[E] There is no control file!\n";
		exit(1);
	}

	symlink($codeml, $codemlpgm);
	symlink($hyphy, $hyphypgm);

	cleanErrors(@ctlnames);

	loop();

	cleanSymlinks();
}

main();
