#!/usr/bin/env perl
# doLABELlevel -  Version 1.0
# Performs LABEL formatting, prediction, and output.

# Copyright (C) 2012, Centers for Disease Control & Prevention
# Author: Samuel S. Shepard (vfn4@cdc.gov)
#
# GPL version 3
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

####################
# PROCESS parameters
use File::Basename;
use Getopt::Long;
use Cwd 'abs_path';
use File::Copy;

GetOptions(
		'retrain|R' => \$retrain,
		'score-field|F=i' => \$scoreField,
		'null-model|N=s' => \$nullModel, 
		'group|G=s' => \$group,
		'module|M=s' => \$module,
		'info|I=s' => \$infoFile,
		'label|L=s' => \$labelFile,
		'training|T=s' => \$trainingFile,
		'data-filter|D=s' => \$dataFilter,
		'shogun-path|S=s' => \$shogunPath,
		'best-score|B' => \$bestScoreMode,
		'custom-reverse|C=s' => \$custRevPath,
		'local|Q' => \$local
	  );

if ( scalar(@ARGV) < 2 ) {
	$message = "Usage:\n\tperl $0 <PPATH> <TNPATH> [FASTA ...] [OPTIONS]\n";
	$message .= "\t\t-F|--score-field <#>\t\tSpecified the scoring field to extract.\n";
	$message .= "\t\t-G|--group <STR>\t\tCurrent group.\n";
	$message .= "\t\t-M|--module <STR>\t\tModule (base) of the lineage assignment operation.\n";
	$message .= "\t\t-L|--label <FILE>\t\tFile for the training labels.\n";
	$message .= "\t\t-T|--training <FILE>\t\tFile for the training data.\n";
	$message .= "\t\t-I|--info <FILE>\t\tFile for the training information, including taxa.\n";
	$message .= "\t\t-D|--data-filter <FILE>\t\tFile for the data filter threshold (module-specific).\n";
	$message .= "\t\t-S|--shogun-path <PATH>\t\tPath to the shogun binary (OS-specific or linked).\n";
	$message .= "\t\t-R|--retrain\t\t\tSpecifies the classifier file does not exists.\n";
	$message .= "\t\t-C|--custom-reverse <path>\tUse custom reverse corrected viterbi model. Path contains equivalent pHMMs.\n";
	die($message."\n");
}

$ppath = $ARGV[0];
$tpath = $ARGV[1];
$PROGRAM = basename($0,'.pl');

if ( !defined($shogunPath) ) {
	$shogunPath = dirname(abs_path($0)). '/shogun';
}

if ( !defined($bestScoreMode) ) {
	$bestScoreMode = 0;
} else {
	$bestScoreMode = 1;
}


if ( !defined($module) ) {
	if ( abs_path($tpath) =~ /LABEL_RES\/training_data\/(.+?)(\/|\Z)/ ) {
		$module = $1;
	} else {
		die("$PROGRAM ERROR: Could not resolve module.\n");
	} 	
}

if ( !defined($group) ) {
	if ( abs_path($tpath) =~ /LABEL_RES\/training_data(\/.+)\/?/ ) {
		$predPath = $1;
	} else {
		die("$PROGRAM ERROR: Could not resolve prediction path/group.\n");
	} 	

} else {
	$predPath='/'.$module;
	if ( $group ne $module ) {
		$predPath .= "/$group";
	}
}	

if ( !defined($infoFile) ) {
	$infoFile = $tpath.'/info.dat';
}

if ( !defined($labelFile) ) {
	$labelFile = $tpath.'/labels.dat';
}

if ( !defined($classifier) ) {
	$classifier = $tpath.'/classifier.dat';
}

$local = defined($local) ? 1 : 0;

open( IN, '<', $labelFile ) or die("Cannot open training labels file ($labelFile).\n");
$line = <IN>; close(IN); chomp($line);
@labels = split(' ', $line);
foreach $label ( @labels ) {
	$classes{$label}=1;
}
$numClasses = scalar(keys(%classes));
close(IN);

if ( !defined($trainingFile) ) {
	$trainingFile = $tpath.'/training.dat';
}

if ( !defined($dataFilter ) ) {
	$dataFilter = 0;
} else {
	open(IN,'<',$dataFilter) or die("Cannot open $dataFilter.\n");
	$/ = "\n";
	$line = <IN>;
	chomp($line);
	($threshold,$maxLength) = split(' ', $line);
	close(IN);
	$dataFilter = 1;
	%filterByID = ();
}

# scoring field
if ( !defined($scoreField) ) {
	$scoreField = 3;
} else {
	if ( $scoreField > 0 ) {
		$scoreField--;
	} else {
		die("$PROGRAM ERROR: Score field must be a positive integer.\n");
	}
}

%nulls = %scores = %lengths = ();
# Process null model data
if ( defined($nullModel) ){
	$subtractNull = 1;
	open(NULL,'<',$nullModel) or die("$PROGRAM ERROR: cannot open $nullModel for reading.\n");
	$header = <NULL>;
	while ( $line = <NULL> ) {
		chomp($line);
		@values = split("\t", $line);
		$nulls{$values[0]} = $values[$scoreField];
	}
	close(NULL);
} else {
	$subtractNull = 0;
}

# Process taxa from training INFO.
$/ = "\n\n";
open( INFO, '<', $infoFile) or die("$0 ERROR: Cannot open $infoFile.\n");
$record = <INFO>;
close(INFO);

@lines = split("\n", $record);
%labelsByTaxa = %taxaByLabels = ();
@taxa = @labels = ();
for($i = 0;$i < scalar(@lines); $i++ ) {
	($label, $taxon) = split("\t", $lines[$i]);
	$taxon =~ tr/\/\\/-/;
	$labelsByTaxa{$taxon} = $label;
	$taxaByLabels{$label} = $taxon;
	$taxa[$i] = $taxon;
	$labels[$i] = $label;
}
$numLabels = $i;

if ( $numClasses != $numLabels ) {
	die("$PROGRAM ERROR: Number classes ($numClasses) in $labelFile different ($numLabels) from in .../info.dat).\n");
}

if ( defined($custRevPath) ) {
	$xrev = 1;
} else {
	$xrev = 0;
}

# Process pHMM output
$/ = "\n"; %bestLabelByID = (); $theScore = 0;
for($i = 0; $i < $numLabels; $i++ ) {
	$file = $ppath.'/'.$taxa[$i].'_hmm.tab';
       	open( IN, '<', $file ) or die("$PROGRAM ERROR: Cannot open $file.\n");
	if ( $xrev ) {
		@revValues = ();
		%revScores = ();
		$revFile = $custRevPath . '/' . $taxa[$i].'_hmm.tab';
		open(REV,'<',$revFile) or die("Cannot open $revFile\n");
		$revHeader = <REV>;
		while($revLine =<REV> ) {
			chomp($revLine);
			@revValues = split("\t",$revLine);
			$revScores{$revValues[0]} = $revValues[$scoreField];
		}	
		close(REV);
	}	

	$header = <IN>;
	$label = $labels[$i];
	while ( $line = <IN> ) {
		chomp($line);
		@values = split("\t", $line);
		$id = $values[0];
		$lengths{$id} = $values[1];
		if ( $subtractNull ) {
			$scores{$id}{$label} = sprintf("%.2f",$values[$scoreField] - $nulls{$id});
			$theScore = $values[2]; # not subtracted because neighborhood doesn't work
		} elsif ( $xrev ) {
			$theScore = $scores{$id}{$label} = sprintf("%.2f",$values[$scoreField] - $revScores{$id});
		} else {
			$scores{$id}{$label} = $values[$scoreField];
			$theScore = $values[$scoreField];
		}

		if ( $dataFilter || $bestScoreMode ) {
			if ( !defined($filterByID{$id}) ) {
				$filterByID{$id} = $theScore;
				if ( $bestScoreMode ) {
					$bestLabelByID{$id} = $label;
				}
			} elsif ( $theScore < $filterByID{$id} ) {
				$filterByID{$id} = $theScore;
				if ( $bestScoreMode ) {
					$bestLabelByID{$id} = $label;
				}
			}
		}
	}
	close(IN);
}
@sortedIDs = sort(keys(%scores)); 
$Ni = scalar(@sortedIDs);


# OUTPUT test data
$d = " "; 
open( DAT, '>', $ppath.'/test.dat') or die("$0 ERROR:Cannot open $ppath/test.dat for writing\n");
for( $i = 0; $i < $Ni; $i++ ) {
	print DAT $scores{$sortedIDs[$i]}{$labels[0]};
	foreach( $h = 1; $h < $numLabels; $h++ ) {
		print DAT $d,$scores{$sortedIDs[$i]}{$labels[$h]};
	}
	print DAT "\n";
}
close(DAT); 

if ( $local ) {
	copy($trainingFile,"$ppath/training.dat") or die "Copy failed: $!";
	copy($labelFile,"$ppath/labels.dat") or die "Copy failed: $!";
	copy($classifier,"$ppath/classifier.dat") or die "Copy failed: $!";
	$trainingFile = 'training.dat';
	$labelFile = 'labels.dat';
	$classifier = 'classifier.dat';
}

if ( ! $bestScoreMode ) {
	# TYPE of kernel
	if ( $numClasses > 2 ) {
		$svm_type = 'GMNPSVM';
		$kernel = 'POLY REAL 20 20 1 1';
	} else {
		$svm_type = 'LIBSVM';
		$kernel = 'POLY REAL 20 20 1 1';
	}

	# BEGIN shogun script
	$sg = '';
	$sg .=	"set_kernel $kernel\n";
	$sg .=	"set_features TRAIN $trainingFile\n";
	$sg .=	"set_labels TRAIN $labelFile\n";
	$sg .=	"new_classifier $svm_type\n";

	# BRANCH for re-training
	if ( defined($retrain) ) {
		$sg .=	"train_classifier\n";
		$sg .=	"save_classifier $classifier\n";
	} else {
		$sg .=	"load_classifier $classifier $svm_type\n";
	}

	# FINISH and execute shogun script
	$testResults = "test_results.dat.tmp";
	if ( $local ) {
		$sg .=	"set_features TEST test.dat\n";
		$sg .=	"$testResults = classify\n";
		$sg .=  "! pwd";
		$out = `cd "$ppath" && echo "$sg" | "$shogunPath"`;
	} else {
		$sg .=	"set_features TEST $ppath/test.dat\n";
		$sg .=	"$ppath/$testResults = classify\n";
		$out = `echo "$sg" | "$shogunPath"`;
	}

	# READ shogun results.
	$/ = "\n";
	open(IN, '<', "$ppath/$testResults") or die("$0 ERROR: result file \"$ppath/$testResults\" was not written by SHOGUN.\n");
	$line = <IN>; chomp($line);
	@results = split(/ /, $line);
	close(IN);

} else {
	for($i = 0; $i < $Ni;$i++) {
		$results[$i] = $bestLabelByID{$sortedIDs[$i]};
	}
}

open( TXT, '>', $ppath.'/LEVEL_result.txt' ) or die("$PROGRAM ERROR: Cannot open LEVEL_result.txt\n");
open( TAB, '>', $ppath.'/LEVEL_result.tab' ) or die("$PROGRAM ERROR: Cannot open LEVEL_result.tab\n");
open( TRACE, '>', $ppath.'/LEVEL_trace.tab' ) or die("$PROGRAM ERROR: Cannot open LEVEL_trace.tab\n");

# OUTPUT the results
%lineagesByID = ();
print TXT sprintf("%80s\t%20s\t%s\n", 'ID','CLASSIFICATION','PREDICTION PATH' );
print TRACE 'ID',"\t",'ANNOTATION',"\t",'PREDICTION';
for( $h = 0; $h < $numLabels; $h++ ) {
	print TRACE "\t",'HMM_',$taxa[$h];
}
print TRACE "\n";
for($i=0; $i < $Ni;$i++ ) {
	$id = $sortedIDs[$i];
	if ( $id =~ /{([^{}]*)}\s*$/ ) {
		$annot = $1;
	} else {
		$annot = 'UNKNOWN';
	}

	if ( $numClasses > 2 ) {
		$pred = $taxaByLabels{int($results[$i])};
	} else {
		$pred = $taxaByLabels{sgn($results[$i])};
	}

	if ( $dataFilter ) {
		# There is a bound to the normalization to ensure appending regions don't unfairly affect the stats.
		if ( $lengths{$id} > $maxLength ) {
			$filterByID{$id} /= $maxLength;
		} else {
			$filterByID{$id} /= $lengths{$id};
		}

		if ( $filterByID{$id} > $threshold ) {
			$pred = 'UNRECOGNIZABLE';
		}
	}

	print TRACE $id,"\t",$annot,"\t",$pred;
	for( $h = 0; $h < $numLabels; $h++ ) {
		print TRACE "\t",$scores{$id}{$labels[$h]};
	}
	print TRACE "\n";

	print TXT sprintf("%80s\t%20s\t%s\n",$id,$pred,$predPath);
	print TAB $id,"\t",$pred,"\n";
	$lineagesByID{$id} = $pred;
}
close(TXT);
close(TAB);
close(TRACE);

$/ = ">"; 
for($i = 2;$i < scalar(@ARGV); $i++ ) {
	open(SEQ, '<', $ARGV[$i]) or die("$0 ERROR: Cannot open $ARGV[$i].\n");
	open(OUT, '>', $ppath.'/'.$ARGV[$i].'.predictions') or die("$0 ERROR: Cannot open $ARGV[$i].predictions\n");
	while($record = <SEQ> ) {
		chomp($record);
		@lines = split(/\r\n|\n|\r/, $record);
		$id = $header = shift(@lines);
		$sequence = join('',@lines);

		if ( length($sequence) == 0 ) {
			next;
		}

		$header .= $delim.'{PRED:'.$lineagesByID{$id}.'}';
		print OUT '>',$header,"\n",$sequence,"\n";
	}
	close(SEQ);
	close(OUT);
}

# FNC - sign
sub sgn($) {
	my $a = shift;
	if ( $a > 0 ) {
		return 1;
	} elsif ( $a < 0 ) {
		return -1;
	} else {
		$a = substr($a,0,1);
		if ( $a eq '-' ) {
			return -1;
		} elsif ( $a =~ /\d/ ) {
			return 1;
		} else {
			print STDERR "$0 WARNING: unexpected shogun return class.";
			return 0;
		}
	}
}

# FNC - PARSE taxa into its short form.  Deprecated.
sub parseTaxa($) {
	my $taxon = $_[0];
	my @parts = ();
	
	$taxon =~ tr/\//-/;
	$taxon =~ tr/+//d;
	if ( length($taxon) > 10 ) {
		@parts = split('_',$taxon);
		if ( scalar(@parts) > 1 ) {
			$taxon = $parts[0].'-'.$parts[$#parts];
		} else {
			$taxon = substr($taxon,0,1).'-'.substr($taxon,-1);
		}
	}

	return $taxon;
}
