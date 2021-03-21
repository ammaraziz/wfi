#!/usr/bin/env perl
# doSVM -  Version 1.0
# Generates a shogun cmdline_static compatible script using input params.
# Outputs a revised result file using ID information.
#
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
# doSVM.py GROUP TEST RESULT TRAIN LABEL SHOGUN
# doSVM.py GROUP TEST RESULT TRAIN LABEL SHOGUN CLASSIFIER
if ( scalar(@ARGV) == 7 ) {
	$retrain = 1;
} elsif ( scalar(@ARGV) > 8 || scalar(@ARGV) < 7 ) {
	die("Usage:\n\tperl $0 <GROUP> <TEST.dat> <TEST_ID.dat> <OUTPUT> <TRAIN.dat> <LABEL.dat> <shogun> [CLASSIFIER]\n");
}


open( IN, '<', $ARGV[5] ) or die("Cannot open labels file ($ARGV[5]).\n");
$line = <IN>; close(IN); chomp($line);
@labels = split(' ', $line);
foreach $label ( @labels ) {
	$classes{$label}=1;
}
$numClasses = scalar(keys(%classes));
close(IN);

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
$sg .=	"set_features TRAIN $ARGV[4]\n";
$sg .=	"set_labels TRAIN $ARGV[5]\n";
$sg .=	"new_classifier $svm_type\n";

# BRANCH for re-training
if ( $retrain ) {
#	$classifier = $ARGV[0] . '_classifier.dat';
	$classifier = 'classifier.dat';
	$sg .=	"train_classifier\n";
	$sg .=	"save_classifier $classifier\n";
} else {
	$sg .=	"load_classifier $ARGV[7] $svm_type\n";
}

# FINISH and execute shogun script
$sg .=	"set_features TEST $ARGV[1]\n";
$outfile = $ARGV[3] . '.tmp';
$sg .=	"$outfile = classify\n";
$out = `echo "$sg" | $ARGV[6]`;

#print STDERR $out,"\n";

# READ shogun results.
open(IN, '<', "$outfile") or die("$0 ERROR: result file $outfile was not written by SHOGUN.\n");
$/ = "\n";
$line = <IN>;
chomp($line);
@results = split(/ /, $line);
close(IN);


# READ the ID file.
@IDs = ();
open(ID_FILE, '<', $ARGV[2] ) or die("$0 ERROR: Cannot open $ARGV[2] for reading.\n"); 
while( $line = <ID_FILE> ) {
	chomp($line);
	($number, $ID) = split(/\t/, $line);
	$IDs[$number-1] = $ID;	
}
close(ID_FILE);


# WRITE revised results using ID information.
open(OUT, '>', $ARGV[3] ) or die("$0 ERROR: Cannot open $ARGV[3] for writing.\n");
$N = scalar(@results);
for($i = 0; $i < $N; $i++ ) {
	if ( defined($IDs[$i]) ) {
		if ( $numClasses > 2 ) {	
			print OUT sprintf("%d\t%d\t%s\t%s\n",($i+1),int($results[$i]),$ARGV[0],$IDs[$i]);
		} else {
			print OUT sprintf("%d\t%d\t%s\t%s\n",($i+1),sgn($results[$i]),$ARGV[0],$IDs[$i]);
		}
	} else {
		die("$0 ERROR: Result $i without ID.\n");
	}
}
close(OUT);
####################

sub sgn($) {
	my $a = shift;
	if ( $a > 0 ) {
		return 1;
	} elsif ( $a < 0 ) {
		return -1;
	} else {
		return 0;
	}
}
