#!/usr/bin/env perl
# buildSvmDataset -  Version 1.1
# Formats data into matrices for the SVM training.

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

GetOptions(	'shorten-taxa|S' => \$SHORTEN, 'transpose|T' => \$transpose, 'null-model|N=s' => \$nullModel );
if ( scalar(@ARGV) < 2 ) {
	$message = "Usage:\n\tperl $0 [OPTIONS] <output_dir> <scores.tab> ... <scoresN.tab>\n";
	$message .= "\t\t-S|--shorten-taxa\tShorten taxa name if over 10 characters and not a c-* or outlier pattern.\n";
	die($message);
}

if ( !defined($SHORTEN) ) {
	$SHORTEN = 0;
} else {
	$SHORTEN = 1;
}

$scoreField=3;
if ( defined($nullModel) ){
	$subtractNull = 1;
	open(NULL,'<',$nullModel) or die("`basename $0` ERROR: cannot open $nullModel for reading.\n");
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

# PROCESS output files
$prefix = shift(@ARGV);
open( LAB, '>', $prefix.'/labels.dat') or die("$0 ERROR: Cannot open labels.dat\n");
open( DAT, '>', $prefix.'/training.dat') or die("$0 ERROR: Cannot open training.dat\n");
open( INFO, '>', $prefix.'/info.dat') or die("$0 ERROR: Cannot open info.dat\n");

# GATHER data
$/ = "\n"; %scores = (); %lengths = ();
%labelsByTaxa = ();
%taxaByLabels = ();
$taxon = '';
foreach $file ( @ARGV ) {
	$file_taxon = basename($file);
	if ( $file_taxon =~ /^(.+?)_hmm/ ) {
		$file_taxon = $1;
	} else {
		die("Malformed file taxa ($file_taxon).\n");
	}

	# PROCESS labels from first appearance
       	open( IN, '<', $file ) or die("$0 ERROR: Cannot open $file.\n");
	$header = <IN>;
	@lines = <IN>;
	chomp(@lines);
	close(IN);

	$found = 0;
	foreach $line ( @lines ) {
		@values = split("\t", $line);
		$id = $values[0];
		if ( $file =~ m/outlier/g ) {
			$taxon = 'outlier';
		} elsif ( $id =~ m/\{(.+)\}$/ ) {
			$taxon = $1;
		} else {
			$taxon = 'outlier';
		}

		if ( parseTaxa($taxon) eq $file_taxon) {
			$found = 1;
			last;
		}
	}

	if ( !$found ) {
		print STDERR "WARNING! Setting profile to $file_taxon. May mismatch sequences.\n";
		$taxon = $file_taxon;
	}

	# FINISH remaining records
	foreach $line ( @lines ) {
		@values = split("\t", $line);
		$id = $values[0];
		if ( $subtractNull ) {
			$scores{$id}{$taxon} = sprintf("%.2f",$values[$scoreField] - $nulls{$id});
		} else {
			$scores{$id}{$taxon} = $values[$scoreField];
		}
		$lengths{$id} = $values[1];
	}
	$labelsByTaxa{$taxon} = -2;
}

# OUTPUT files for SVM test phase
$d = " "; 
@sortedIDs = sort(keys(%scores));
@taxa = sort(keys(%labelsByTaxa));
$Nt = scalar(@taxa);
if ( $Nt > 2 ) {
	for($label = 0; $label < $Nt; $label++ ) {
		$labelsByTaxa{$taxa[$label]} = $label;
		$taxaByLabels{$label} = $taxa[$label];
	}
} else {
	$labelsByTaxa{$taxa[0]} = -1;
	$taxaByLabels{'-1'} = $taxa[0];
	$labelsByTaxa{$taxa[1]} = 1;
	$taxaByLabels{'1'} = $taxa[1];
}
@labels = sort {$a <=> $b} keys(%taxaByLabels);

if ( defined($transpose) ) {
	for( $i = 0; $i < scalar(@sortedIDs); $i++ ) {
		print DAT $scores{$sortedIDs[$i]}{$taxa[0]};
		for( $h = 1; $h < $Nt; $h++ ) {
			print DAT $d,$scores{$sortedIDs[$i]}{$taxa[$h]};
		}
		print DAT "\n";
	}
	for( $h = 0; $h < $Nt; $h++ ) {
		print INFO $labels[$h],"\t",$taxaByLabels{$labels[$h]},"\n";
	}
} else {
	for( $h = 0; $h < $Nt; $h++ ) {
		print DAT $scores{$sortedIDs[0]}{$taxa[0]};
		for( $i = 1; $i < scalar(@sortedIDs); $i++ ) {
			print DAT $d,$scores{$sortedIDs[$i]}{$taxa[$h]};
		}
		print DAT "\n";
		print INFO $labels[$h],"\t",$taxaByLabels{$labels[$h]},"\n";
	}
}
close(DAT); 


print INFO "\n";
$i = 1;
foreach $id ( @sortedIDs ) {
	print INFO $i,"\t",$id,"\n";
	$i++;
}
close(INFO);

for( $i = 0; $i < scalar(@sortedIDs); $i++ ) {
	$id = $sortedIDs[$i];
	if ( $id =~ m/\{(.+)\}$/ ) {
		$label = $labelsByTaxa{$1};
	} else {
		$label = $labelsByTaxa{'outlier'};
	}

	if ( $label eq '' ) { die("$0 ERROR: Missing data.  See $id.\n"); }
	if ( $i < $#sortedIDs ) {
		print LAB $label,$d;
	} else {
		print LAB $label,"\n";
	}
}
close(LAB);

# FNC - PARSE taxa into its short form.  Deprecated.
sub parseTaxa($) {
	my $taxon = $_[0];
	my @parts = ();
	
	$taxon =~ tr/\//-/;
	$taxon =~ tr/+//d;
	if ( $SHORTEN && length($taxon) > 10 && $taxon !~ /outlier/ && $taxon !~ /^c-/ ) {
		@parts = split('_',$taxon);
		if ( scalar(@parts) > 1 ) {
			$taxon = $parts[0].'-'.$parts[$#parts];
		} else {
			$taxon = substr($taxon,0,1).'-'.substr($taxon,-1);
		}
	}

	return $taxon;
}
####################
