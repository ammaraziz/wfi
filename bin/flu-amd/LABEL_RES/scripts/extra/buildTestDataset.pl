#!/usr/bin/env perl
# buildTestDataset -  Version 1.1
# Formats data into a matrix for the SVM testing phase.

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

GetOptions(	'shorten-taxa|S' => \$SHORTEN, 'transpose|T' => \$transpose, 'score-field|F=i' => \$scoreField, 'null-model|N=s' => \$nullModel );
if ( scalar(@ARGV) < 3 ) {
	$message = "Usage:\n\tperl $0 <info_file> <output_prefix> <scores.tab> ... <scoresN.tab>\n";
	$message .= "\t\t-S|--shorten-taxa\tShorten taxa name if over 10 characters and not a c-* or outlier pattern.\n";
	$message .= "\t\t-F|--score-field <#>\tScore field to build matrix.\n";
	die($message);
}
$PROGRAM = basename($0,'.pl');

if ( !defined($SHORTEN) ) {
	$SHORTEN = 0;
} else {
	$SHORTEN = 1;
}

if ( !defined($scoreField) ) {
	$scoreField = 3;
} else {
	if ( $scoreField > 0 ) {
		$scoreField--;
	} else {
		die("$PROGRAM ERROR: Score field must be a positive integer.\n");
	}
}
$info = shift(@ARGV);
$prefix = shift(@ARGV);

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


# GATHER profile scoring data
open( DAT, '>', $prefix.'.dat') or die("$0 ERROR:Cannot open ${prefix}.dat\n");
open( INFO, '<', $info) or die("$0 ERROR: Cannot open $info.\n");
$/ = "\n\n";
$header = <INFO>;
close(INFO);
$i = 0;
@lines = split("\n", $header);
%scores = ();
%labelsByTaxa =();
@taxa = @labels = ();
foreach $line ( @lines ) {
	($label, $taxon) = split("\t", $line);
	$taxon =~ tr/\/\\/-/;
	$labelsByTaxa{$taxon} = $label;
	$taxa[$i] = $taxon;
	$labels[$i] = $label;
	$i++;
}

$/ = "\n"; $numberProfiles = 0;
foreach $file ( @ARGV ) {
	#process label
       	open( IN, '<', $file ) or die("$0 ERROR: Cannot open $file.\n");
	$header = <IN>;
	$file = basename($file);
	if ( $file =~ /(.+?)_hmm.*?\.tab/ ) {
		$taxon = $1;
	} else {
		die("$0 ERROR: Malformed file name ($file).\n");
	}	
	$label = $labelsByTaxa{$taxon};
	while ( $line = <IN> ) {
		chomp($line);
		@values = split("\t", $line);
		$id = $values[0];
		$lengths{$id} = $values[1];
		if ( $subtractNull ) {
			$scores{$id}{$label} = sprintf("%.2f",$values[$scoreField] - $nulls{$id});
		} else {
			$scores{$id}{$label} = $values[$scoreField];
		}
	}
	close(IN);
	$numberProfiles++;
}

$numTrainingProfiles = scalar(@labels);
if ( $numTrainingProfiles != $numberProfiles ) {
	die("$0 ERROR: Number of training profiles ($numTrainingProfiles) different from test data profiles ($numberProfiles).\n");
}
$d = " "; 
@sortedIDs = sort(keys(%scores)); 
$Ni = scalar(@sortedIDs);

if ( defined($transpose) ) {
	for( $i = 0; $i < $Ni; $i++ ) {
		print DAT $scores{$sortedIDs[$i]}{$labels[0]};
		foreach( $h = 1; $h < $numberProfiles; $h++ ) {
			print DAT $d,$scores{$sortedIDs[$i]}{$labels[$h]};
		}
		print DAT "\n";
	}
	for( $h = 0; $h < $numberProfiles; $h++ ) {
		print INFO $labels[$h],"\t",$taxa[$h],"\n";
	}
} else {
	for( $h = 0; $h < $numberProfiles; $h++ ) {
		print DAT $scores{$sortedIDs[0]}{$labels[$h]};
		for( $i = 1; $i < $Ni; $i++ ) {
			print DAT $d,$scores{$sortedIDs[$i]}{$labels[$h]};
		}
		print DAT "\n";
	}
}
close(DAT); 


# OUT data for the SVM and trace matrix
open( ID, '>', $prefix.'_IDs.dat') or die("$0 ERROR: Cannot open ${prefix}_IDs.dat\n");
open( TSV, '>', $prefix.'.tsv') or die("$0 ERROR: Cannot open ${prefix}.tsv\n");
$k = 1;

print TSV 'ID',"\t",'ANNOTATION';
for( $h = 0; $h < $numberProfiles; $h++ ) {
	print TSV "\t",'HMM_',$taxa[$h];
}
print TSV "\n";

foreach $id ( @sortedIDs ) {
	print ID $k,"\t",$id,"\n";
	if ( $id =~ /{(.+?)}/ ) {
		$annot = $1;
	} else {
		$annot = 'Unknown';
	}
	print TSV $id,"\t",$annot;
	for( $h = 0; $h < $numberProfiles; $h++ ) {
		print TSV "\t",$scores{$id}{$labels[$h]};
	}
	print TSV "\n";

	$k++;
}
close(ID);
close(TSV);
