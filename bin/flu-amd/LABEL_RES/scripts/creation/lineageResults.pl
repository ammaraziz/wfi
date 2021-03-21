#!/usr/bin/env perl
# lineageResults -  Version 1.0
# Conveys the prediction results for the current prediciton level.
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
if ( scalar(@ARGV) < 5 ) {
	die("Usage:\n\tperl $0 <TRAINING_INFO> <GROUP> <prediction_results> <out_prefix> <prot> [<input1.fasta>... <inputN.fasta>]\n");
}
open( OUT, '>', $ARGV[3].'_result.txt' ) or die("$0 ERROR: Cannot open $ARGV[3]_result.txt\n");
open( TAB, '>', $ARGV[3].'_result.tab' ) or die("$0 ERROR: Cannot open $ARGV[3]_result.tab\n");

$group = $ARGV[1];
$delim = '|';

# GATHER lineage/clade labels.
open( LIN, '<', $ARGV[0] ) or die("$0 ERROR: Cannot open $ARGV[0].\n");
$/ = "\n\n"; %lineages = (); %taxa = (); #$gene = '';
$record = <LIN>;
chomp($record);
close(LIN);
@lines = split(/\r\n|\n|\r/, $record);
foreach $line ( @lines ) {
	($label, $lineage) = split("\t", $line);
	$lineages{$group}{$label} = $lineage;
	$taxa{$lineage} = 1;
}


open( RES, '<', $ARGV[2] ) or die("$0 ERROR: Could not open $ARGV[2].\n");
$/ = "\n"; %labelBySequence = ();
while( $line = <RES> ) {
	chomp($line);
	($number, $label, $group, $id) = split("\t", $line);


	$labelBySequence{$id}{$group} = $label;
}
close(RES);
@sorted = sort(keys(%labelBySequence));

# OUTPUT the results
%lineagesByID = ();
print OUT sprintf("%80s\t%20s\t%s\n", 'VIRUS STRAIN','CLADE','PREDICTION PATH' );
foreach $id ( @sorted ) {
	$lineage = $lineages{$group}{$labelBySequence{$id}{$group}};
	$predictionPath='/'.$ARGV[4];
	if ( $group ne $ARGV[4] ) {
		$predictionPath .= "/$group";
	}
	print OUT sprintf("%80s\t%20s\t%s\n",$id,$lineage,$predictionPath);
	print TAB $id,"\t",$lineage,"\n";
	$lineagesByID{$id} = $lineage;
}
close(OUT);
close(TAB);

$/ = ">"; 
for($i = 5;$i < scalar(@ARGV); $i++ ) {
	open(SEQ, '<', $ARGV[$i]) or die("$0 ERROR: Cannot open $ARGV[$i].\n");
	open(OUT, '>', $ARGV[$i].'.predictions') or die("$0 ERROR: Cannot open $ARGV[$i].predictions\n");
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
