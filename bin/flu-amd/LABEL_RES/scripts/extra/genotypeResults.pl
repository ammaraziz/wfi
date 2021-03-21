#!/usr/bin/env perl
# genotypeResults - Version Alpha
# Used in experimental gene ensemble genotyping system for final genotyping results.
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
if ( scalar(@ARGV) < 2 ) {
	die("Usage:\n\tperl $0 <genotype.table> <input1.tab> [... <inputN.tab>]\n");
}

$/ = "\n"; %genotypes = ();
open ( GEN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0].\n");
$line = <GEN>; chomp($line);
@fieldNames = split("\t", $line );
while ( $line = <GEN> ) {
	chomp($line);
	@lineages = split("\t", $line );
	for($i = 1; $i < scalar(@lineages); $i++ ) {
		$genotypes{uc($fieldNames[$i])}{uc($lineages[$i])}{$lineages[0]}=1;
	}
}
close(GEN);


for($i = 1;$i < scalar(@ARGV); $i++ ) {
	open( RES, '<', $ARGV[$i] ) or die("Cannot open $ARGV[$i].\n");
	while( $record = <RES> ) {
		chomp($record);
		($ID, $lineage) = split("\t", $record);
		if ( $ID =~ /(.+?)\|([^\|]+?)(_?{.+?})?$/ ) {
			$ID = $1;
			$prot = $2;
			$lineageByKeyProt{$ID}{$prot} = $lineage;
		} else {
			die("$0 ERROR: '$ID' is invalid.\n");
		}
	}
	close(RES);
}

$applicableProts == 0; $numGTs = 0;
shift(@fieldNames); $totalProts = scalar(@fieldNames);
@validProts{@fieldNames} = ();
@allProts{@fieldNames} = ();
foreach $ID ( keys(%lineageByKeyProt) ) {
	@prots = keys(%{$lineageByKeyProt{$ID}});
	%gtList = (); %prevList = (); $applicableProts = 0;
	for( $i=0; $i < scalar(@prots); $i++ ) {
		$prot = $prots[$i];
		$allProts{$prot}++;
		%gtList = ();
		if ( !exists($validProts{$prot}) ) {
			next;
		}
		$lineage = $lineageByKeyProt{$ID}{$prot};
		@gtArray = keys(%{$genotypes{uc($prot)}{uc($lineage)}}); 
		foreach $genotype ( @gtArray ) {
			if ( $i == 0 || exists($prevList{$genotype}) ) {
				$gtList{$genotype}=1;

			}
		}
#		print $prot,"\t",scalar(@gtArray),"\t",scalar(keys(%gtList)),"\n";
		%prevList = %gtList; $applicableProts++;
	}
	@possibleGTs = sort(keys(%gtList));
#	print $applicableProts,"\t",$totalProts,"\n";

	# We have the required number for the genotyping system
	$numGTs = scalar(@possibleGTs);
	if ( $applicableProts == $totalProts ) {
		if ( $numGTs == 1 ) {
			$genotypeByID{$ID} = $possibleGTs[0];
		} elsif ( $numGTs == 0 ) {
			$genotypeByID{$ID} = 'NEW_GENOTYPE';
		} else {
			$genotypeByID{$ID} = join(";",@possibleGTs);
		}
	} else {
		if ( $applicableProts == 0 ) {
			$genotypeByID{$ID} = 'NO_DATA';
		} elsif ( $numGTs != 0 ) {
			$genotypeByID{$ID} = 'MISSING_DATA:'.join(";",@possibleGTs);
		} else {
			$genotypeByID{$ID} = 'UNKNOWN_GENOTYPE';
		}
	}
}

@allProtsS = sort(keys(%allProts));
print "ID\tGenotype";
foreach $p ( @allProtsS ) { print "\t$p"; }
print "\n";

foreach $ID ( sort { $genotypeByID{$a} cmp $genotypeByID{$b} } keys(%genotypeByID) ) {
	print $ID,"\t",$genotypeByID{$ID};
	foreach $p ( @allProtsS ) {
		if ( $lineageByKeyProt{$ID}{$p} ne "" ) {
			print "\t",$lineageByKeyProt{$ID}{$p};
		} else {
			print "\tn/a";
		} 
	}
	print "\n";
			
}
