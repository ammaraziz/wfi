#!/usr/bin/env perl
# sortByAnnotation - Version 1.0
# Sort fasta sequences by their taxa annotations.
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
if ( scalar(@ARGV) != 1 ) {
	die("Usage:\n\t$0 <file.fas>\n");
}


open( IN, '<', $ARGV[0] ) or die("$0 ERROR: Cannot open $ARGV[0].\n");

# READ fasta sequences and their annotations
$/ = ">"; %taxa = (); 
while($record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = join('',@lines);
	eval($strip);

	if ( $header =~ /{(.*)}$/ ) {
		$taxon = $1;
	}
	
	if ( length($sequence) == 0 ) {
		next;
	}
	
	$taxa{$taxon} = 1;
	$sequences{$taxon}{$header} = $sequence;
}
close(IN);

# PRINT sorted sequences
@sorted = sort( keys(%taxa) );
foreach $s ( @sorted ) {
	foreach $t ( keys(%{$sequences{$s}}) ) {
		print '>',$t,"\n";
		print $sequences{$s}{$t},"\n";
	}
}
####################
