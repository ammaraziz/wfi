#!/usr/bin/env perl
# stripSequences - Version 1.0
# Strips unwanted chracters from the fasta sequence (such as alignment characters).
# May also "fix" the header via trimming, underlining spaces & removing commas/apostrophes.
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
use Getopt::Long;
GetOptions( 'fix-header|F' => \$fixHeader, 'strip-lower|L' => \$stripLower );
if ( (defined($stripLower) && scalar(@ARGV) != 1) || (!defined($stripLower) && scalar(@ARGV) != 2) ) {
	die("Usage:\n\t$0 [-F] <file.fas> {-L|<quoted_characters_to_delete>}\n");
}


open( IN, '<', $ARGV[0] ) or die("$0 ERROR: Cannot open $ARGV[0].\n");

# PREPARE the strip deletion
$strip = '$sequence =~ tr/'.$ARGV[1].'//d;';

# PROCESS fasta data
$/ = ">";
while($record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	if ( $fixHeader ) {
		$header =~ s/^\s*(.*?)\s*$/\1/;
		$header =~ s/[\s:]/_/g;
		$header =~ tr/',//d;
	}
	$sequence = join('',@lines);
	if ( $stripLower ) {
		$sequence =~ tr/[a-z]//d;
	} else {
		eval($strip);
	}

	if ( length($sequence) == 0 ) {
		next;
	} else {
		print '>',$header,"\n";
		print $sequence,"\n";
	}
}
close(IN);
####################
