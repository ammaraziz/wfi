#!/usr/bin/env perl
# removeInsertionStates - Version 1.0
# Removes insertion states from align2model file for quick alignment.
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
if ( scalar(@ARGV) < 1 ) {
	die("Usage:\n\t$0 <file1.fa> <...>\n");
}

# PROCESS fasta data
$/ = ">";
while($record = <> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = join('',@lines);

	if ( length($sequence) == 0 ) {
		next;
	} else {
		$sequence =~ tr/[a-z.]//d;
		print '>',$header,"\n";
		print $sequence,"\n";
	}
}
close(IN);
####################
