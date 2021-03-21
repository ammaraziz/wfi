#!/usr/bin/env perl
# parseScores -  Version 1.0
# Converts SAM score distance files to a tab-delimited format.
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
	die("Usage:\n\tperl $0 <scores.dist ...>\n");
}
$/ = "\n"; @headers = @fields = ();
while ( $line = <> ) {
	chomp($line);
	if ( $line =~ m/^% Sequence ID/ ) {
		$line =~ tr/%//d;
		$line = trim($line);
		@headers  = split(/\s{2,}/, $line );
		print $headers[0],"\t",$headers[1],"\t",$headers[2],"\t",$headers[3],"\n";
		last;
	} elsif ( $line !~ m/^%/ ) {
		@fields = split(/\s+/, $line );
		print $fields[0],"\t",$fields[1],"\t",$fields[2],"\t",$fields[3],"\n";
	}
}

while ( $line = <> ) {
	chomp($line);
	if ( $line !~ m/^%/ ) {
		@fields = split(/\s+/, $line );
		print $fields[0],"\t",$fields[1],"\t",$fields[2],"\t",$fields[3],"\n";
	}
}

# FNC - Trim function.
# Removes whitespace from the start and end of the string
 sub trim($) {
 	my $string = shift;
	$string =~ /^\s*(.*?)\s*$/;
 	return $1;
}
####################
