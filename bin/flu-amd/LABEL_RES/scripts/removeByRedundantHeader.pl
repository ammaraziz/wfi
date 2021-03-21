#!/usr/bin/env perl
# removeByRedundantHeader -  Version 1.0
# Remove sequences from a FASTA given the header is entirely redundant.
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
GetOptions(
		'underline-spaces|U'=>\$underline
		);

if ( -t STDIN and not @ARGV ) {
	die("Usage:\n\tperl $0 <file.fa> [-U]\n");
}

# PROCESS fasta data
$/ = ">"; %headers = ();
while( $record = <> ) {
	chomp($record);
	@lines = split(/\n/, $record);
	$header = shift(@lines);
	$sequence = lc(join('',@lines));
	if ( $underline ) { $header =~ tr/ /_/d; }

	if ( length($sequence) == 0 ) {
		next;
	} elsif ( !exists($headers{$header}) )  {
		$headers{$header} = $sequence;
	}
}

# OUTPUT non-redundant from the hash
foreach $header ( keys(%headers) ) {
	print '>',$header,"\n";
	print $headers{$header},"\n";
}
####################
