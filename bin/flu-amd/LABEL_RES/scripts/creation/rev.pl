#!/usr/bin/env perl
# Sam Shepard
# rev.pl 
# Reverse complements the data.
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

use Getopt::Long;
GetOptions(	'fastQ|Q'=> \$fastQ, 'reverse-only|R' => \$reverseOnly, 'single-line-mode|L' => \$singleLine );

if ( scalar(@ARGV) != 1 ) {
	die("Usage:\n\tperl $0 <input.fasta> [-Q] [-R] [-L]\n");
}

if ( defined($reverseOnly) ) {
	$takeComplement = 0;
} else {
	$takeComplement = 1;
}


if ( $singleLine ) {
	$sequence = $ARGV[0];
	$sequence = reverse( $sequence );
	if ( $takeComplement ) {
		$sequence =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
	}
	print $sequence,"\n";
	exit;
}

if ( $fastQ ) {
	$/ = "\n";
} else {
	$/ = ">";
}

while( $record = <> ) {
	chomp($record);
	if ( $fastQ ) {
		$header = $record;
		$sequence=<IN>; chomp($sequence);
		$junk=<IN>; chomp($junk);
		$quality=<IN>; chomp($quality);
		$quality=reverse($quality);
	} else {
		@lines = split(/\r\n|\n|\r/, $record);
		$header = shift(@lines);
		$sequence = lc(join('',@lines));
	}

	$length = length($sequence);
	if ( $length == 0 ) {
		next;	
	}
	$sequence = reverse( $sequence );

	if ( $takeComplement ) {
		$sequence =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
	}

	if ( $fastQ ) {
		print $header,"\n",$sequence,"\n",$junk,"\n",$quality,"\n";
	} else {
		print '>',$header,"\n",$sequence,"\n";
	}
}
