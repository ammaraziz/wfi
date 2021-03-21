#!/usr/bin/env perl
# partitionTaxa - Version 1.0
# Partitions sequences by named taxa in the {taxa_name} format.
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
		'prefix|P=s' => \$prefix,
		'suffic|S=s' => \$suffix
		);

if ( scalar(@ARGV) != 2 ) {
	$message 	= "Usage:\n\tperl $0 <input.fasta> <directory_path_output> [OPTIONS]\n";
	$message 	.= "\t\t-P|--prefix\tPrefix for output files.\n";
	$message	.= "\t\t-S|--suffix\tSuffix for output files (before extension).\n";
	die($message);
}

open(IN, '<', $ARGV[0] ) or die("$0 ERROR: Cannot open $ARGV[0].\n");

if ( substr($ARGV[1], -1) eq '/' ) {
	$path = $ARGV[1];
} else {
	$path = $ARGV[1].'/';
}
@parts = split(/\./, $ARGV[0]);
if ( scalar(@parts) == 1 ) {
	$ext = '';
} else {
	$ext = '.'.$parts[$#parts];
}

# PROCESS fasta data
$/ = ">"; %taxa = (); %taxaHandles = ();
while($record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = join('',@lines);
	$length = length($sequence);

	if ( $length == 0 ) {
		next;
	}

	if ( $header =~ /{([^{}]*)}$/ ) {
		$taxon = $1;
		$taxon =~ tr/\//-/;
#		$taxon =~ tr/+//d;
#		if ( length($taxon) > 10 && $taxon !~ /-like/ && $taxon !~ /outlier/i ) {
#			@parts = split('_',$taxon);
#			if ( scalar(@parts) > 1 ) {
#				$taxon = $parts[0].'-'.$parts[$#parts];
#			} else {
#				$taxon = substr($taxon,0,1).'-'.substr($taxon,-1);
#			}
#		}
	} else {
		$taxon = 'unknown';
	}

	if ( !exists($taxa{$taxon}) ) {
		if ( defined($prefix) ) {
			$fullname=$path.$prefix.$taxon;
		} else {
			$fullname=$path.$taxon;
		}

		if ( defined($suffix) ) {
			$fullname .= $suffix.$ext;
		} else {
			$fullname .= $ext;
		}

		open($taxaHandles{$taxon},'>', $fullname) or die("$0 ERROR: Cannot open $fullname.\n");
		$handle = $taxaHandles{$taxon};
		print $handle '>',$header,"\n",$sequence,"\n";
		$taxa{$taxon} = 1;
	} else {
		$taxa{$taxon}++;
		$handle = $taxaHandles{$taxon};
		print $handle '>',$header,"\n",$sequence,"\n";
	}
}
close(IN);

foreach $handle ( keys( %taxaHandles ) ) {
	close($handle);
}
####################
