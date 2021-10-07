#!/usr/bin/env perl
# Generate Fragment Hash  -  Version 1.0
# Helps generate hashes for outlier exceptions.
#
# Copyright (C) 2019, Centers for Disease Control & Prevention
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
use Digest::SHA qw(sha1_hex);
use strict; use warnings;

use Getopt::Long;
my $cladePurity;
GetOptions(	'purity-of-clade|C=s' =>\$cladePurity );

my $message;
if ( scalar(@ARGV) != 3 ) {
	$message = "\nUsage:\n\tperl $0 <input_fasta> <start> <stop>\n";
	die($message."\n");
}

my @lines = ();
my $header = '';
my $sequence = '';
my $seq2 = '';
my $hash = '';
my $annot = '';
my $length = 0;

my ($start,$stop) = ($ARGV[1],$ARGV[2]);
my ($pos,$L);
if ( defined($start) && defined($stop) ) {
	if ( $stop < $start ) {
		$pos = int($stop) - 1;
		$L = $start - $stop + 1;
	} else {
		$pos = int($start) - 1;
		$L = $stop - $start + 1;
	}
	if ( $pos < 0 ) { $pos = 0; }
}

$/ = '>'; my $record = '';
my %hash = ();
my %clade_hash = ();
my %target_clade = ();
my %hash_seq = ();
open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while ( $record = <IN>  ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = join('',@lines);
	$length = length($sequence);
	$annot = 'UNCLASSIFIED';
	if ( $length == 0 ) { next; }

	$seq2 = uc(substr($sequence,$pos,$L));
	if ( length($seq2) != $L ) { 
		next;
	} else {	
		$seq2 =~ tr/ :.~-//d;
		if ( $seq2 ne '' ) {
			$length = length($seq2);
			$hash = sha1_hex($seq2);
			if ( $header =~ /\{([^}]+)\}\s*$/ ) {
				$annot = $1;
			}
			$hash{$header} = $hash."\t".$annot."\t".$length;
			if ( defined($cladePurity) ) {
				if ( $cladePurity eq $annot ) {
					$target_clade{$hash}++;
					$hash_seq{$hash} = $seq2;
				} else {
					$clade_hash{$hash}{$annot}++;
				}
			}
		}
	}
}
close(IN);

if ( defined($cladePurity) ) {
	foreach my $h ( sort { $target_clade{$b} <=> $target_clade{$a} } keys(%target_clade) ) {
		print STDOUT "target\t",$cladePurity,"\t",$h,"\t",$target_clade{$h},"\n";
		if ( defined($clade_hash{$h}) ) {
			foreach my $clade ( sort { $clade_hash{$h}{$b} <=> $clade_hash{$h}{$a} } keys(%{$clade_hash{$h}}) ) {
				print STDOUT "\t",$clade,"\t",$h,"\t",$clade_hash{$h}{$clade},"\n";
			
			}
		}
		print STDOUT "$hash_seq{$h}\n";
	}
} else {
	foreach my $id ( sort(keys(%hash)) ) {
		print STDOUT $id,"\t",$hash{$id},"\n";
	}
}
