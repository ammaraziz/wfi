#!/usr/bin/env perl
# annotateByFragment -  Version 1.0
# Annotates sequences containiner N-mer without storing fragment itself.
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


my $message;
if ( scalar(@ARGV) != 2 ) {
	$message = "\nUsage:\n\tperl $0 <input_fasta> <hash_file>\n";
	die($message."\n");
}

$/ = "\n";
my %fragments = ();
my $line = ();
my ($key,$annot,$length) = ('','',0);
open(IN,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
while($line=<IN>) {
	chomp($line);
	($key,$annot,$length) = split("\t",$line);
	if ( $annot ne '' && $key ne '' && $length > 0) {
		$fragments{$length}{$key} = $annot;
	}
}
close(IN);

if ( scalar(keys(%fragments)) <= 0 ) {
	print STDERR "No fragments read.\n";
	exit(0);
} 

$/ = '>';
open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
my @records = <IN>; 
chomp(@records);
close(IN);

my @lines = ();
my $header = '';
my $sequence = '';
my $seq2 = '';
my $hash = '';
my $found = 0;
my @fragment_lengths = sort { $b <=> $a } keys(%fragments);

#open(OUT,'>',$ARGV[0]) or die("Cannot open $ARGV[0] for writing.\n");
my $fragment = '';
foreach my $record ( @records ) {
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = join('',@lines);
	$length = length($sequence);
	if ( $length == 0 ) { next; }

	$seq2 = uc($sequence);
	$seq2 =~ tr/ :.~-//d;
	$found = 0;
	if ( $seq2 ne '' ) {
		foreach my $L ( @fragment_lengths ) {
			if ( $L <= $length ) {
				foreach my $p ( 0..($length-$L) ) {
					$fragment = substr($seq2,$p,$L);
					$hash = sha1_hex($fragment);
					if ( defined($fragments{$L}{$hash}) ) {
						print STDOUT $header,"\t",$fragments{$L}{$hash},"\t",$hash,"\n";
						last;
					}		
				}
				if ( $found ) { last; }
			}
		}
	}

	if ( !$found ) {
		print '';
		#print OUT '>',$header,"\n",$sequence,"\n";
	}
}
#close(OUT);
