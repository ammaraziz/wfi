#!/usr/bin/env perl
# Samuel Shepard - 9.2014
# Take interleaved sampling throughout a FASTA file.
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

use File::Basename;
use Getopt::Long;
GetOptions(	'groups|G=i'=> \$numberGroups,
		'fraction|F=i' => \$fraction,
		'fastQ|Q' => \$fastQ,
		'by-read-pairs|P' => \$byReadPairs,
		'read-zipped|Z' => \$readZipped,
		'underscore-header|U' => \$underscoreHeader,
       		'extension|X:s' => \$extension	
);

if (  scalar(@ARGV) < 2 ) {
	$message = "\n$0 <input.fasta> <out_prefix> [-G <#groups>|-F <denom-fraction>] [OPTIONS]\n";
	$message .= "\t-F|--fraction POSITIVE_NUMBER\t\tFraction of dataset, using denominator D: 1/D.\n";
	$message .= "\t-G|--groups POSITIVE_NUMBER\t\tNumber of datasets required.\n";
	$message .= "\t-X|--extension\t\t\t\tExtension for output samplings.\n";
	$message .= "\t-Q|--fastQ\t\t\t\tFastQ format for input and output.\n";
	$message .= "\t-P|--by-read-pairs\t\t\tFastQ format for IN/OUT, interleave by read molecular ID (implies -Q).\n";
	die($message."\n");
}
$PROGRAM = basename($0,'.pl');

if ( $byReadPairs ) {
	$fastQ = 1;
}

if ( $fastQ ) {
	$extension = 'fastq';
}

if ( defined($fraction) && defined($numberGroups) ) {
	die("$PROGRAM ERROR: specify Fraction OR the number of Groups.\n");
} elsif ( defined($numberGroups) ) {
	if ( $numberGroups < 1 ) {
		die("ERROR: The number of groups must be more than zero.\n");
	} elsif( $numberGroups > 999 ) {
		print STDERR "$PROGRAM WARNING: groups currently capped to 200.\n";
		$numberGroups = 999;
	}
	$fraction = 0;
} else {
	$numberGroups = $fraction;
	if ( $numberGroups < 2 ) {
		die("$PROGRAM ERROR: The denominator must be more than one.\n");
	}
	$fraction = 1;
}

if ( $fraction ) {
	if ( $numberGroups =~ /1$/ ) {
		$suffix = 'st';
	} elsif ( $numberGroups =~ /2$/ ) {
		$suffix = 'nd';
	} elsif ( $numberGroups =~ /3$/ ) {
		$suffix = 'rd';
	} else {
		$suffix = 'th';
	}

	$filename = sprintf("%s_%d%s",$ARGV[1],$numberGroups,$suffix);
	if ( defined($extension) ) {
		$filename .= '.'.$extension; 
	} else {
		$filename .= '.fasta';
	}
	open($handle, '>', $filename ) or die("$PROGRAM ERROR: cannot open $filename\n");
} else {
	@handles = @count = ();
	for($i = 0;$i < $numberGroups; $i++ ) {
		$filename = sprintf("%s_%04d",$ARGV[1],($i+1));
		if ( defined($extension) ) {
			$filename .= '.'.$extension; 
		} else {
			$filename .= '.fasta';
		}
		open( $handles[$i], '>', $filename ) or die("$PROGRAM ERROR: Cannot open $filename\n");
		$files[$i] = $filename;
		$count[$i] = 0;
	}
}

# process parameters
chomp(@ARGV);
if ( $readZipped ) {
	open( IN, "zcat $ARGV[0] |" ) or die("$PROGRAM ERROR: Could not open $ARGV[0].\n");
} else {
	open( IN, '<', $ARGV[0] ) or die("$PROGRAM ERROR: Could not open $ARGV[0].\n");
}
$id = 0;
if ( $fastQ ) {
	$/ = "\n"; 
	if ( $byReadPairs ) {
		%indexByMolID = ();
		$REgetMolID = qr/@(.+?)[_ ][123]:.+/;
		while($hdr=<IN>) {
			$seq=<IN>;
			$junk=<IN>;
			$quality=<IN>; chomp($quality);
			if ( $hdr =~ $REgetMolID ) {
				$molID = $1;
				if ( defined($indexByMolID{$molID}) ) {
					$index = $indexByMolID{$molID};
				} else {
					$index = $id % $numberGroups;
					$indexByMolID{$molID} = $index;
					$id++;
					$count[$index]++;
				}
			} else {
				die("Irregular header for fastQ read pairs.\n");
			}
						
			if ( !$fraction ) {
				$handle = $handles[$index];
				print $handle $hdr,$seq,$junk,$quality,"\n";
			} elsif( $index == 0 ) {
				print $handle $hdr,$seq,$junk,$quality,"\n";
			}
		}
	} else {
		while($hdr=<IN>) {
			$seq=<IN>;
			$junk=<IN>;
			$quality=<IN>; chomp($quality);

			$index = $id % $numberGroups;
			$id++;
			$count[$index]++;
			if ( !$fraction ) {
				$handle = $handles[$index];
				print $handle $hdr,$seq,$junk,$quality,"\n";
			} elsif( $index == 0 ) {
				print $handle $hdr,$seq,$junk,$quality,"\n";
			}
		}
	}
} else {
	$/ = ">";
	while( $record = <IN> ) {
		chomp($record);
		@lines = split(/\r\n|\n|\r/, $record);
		$header = shift(@lines);
		if ( defined($underscoreHeader) ) { $header =~ tr/ /_/; }
		$sequence = lc(join('',@lines));

		$length = length($sequence);
		if ( $length == 0 ) {
			next;	
		}

		$index = $id % $numberGroups;
		$id++;
		$count[$index]++;
		if ( !$fraction ) {
			$handle = $handles[$index];
			print $handle '>',$header,"\n",$sequence,"\n";
		} elsif( $index == 0 ) {
			print $handle '>',$header,"\n",$sequence,"\n";
		}
	}
}
close(IN);
if ( $fraction ) {
	close($handle);
	print "\n Total\t  Got\tSample Name\n";
	print '----------------------------------------------------',"\n";
	printf("%6d\t%5d\t%s\n",$id,$count[0],$filename);
	print '----------------------------------------------------',"\n";
} else {
	foreach $handle (@handles) {
		close($handle);
	}
	print "\n Total\t  Got\tSample Name\n";
	print '----------------------------------------------------',"\n";
	for($i = 0;$i < $numberGroups;$i++) {
		printf("%6d\t%5d\t%s\n",$id,$count[$i],$files[$i]);
	}
	print '----------------------------------------------------',"\n";
}


