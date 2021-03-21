#!/usr/bin/env perl
# fastaExtractor -  Version 1.0
# Extract FASTA sequences given a line-delimited list of sequences.
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
GetOptions( 'delimiter|D:s'=> \$delim, 'field|F:i' => \$field, 'deep-search|S' => \$deepSearch, 'trim-field|T' => \$trimField, 'reverse-extraction|R' => \$reverse, 'ignore-case|C' =>\$noCase, 'ignore-annotation|A' => \$ignoreAnnotation, 'not-found|N' => \$notFound, 'not-used|U' => \$notUsed );

if ( scalar(@ARGV) == 1 && not -t STDIN ) {
	$handle = 'STDIN';
} elsif ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <file.fasta> <line-delimited.txt>\n";
	$message .= "\t\t-T|--trim-field\t\t\t\tTrim field of white space.\n";
	$message .= "\t\t-D|--delimiter CHARACTER\t\tSingle character delimiter within FASTA header. Default is TAB.\n";
	$message .= "\t\t-F|--field POSITIVE_NUMBER\t\tHeader field to use for pairing records. Default is line.\n";
	$message .= "\t\t-A|--ignore-annotation\t\t\tIgnores annotations for pairings.\n";
	$message .= "\t\t-R|--reverse-extraction\t\t\tTurns the ID file into an exclusion list rather than an inclusion one.\n";
	$message .= "\t\t-C|--ignore-case\t\t\tIgnores header case.\n";
	$message .= "\t\t-N|--not-found\t\t\t\tShow data not found in list.\n";
	$message .= "\t\t-U|--not-used\t\t\t\tShow data not used in FASTA.\n";
	die($message."\n");
} else {
	open(IN, '<',$ARGV[1]) or die("$0 ERROR: Cannot open $ARGV[1].\n");
	$handle = 'IN';
}

if ( !defined($field) && !defined($delim) ) {
	$delimited = 0;
} else {
	if ( !defined($field) ) {
		$field = 1;
	} elsif($field == 0 ) {
		die("$0 ERROR: field must be specified.\n");
	} elsif($field < 0) {
		die("$0 ERROR: field must be a positive number.\n");
	} else {
		$field -= 1;
	}

	if ( !defined($delim) ) {
		$delim = "\t";
	} elsif( $delim eq '' ) {
		die("$0 ERROR: No delimiter argument detected.\n");
	} elsif( length($delim) > 1 ) {
		die("$0 ERROR: single character delimiter expected instead of '$delim'.\n");
	}
	$delimited = 1;
}

# READ strain information
$/ = "\n"; %ids = (); %found = ();
while ( $line = <$handle> ) {
	chomp($line);

	if ( $noCase ) {
		$line = lc($line);
	}

	if ( defined($ignoreAnnotation) && $line =~ /_?\{.+?\}/ ) {
		$line =~ s/_?\{.+?\}//;
	}

	if ( $delimited ) {
		if ( index($line,$delim) == -1 && $delim ne '' ) {
			die("$0 ERROR: delimiter '$delim' not found in header.\n");
		}

		if ( $delim eq '|' ) {
			@fields = split( /\|/, $line );
		} elsif ( $delim ne "\t" ) {
			@fields = split( /\Q$delim\E/, $line );
		} else {
			@fields = split( /\t/, $line );
		}
		if ( $trimField ) {
			$ids{trim($fields[$field])} = trim($fields[$field]);
		} else {
			$ids{$fields[$field]} = $fields[$field];
		}
	} else {
		if ( $trimField ) {
			$ids{trim($line)} = trim($line);
		} else {
			$ids{$line} = $line;
		}
	}

}
close($handle);

$tmp = '';
if ( $notFound ) {
	foreach $id ( keys(%ids) ) {
		$found{$id} = 1;
		$tmp = $id;
	}
}

# READ fasta data
open(IN,'<', $ARGV[0]) or die("$0 ERROR: Cannot open $ARGV[0].\n");
$/ = ">"; @keys = keys(%ids);
while($record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$ID = $header = shift(@lines);
	$sequence = join('',@lines);
	
	if ( $noCase ) {
		$ID = lc($header);
	}

	if ( defined($ignoreAnnotation) && $ID =~ /_?\{.+?\}/ ) {
		$ID =~ s/_?\{.+?\}//;
	}

	$result = headerInDB(\%ids,\@keys,$ID);
	if ( length($sequence) == 0 ) {
		next;
	} elsif( ! $result != ! $reverse ) {
		if ( $notFound ) {
			delete $found{$result};
		}
		print '>',$header,"\n",$sequence,"\n";
	} else {
		if ( $notUsed ) {
			print STDERR "Scrubbed:\t",$header,"\n";
		}
		next;
	}

}
close(IN);

if ( $notFound ) {
	foreach $id ( keys(%found) ) {
		print STDERR $id;
	}
}


# FNC - search for a header in the header hash database
sub headerInDB(\%\@$) {
	my $ids = $_[0];
	my $keys = $_[1];
	my $header = $_[2];
	my $id = '';

	if ( exists($ids->{$header}) ) {
		return $header;
	}

	if ( $deepSearch ) {
		foreach $id ( @{$keys} ) {
			if ( $header =~ /\Q$id\E/ ) {
				return $id;
			}
		}
	}
	return 0;
}

# FNC - trim function.
# Removes whitespace from the start and end of the string
 sub trim($) {
 	my $string = shift;
	$string =~ /^\s*(.*?)\s*$/;
 	return $1;
}
####################
