#!/usr/bin/env perl
# reviseTaxa -  Version 1.1
# Allows taxa annotation manipulations for fasta files.
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
use File::Basename;
GetOptions(
		'find|F:s'=> \$find,
		'in-place|I'=>\$inPlace,
		'replace|R:s' => \$replace,
		'join-to-end|J=s' => \$joinAnnot,
		'confirm-prediction|C' => \$confirm,	
		'delete-prev|D' => \$deletePrev,
		'delete-single|S'=> \$deleteSingle,
		'add-annot|A=s' => \$addAnnot,
		'add-replace|L=s' => \$addReplace,
		'fuzzy-match|Z' => \$fuzzyMatch,
		'order-mode|O:s' => \$orderMode,
		'match-files|M' => \$matchFiles,
		'previous-infix|N' => \$prevInfix,
		'prefix|P=s' => \$prefix,
		'suffix|X=s' => \$suffix
	  );

if ( (scalar(@ARGV) != 1 && !$matchFiles) || ($matchFiles && scalar(@ARGV) < 2) ) {
	$message = "Usage:\n\tperl $0 <input.fasta> [OPTIONS] [-M <file1 file2 ...>]\n";
	$message .= "\t\t--confirm-prediction|-C\t\tConfirm predicted annotations.\n";
	$message .= "\t\t--delete-prev|-D\t\tDelete previous annotations where there are two.\n";
	$message .= "\t\t--delete-single|-S\t\tDelete previous annotations where there is one, filters AFTER delete-prev.\n";
	$message .= "\t\t--find|-F <TEXT>\t\tSelect sequences including this annotation.\n";
	$message .= "\t\t--replace|-R <TEXT>\t\tReplace the annotation with TEXT.\n";
	$message .= "\t\t--add-annot|-A <FILE>\t\tAdd annotations based on tab delimited file (ID\tANNOT).\n";
	$message .= "\t\t--fuzzy-match|-Z\t\tSearches for IDs in FASTA (-A option), matching if header contains the ID.\n";
	$message .= "\t\t--order-mode|-O <out.file>\tAnnotation with ordinals for truncated names.\n";
	$message .= "\t\t--match-files|-M <file1 ...>\tOutput filenames containing annotation names in the input file.\n";
	$message .= "\t\t--previous-infix|N\t\tSuffix and prefix only apply to a 'previous' or first of a doublet annotation.\n";
	$message .= "\t\t--prefix|-P <TEXT>\t\tPrefix for matching filenames with annotations OR adds prefix to header.\n";
	$message .= "\t\t--suffix|-X <TEXT>\t\tSuffix for matching filenames with annotations OR adds suffix to header.\n";
	$message .= "\t\t--in-place|-I\t\t\tRevise files in-place (could be dangerous).\n";
	$message .= "\t\t--join-to-end|-J <TEXT>\t\tJoin annotation to end of the header (similar to add but without a file).\n";

	die($message);
}

if ( defined($orderMode) ) {
	open(ORD, '>', $orderMode ) or die("$0 ERROR: Cannot open $orderMode.\n");
}

if ( $addReplace ) { $addAnnot = $addReplace; }
if ( $addAnnot ) {
	open(ANNOT,'<', $addAnnot ) or die("$0 ERROR: Cannot open $addAnnot.\n");
	%annotMap = (); $/ = "\n";
	while( $line = <ANNOT> ) {
		chomp($line);
		($key, $value) = split(/\t/, $line);
		$annotMap{uc($key)} = $value;
	}
	close(ANNOT);
	@annotIDs = keys(%annotMap);
}

$/ = ">";
open( IN, '<', $ARGV[0] ) or die("$0 ERROR: Cannot open $ARGV[0] for reading.\n");
@records = <IN>;
close(IN);
%match = ();

if ( $inPlace ) {
	open( OUT, '>', $ARGV[0] ) or die("$0 ERROR: Cannot open $ARGV[0] for writing.\n");
}

# PROCESS stored fasta information
foreach $record ( @records ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$id = shift(@lines);
	$sequence = join('',@lines);

	if ( length($sequence) == 0 ) {
		next;
	}

	# Add annotations from a file
	if ( $addAnnot ) {
		if ( $newID = headerInDB(\%annotMap,\@annotIDs,uc($id)) ) {
			$newID = $annotMap{uc($newID)};
			if ( $addReplace && $id =~ m/{([^{}]*)}\s*$/ ) {
				$id =~ s/{([^{}]*)}\s*$/{$newID}/;

			} else {
				$id = $id."{$newID}";
			}
		}
	}
	
	# join to the end
	if ( defined($joinAnnot) ) {
		$id .= '{'.$joinAnnot.'}';
	}

	# Get rid of a previous annotation.
	if ( $deletePrev ) {
		$id =~ s/{.+?}.*?{(.+?)}/{$1}/;
	}

	if ( $deleteSingle ) {
		$id =~ s/_?{.+?}//;
	}


	# Confirm prediction if they exist.
	if ( $confirm ) {
		$id =~ s/{PRED:(.+?)}$/{$1}/;
	}

	# Find and replace.
	if ( defined($find) && defined($replace) ) {
		$id =~ s/\Q{$find}\E/{$replace}/;

	# Always replace.
	} elsif ( defined($replace ) ) {
		if ( $id !~ s/{([^{}]*)}\s*$/{$replace}/ ) {
			$id .= "{$replace}";
		}
	}
	
	if ( $orderMode ) {
		$id =~ /{(.+?)}/;
		$annot = $1;
		if ( !defined( $count{$annot} ) ) {
			$count{$annot} = 1;
		} else {
			$count{$annot}++;
		}


		print ORD $id,"\t",$count{$annot},'_',$annot,"\n";
		$id = $count{$annot}.'_'.$annot;
	}

	# Check for prefix or suffix
	if ( (defined($prefix) || defined($suffix)) && !$matchFiles ) {
		if ( defined($prevInfix) ) {
			$id =~ s/{(.+?)}\{/{$prefix$1$suffix}{/;
		} else {
			$id =~ s/{(.+?)}/{$prefix$1$suffix}/;
		}
	}


	if ( $matchFiles ) {
		if ( $id =~ /{(.+?)}/ ) {
			$match{$1} = 1;
		}
	} elsif ( $inPlace ) {
		print OUT '>',$id,"\n",$sequence,"\n";
	} else {
		print '>',$id,"\n",$sequence,"\n";
	}
}
close(OUT);
if ( defined($orderMode) ) {
	close(ORD);
}

if ( $matchFiles ) {
	for($i = 1; $i < scalar(@ARGV); $i++ ) {
		$filename = basename($ARGV[$i]);
		
		foreach $annot ( keys(%match) ) {
			if ( $filename =~ /^$prefix$annot$suffix/i ) {
				print $ARGV[$i]," ";
				last;
			}
		}
	}
}

#
# FNC - search for a header in the header hash database
sub headerInDB(\%\@$) {
	my $ids = $_[0];
	my $keys = $_[1];
	my $header = $_[2];
	my $id = '';

	if ( exists($ids->{$header}) ) {
		return $header;
	}

	if ( $fuzzyMatch ) {
		foreach $id ( @{$keys} ) {
			if ( $header =~ /\Q$id\E/ ) {
				return $id;
			}
		}
	}
	return 0;
}
