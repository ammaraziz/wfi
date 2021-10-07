#!/usr/bin/env perl
# buildDataMatrix -  Version 1.1
# Formats data into a matrix for input to a graphics program.

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
use File::Basename;
use Getopt::Long;

#TO-DO: Make normalized mode for data across all columns.

GetOptions(	'lengths|L' => \$addLengths,'shorten-taxa|S' => \$SHORTEN, 'score-field|F=i' => \$scoreField,
       		'null-model|N=s' => \$nullModel, 'column-minimum|M' => \$colMin, 'custom-reverse|C=s' => \$custRevPath
);
if ( scalar(@ARGV) < 2 ) {
	$message = "Usage:\n\tperl $0 <output_file> <scores.tab> ... <scoresN.tab>\n";
	$message .= "\t\t-S|--shorten-taxa\t\tShorten taxa name if over 10 characters and not a c-* or outlier pattern.\n";
	$message .= "\t\t-F|--score-field <#>\t\tScore field to build matrix (1-N, default 4).\n";
	$message .= "\t\t-L|--lengths\t\t\tAdd the length field to the table.\n";
	$message .= "\t\t-N|--null-model <null.mod>\tSpecify a null file.\n";
	$message .= "\t\t-M|--column-minimum\t\tCalculate a column for the vector minimum.\n";
	$message .= "\t\t-C|--custom-reverse <path>\tUse custom reverse corrected viterbi model. Path contains equivalent pHMMs.\n";
	die($message);
}

if ( !defined($SHORTEN) ) {
	$SHORTEN = 0;
} else {
	$SHORTEN = 1;
}
%scores = (); %lengths = ();
$outfile = shift(@ARGV);
%labelsByTaxa = %taxaByLabels = ();
$taxon = ''; $label = 0;

if ( !defined($scoreField) ) {
	$scoreField = 3;
} else {
	if ( $scoreField > 0 ) {
		$scoreField--;
	} else {
		die("`basename $0` ERROR: Score field must be a positive integer.\n");
	}
}

if ( defined($nullModel) ){
	$subtractNull = 1;
	open(NULL,'<',$nullModel) or die("`basename $0` ERROR: cannot open $nullModel for reading.\n");
	$header = <NULL>;
	while ( $line = <NULL> ) {
		chomp($line);
		@values = split("\t", $line);
		$nulls{$values[0]} = $values[$scoreField];
	}
	close(NULL);
} else {
	$subtractNull = 0;
}

# GATHER tabular data
if ( defined($colMin) ) {
	%minByID = ();
}
$/ = "\n";

if ( defined($custRevPath) ) {
	$xrev = 1;
} else {
	$xrev = 0;
}

foreach $file ( @ARGV ) {
	# PROCESS labels & first record
       	open( IN, '<', $file ) or die("$0 ERROR: Cannot open $file.\n");

	$header = <IN>;
	$line = <IN>;
	chomp($line);
	@values = split("\t", $line);
	$file = basename($file);
	if ( $file =~ /(.+?)_hmm.*?\.tab/ ) {
		$taxon = $1;
	} else {
		die("$0 ERROR: Malformed file name ($file).\n");
	}

	if ( $xrev ) {
		@revValues = ();
		%revScores = ();
		$revFile = $custRevPath . '/' . $file;
		open(REV,'<',$revFile) or die("Cannot open $revFile\n");
		$revHeader = <REV>;
		while($revLine =<REV> ) {
			chomp($revLine);
			@revValues = split("\t",$revLine);
			$revScores{$revValues[0]} = $revValues[$scoreField];
		}	
		close(REV);
	}	

	$id = $values[0];
	$lengths{$id} = $values[1];
	if ( $subtractNull ) {
		$scores{$id}{$label} = sprintf("%.2f",$values[$scoreField] - $nulls{$id});
	} elsif ( $xrev ) {
		$scores{$id}{$label} = sprintf("%.2f",$values[$scoreField] - $revScores{$id});
	} else {
		$scores{$id}{$label} = $values[$scoreField];
	}
	$labelsByTaxa{$taxon} = $label;
	$taxaByLabels{$label} = $taxon;

	if ( $colMin ) {
		if ( !defined($minByID{$id}) || $scores{$id}{$label} < $minByID{$id} ) {
			$minByID{$id} = $scores{$id}{$label};
		}
	}

	# FINISH remaining records
	while ( $line = <IN> ) {
		chomp($line);
		@values = split("\t", $line);

		$id = $values[0];
		$lengths{$id} = $values[1];
		if ( $subtractNull ) {
			$scores{$id}{$label} = sprintf("%.2f",$values[$scoreField] - $nulls{$id});
		} elsif ( $xrev ) {
			$scores{$id}{$label} = sprintf("%.2f",$values[$scoreField] - $revScores{$id});
		} else {
			$scores{$id}{$label} = $values[$scoreField];
		}

		if ( $colMin ) {
			if ( !defined($minByID{$id}) || $scores{$id}{$label} < $minByID{$id} ) {
				$minByID{$id} = $scores{$id}{$label};
			}
		}
	}
	close(IN);
	$label++;
}
$numberProfiles = $label;
@taxa = sort(keys(%labelsByTaxa));
@labels = sort { $a <=> $b } keys(%taxaByLabels);

# SANITY check
$numTrainingProfiles = scalar(@labels);
if ( $numTrainingProfiles != $numberProfiles ) {
	die("$0 ERROR: Number of training profiles ($numTrainingProfiles) different from test data profiles ($numberProfiles).\n");
}

# OUTPUT data matrix
open( TSV, '>', $outfile) or die("$0 ERROR: Cannot open $outfile for writing.\n");
print TSV 'ID',"\t",'ANNOTATION';
for( $h = 0; $h < $numberProfiles; $h++ ) {
	print TSV "\t",'HMM_',$taxaByLabels{$labelsByTaxa{$taxa[$h]}};
}
if ( $addLengths ) {
	print TSV "\tLength";
}

if( defined($colMin) ) {
	print TSV "\tMinCol";
}
print TSV "\n";

@sorted = sort(keys(%scores));
foreach $id ( @sorted ) {
	if ( $id =~ /{([^{}]*)}\s*$/ ) {
		$annot = $1;
	} else {
		$annot = 'Unknown';
	}
	print TSV $id,"\t",$annot;
	for( $h = 0; $h < $numberProfiles; $h++ ) {
		print TSV "\t",$scores{$id}{$labelsByTaxa{$taxa[$h]}};
	}
	if ( $addLengths ) {
		print TSV "\t",$lengths{$id};
	}

	if ( $colMin ) {
		print TSV "\t",$minByID{$id};	
	}

	print TSV "\n";
}
close(TSV);
####################
