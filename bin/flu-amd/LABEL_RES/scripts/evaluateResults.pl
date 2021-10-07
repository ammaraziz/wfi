#!/usr/bin/env perl
# evaluateResults.pl - Version 1.0
# Generates a shogun cmdline_static compatible script using input params.
# Outputs a revised result file using ID information.
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
GetOptions( 	'terse|T' => \$terse, 'header|H' => \$header, 'grouped|G' => \$grouped, 
		'suppress|S' => \$suppress, 'just-seq-header|J' => \$justSeqHeader,
       		'path-list|P=s' => \$pathList	);

if ( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <results.txt> [OPTIONS]\n";
	$message .= "\t-G|--grouped\t\tResults were grouped using a hierarhical decimal scheme and prefix.\n";
	$message .= "\t-H|--header\t\tData has a header row.\n";
	$message .= "\t-S|--suppress\t\tSuppress printing of incorrect values.\n";
	$message .= "\t-T|--terse\t\tSummary output is just number correct.\n";
	$message .= "\t-J|--just-seq-header\tJust print out sequence headers for incorrect values.\n";
	$message .= "\t-P|--path-list\t\tUses a path list to evaluate correctness.\n";
	die($message."\n");
}

open(IN, '<', $ARGV[0]) or die("Cannot open $ARGV[0].\n");
$correct = $total = 0;
if ( defined($header) ) {
	$header = <IN>;
}

if ( defined($pathList) ) {
	open(PL,'<',$pathList) or die("Cannot open path list $pathList for reading.\n");
	$/ = "\n"; %pathArrayByAnnot = ();
	while($line = <PL>) {
		chomp($line);
		@pieces = split('/', $line);
		$annot = pop(@pieces);
		$junk = shift(@pieces);
		if ( scalar(@pieces) > 1 ) {
			$junk = shift(@pieces);
			$pathArrayByAnnot{$annot} = [@pieces];
		}
	}
	close(PL);
	$pathList = 1;
} else {
	$pathList = 0;
}

%annots = ();
%corByAnno = ();
while( $line = <IN> ) {
	chomp($line);

	$total++;

	@fields = split("\t", $line);
	$fields[0] = trim($fields[0]);
	$pred = trim($fields[1]);
	if ( $fields[0] =~ /{([^{}]*)}\s*$/ ) {
		$anno = $1;
#		if ( !defined($annots{$anno}) ) { $annots{$anno} = 1; } else { $annots{$anno}++; }
	} else {
		$anno = 'unknown';
#		if ( !defined($annots{$anno}) ) { $annots{$anno} = 1; } else { $annots{$anno}++; }
#		if ( $pred =~ /-like$/ || $pred =~ /outlier$/ ) {
#			$correct++;
#			if ( !defined($corByAnno{$anno}) ) { $corByAnno{$anno} = 1; } else { $corByAnno{$anno}++; }
#			next;
#		}
	}

	if ( !defined($annots{$anno}) ) { $annots{$anno} = 1; } else { $annots{$anno}++; }
	# Check annotation against prediction.
	if ( $anno eq $pred ) {
		if ( !defined($corByAnno{$anno}) ) { $corByAnno{$anno} = 1; } else { $corByAnno{$anno}++; }
		$correct++;
		next;
	} elsif ( $grouped && $pred =~ /^c-(.+)$/ ) {
		$pred = $1;
		if ( $anno =~ /^\Q$pred\E/ ) {
			if ( !defined($corByAnno{$anno}) ) { $corByAnno{$anno} = 1; } else { $corByAnno{$anno}++; }
			$correct++;
			next;
		}
	} elsif ( $pathList ) {
		if ( defined($pathArrayByAnnot{$anno}) ) {
			@alternatives = @{$pathArrayByAnnot{$anno}};
			$found = 0;
			foreach $alt ( @alternatives ) {
				if ( $pred eq $alt ) {
					if ( !defined($corByAnno{$alt}) ) { $corByAnno{$alt} = 1; } else { $corByAnno{$alt}++; }
					if ( !defined($annots{$alt}) ) { $annots{$alt} = 1; } else { $annots{$alt}++; }
					if ( $annots{$anno} == 1 ) { delete $annots{$anno} } else { $annots{$anno}--; }

					$correct++;
					$found = 1;
					last;
				}
			}

			if ( $found ) {
				next;
			}
		}
	}

	if ( !defined($suppress) ) {
		if ( $justSeqHeader ) {
			print $fields[0],"\n";
		} else {
			print $line,"\n";
		}
	}
}
close(IN);

if ( $pathList ) {
	foreach $anno ( keys(%annots) ) {
		if ( !defined($corByAnno{$anno}) && defined($pathArrayByAnnot{$anno}) ) {
			@alternatives = @{$pathArrayByAnnot{$anno}};
			foreach $alt ( @alternatives ) {
				if ( defined($corByAnno{$alt}) ) {
					$annots{$alt} += $annots{$anno};
					delete $annots{$anno};
					last;
				}
			}
		}
	}
}


if ( defined($terse) ) {
	print $correct,"\n";
} elsif (!$justSeqHeader) {
	print '--------------------------------------',"\n";
	print sprintf("%10s\t%5s\t%5s\t%6s\n",'CLADE','TP','Num','%Acc');
	print '--------------------------------------',"\n";
	foreach $anno (sort {$a cmp $b} keys(%annots)) {
		print sprintf("%10s\t%5d\t%5d\t%5.1f%%\n",$anno,$corByAnno{$anno},$annots{$anno},(100*$corByAnno{$anno}/$annots{$anno}));
	}
	print '--------------------------------------',"\n";
	print sprintf("%10s\t%5d\t%5d\t%5.1f%%\n",'TOTAL',$correct,$total,(100*$correct/$total));
	print '--------------------------------------',"\n";
}

# FNC - trim function.
# Removes whitespace from the start and end of the string
sub trim($) {
	my $string = shift;
	$string =~ /^\s*(.*?)\s*$/;
	return $1;
}
####################
