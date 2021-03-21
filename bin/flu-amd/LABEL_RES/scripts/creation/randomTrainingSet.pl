#!/usr/bin/env perl
# Samuel Shepard - 12.2012
# Version 1.0

use Getopt::Long;
GetOptions(
		'appendix|A=s' => \$appendixFile,
		'samples|S=s' => \$setSizes, 'prefix|P=s' => \$prefix, 
       		'with-replacement|W' => \$withReplacement,
       		'allow-upsampling|U' => \$allowUpsampling,
	       	'distinct-names|D' => \$distinctNames,
		'verbose|V' => \$verbose,
		'exhaustive|E' => \$exhaustive,
		'ID-file|I' => \$IDfile,
       		'total-coverage|T' => \$totalCoverage	);

if (  scalar(@ARGV) < 2 || !($setSizes xor $exhaustive) ) {
	$message = "\n$0 <data_matrix.tab> <output_dir> {-S <#samples>|-E} [OPTIONS]\n";
       	$message .= "\t-S|--samples #[,#,...,#]\t\tComma-delimited list of number of samples per file. Rotates.\n";
	$message .= "\t-W|--with-replacement\t\t\tSampling with replacement.\n";
	$message .= "\t-U|--allow-upsampling\t\t\tUp-samples if set too small, implies -W option.\n";
	$message .= "\t-D|--distinct-names\t\t\tMakes distinct headers for up-sampled data, implies -U option.\n";
	$message .= "\t-T|--total-coverage\t\t\tEnsure total coverage for up-sampled data, implies -U option.\n";
	$message .= "\t-A|--appendix <file>\t\t\tEnsure IDs listed in appendix are selected.\n";
	$message .= "\t-P|--prefix\t\t\t\tPrefix to output files.\n";
	$message .= "\t-E|--exhaustive\t\t\t\tDo not sample, merely produce a set of all data.\n";
	$message .= "\t-I|--ID-file\t\t\t\tOutput an ID file rather than an INFO file.\n";
	die($message."\n");
}

$delim = '|';
if ( defined($distinctNames) ) {
	$allowUpsampling = 1;
}

if ( defined($allowUpsampling) ) {
	$withReplacement = 1;
}

if ( defined($totalCoverage) ) {
	$withReplacement = 0;
	$allowUpsampling = 1;
}

if ( defined($appendixFile) ) {
	if ( open(APN,'<', $appendixFile ) ) {
		$/ = "\n"; %appendix = ();
		while($line = <APN>){
			chomp($line);
			if (!defined($appendix{$line})) {
				$appendix{$line} = 1;
			} else {
				$allowUpsampling = 1;
				$withReplacement = 1;
				$appendix{$line}++;
			}
		}
		close(APN);
		$appended = 1;
	} else {
		$appended = 0;
	}
}

# process parameters
open( IN, '<', $ARGV[0] ) or die("ERROR: Could not open $ARGV[0].\n");
$outdir = $ARGV[1];
if ( defined($prefix) ) {
	$prefix .= '_';
}


open( DAT, '>', $outdir."/${prefix}training.dat") or die("$0 ERROR: Cannot open training.dat\n");
if ( $IDfile ) { 
	$suffix = 'IDs';
} else {
	$suffix = 'info';
} 
open( INF, '>', $outdir."/${prefix}${suffix}.dat") or die("$0 ERROR: Cannot open info.dat\n");
open( LAB, '>', $outdir."/${prefix}labels.dat") or die("$0 ERROR: Cannot open labels.dat\n");

@samplesetSizes = split(',', $setSizes);
%usedRecords = (); @headers = (); @sequences = ();
%scores = %labelsByTaxa = %idByTaxa = %countByTaxa = ();
$header = <IN>; $id = $taxon = '';
while( $line = <IN> ) {
	chomp($line);
	@values = split("\t", $line);
	$id = shift(@values);
	$taxon = shift(@values);
	$scores{$id} = [@values];

	$labelsByTaxa{$taxon} = 1;
	if ( defined($countByTaxa{$taxon}) ) {
		$idByTaxa{$taxon}[$countByTaxa{$taxon}] = $id;
		$countByTaxa{$taxon}++;
	} else {
		$idByTaxa{$taxon}[0] = $id;
		$countByTaxa{$taxon} = 1;
	}
	$usedRecords{$id} = 0;
}
close(IN);
$totalRequested = 0;
$Ns = scalar(@samplesetSizes);


@taxa = sort(keys(%labelsByTaxa));
$Nt = scalar(@taxa);
if ( $Nt > 2 ) {
	for($label = 0; $label < $Nt; $label++ ) {
		$labelsByTaxa{$taxa[$label]} = $label;
		$taxaByLabels{$label} = $taxa[$label];
	}
} else {
	$labelsByTaxa{$taxa[0]} = -1;
	$taxaByLabels{'-1'} = $taxa[0];
	$labelsByTaxa{$taxa[1]} = 1;
	$taxaByLabels{'1'} = $taxa[1];
}
@labels = sort {$a <=> $b} keys(%taxaByLabels);

if ( defined($verbose) ) {
	print "\nNumber\t Have\t Want\t  Got\tTaxon\n";
	print '----------------------------------------',"\n";
}
if ( $exhaustive ) {
	for($i = 0; $i < $Nt; $i++ ) {	
		$taxon = $taxa[$i];
		@taxonIDs = @{$idByTaxa{$taxon}};
		foreach $theID (@taxonIDs) {
			$labelsByID{$theID} = $labelsByTaxa{$taxon};
			$usedRecords{$theID} = 1;
			$sampleCount++;
		}
	}
} else {
	for($i = 0; $i < $Nt; $i++ ) {	
		$taxon = $taxa[$i];
		$setSize = $samplesetSizes[$i%$Ns];
		$have = $countByTaxa{$taxon};
		$sampleCount = 0;
		@taxonIDs = @{$idByTaxa{$taxon}};
		$Ng = scalar(@taxonIDs);
		$recordsAvailable = $recordsAvailable2 = $Ng;

		if ( $appended ) {
			foreach $theID (@taxonIDs) {
				if ( defined($appendix{$theID}) ) {
					while ( $appendix{$theID} > 0 ) {
						$labelsByID{$theID} = $labelsByTaxa{$taxon};
						$usedRecords{$theID}++;
						$sampleCount++;
						if ( !$withReplacement ) {
							$recordsAvailable--;
						}
						$appendix{$theID}--;
					}
				}
			}	
		}

		while( $sampleCount < $setSize ) {
			$slot 	= int(rand($Ng));
			$theID	= $taxonIDs[$slot];

			if ( $usedRecords{$theID} == 0 || $withReplacement ) {
				$labelsByID{$theID} = $labelsByTaxa{$taxon};
				$usedRecords{$theID}++;
				$sampleCount++;
				if ( !$withReplacement ) {
					$recordsAvailable--;
				}
			}

			if ( $recordsAvailable == 0 ) {
				if ( defined($totalCoverage) ) {
					$recordsAvailable = $recordsAvailable2;
					$withReplacement = 1;
				} else {
					last;
				}
			}
		}

		if ( defined($verbose) ) {
			printf("%6d\t%5d\t%5d\t%5d\t%s\n",($i+1),$have,$setSize,$sampleCount,$taxon);
		}
	}
}
if ( defined($verbose) ) {
	print "\n";
}

if ( !$IDfile ) {
	for( $i = 0; $i < $Nt; $i++ ) {
		print INF $labelsByTaxa{$taxa[$i]},"\t",$taxaByLabels{$labelsByTaxa{$taxa[$i]}},"\n";
	}
	print INF "\n";
}

@sortedIDs = sort(keys(%usedRecords));
$k = 1; @labelBuffer = ();
foreach $id ( @sortedIDs ) {
	$count = $usedRecords{$id};
	if ( $count > 0 ) {
		if ( defined($distinctNames) ) {
			while ( $count > 0 ) {
				$labelBuffer[$k-1] = $labelsByID{$id};
				print INF $k++,"\t",$count,$delim,$id,"\n";
				print DAT join(' ',@{$scores{$id}}),"\n";
				$count--;
			}
		} else {
			$labelBuffer[$k-1] = $labelsByID{$id};
			print INF $k++,"\t",$id,"\n";
			print DAT join(' ',@{$scores{$id}}),"\n";
		}
	}
}
print LAB join(' ', @labelBuffer),"\n";
close(DAT);
close(INF);
close(LAB);
