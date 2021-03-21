#!/usr/bin/env perl
# Extract Random Fragments - Sam Shepard - 9.2012


use Getopt::Long;
GetOptions( 'fragment-length|L=i'=> \$fragLength, 'gap-flanking|G' => \$gapFlanking ); #, 'number-sequences|N=i' => \$numSequences );
if ( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\t$0 <input.fa> [OPTIONS]\n";
	$message .= "\t\t-L|--fragment-length <#>\tLength of fragments to take.\n";
	$message .= "\t\t-G|--gap-flanking\t\tGap the flanking regions to fill alignment.\n";
#	$message .= "\t\t-N|--number-sequences <#>\tRotate over N sequences instead of matching input set.\n";
	die($message);
}

if ( !defined($fragLength) ) {
	$fragLength = 50;
}


$/ = '>';
open(IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0].\n");
while( $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = lc(join('',@lines));

	$length = length($sequence);
	if ( $length == 0 || $length < $fragLength ) {
		next;	
	}


	$N = $length - $fragLength + 1;
	$offset = int(rand($N));
	if ( defined($gapFlanking) ) {
		# number to pad on the 5 prime end, starting at 0
		$leftPad = $offset;
		# number to pad on the 3 prime end, starting at offset + fraglength
		$rightPad = $length - $offset - $fragLength;
		if ( $leftPad > 0 ) {
			substr( $sequence, 0, $leftPad, '-' x $leftPad );
		}

		if ( $rightPad > 0 ) {
			# create BUFFER or, find a replacement method
			substr( $sequence, ($offset+$fragLength), $rightPad, '-' x $rightPad );
		}
		print '>',$header,"\n",$sequence,"\n";
	} else {
		$fragment = substr($sequence, $offset, $fragLength);
		print '>',$header,"\n",$fragment,"\n";
	}
}
close(IN);
