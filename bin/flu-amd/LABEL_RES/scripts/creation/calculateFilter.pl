#!/usr/bin/env perl

use POSIX;

if ( scalar(@ARGV) != 2 ) {
	$message = "\nUsage:\n\tperl $0 <table> <true_positive_label>\n";
	die($message."\n");
}


open(TBL,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
$/ = "\n"; $positive_class = uc($ARGV[1]);
$header = <TBL>; chomp($header);
@headers = split("\t",$header);


$IDf = $annoF = $lenF = -1;
%hmm = (); @hmmF = ();
foreach $index ( 0..$#headers ) {
	$header = uc($headers[$index]);
	if ( $header =~ /HMM/ ) {
		$hmm{$header} = $index;
	} elsif ( $header eq 'LENGTH' ) {
		$lenF = $index;
	} elsif ( $header eq 'ANNOTATION' ) {
		$annoF = $index;
	} elsif ( $header eq 'ID' ) {
		$IDf = $index;
	}
}

if ( $IDf < 0 || $lenF < 0 || $annoF < 0 ) { die("Missing field.\n"); }
@hmmF = values(%hmm); %frags = ();
if ( scalar(@hmmF) < 1 ) { die("Missing hmm field.\n"); }
%max = %min = (); $maxL = 0;
while($line = <TBL>) {
	chomp($line);
	@f = split("\t",$line);
	($ID,$length,$annotation) = ($f[$IDf],$f[$lenF],$f[$annoF]);

	$score = $f[$hmmF[0]] / $length;
	for $i ( @hmmF ) {
		$tmp = $f[$i] / $length;
		if ( $tmp < $score ) {
			$score = $tmp;
		}
	}

	if ( uc($annotation) eq $positive_class ) { 
		$pos=1;
		if ( $maxL < $length ) { 
			$maxL = $length;
		}
	} else {
		$pos=0;
	}

	if ( $ID =~ /\Q..\E/ ) {
		if ( $pos == 0 ) { $frags{$ID} = $score; }
	} else {
		if ( $score < $min{$pos} || !defined($min{$pos}) ) { $min{$pos} = $score; }
		if ( $score > $max{$pos} || !defined($max{$pos}) ) { $max{$pos} = $score; }
	}
}
close(TBL);

$POS = 1; $NEG = 0;
$t1 = mid($max{$POS},$min{$NEG});
@scores = sort { $a <=> $b } values(%frags);
$N = scalar(@scores); $k0 = floor($N/1000); $k1 = ceil($N/1000);
$t2 = mid($scores[$k0],$scores[$k1]);
$maxL = ceil($maxL/100) * 100;


if ( $t1 < $t2 ) {
	$t3 = $t1;
} else {
	$t3 = $t2;
}

print STDERR sprintf("[%.3f,%.3f]+ [%.3f,%.3f]-\n",$min{$POS},$max{$POS},$min{$NEG},$max{$NEG});
print STDERR sprintf("min(%.3f,%.3f) = %.3f\n",$t1,$t2,$t3);
print sprintf("%.3f %d\n",$t3,$maxL);

sub mid($$) {
	my ($x1,$x2) = @_;
	if ( $x1 < $x2 ) {
		return ( $x1 + ($x2-$x1)/2  );
	} else {
		return ( $x2 + ($x1-$x2)/2  );
	}
}
