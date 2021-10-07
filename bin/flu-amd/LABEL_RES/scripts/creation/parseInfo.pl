#!/usr/bin/perl


$/ = "\n\n";
open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0].\n");
$record = <IN>;
$record = <IN>;
close(IN);

$extra = '';
foreach $arg ( @ARGV[1..$#ARGV] ) {
	$extra .= "\t$arg";
}

@lines = split("\n",$record);
foreach $line ( @lines ) {
	($label,$ID) = split("\t",$line);
	@header = split('\|',$ID);
	$ID = join('|',@header[1..$#header]);
	print $ID,$extra,"\n";
	
}
