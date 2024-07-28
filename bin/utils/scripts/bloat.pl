#!/usr/bin/perl

open(FILE,$ARGV[0]) ;
while (<FILE>)
{ @values = split;
  $values[3] += $ARGV[1];
  print "$values[0] $values[1] $values[2] $values[3]\n";
}
close(FILE)
  
