#!/usr/bin/env perl

use strict;
use warnings;

my $genblasta = $ARGV[0];
my $blast = $ARGV[1];

my $evalue = 10;
if ($blast =~ /-e\s+(.+?)(\s|$)/) {
    $evalue = $1;
} elsif ($genblasta =~ /-e\s+(.+?)(\s|$)/) {
    $evalue = $1;
}

$blast =~ s/-e\s+.+?(\s|$)//g;
$genblasta =~ s/-e\s+.+?(\s|$)//g;

if (length($blast) > 0) {
    print "$genblasta -e \"$evalue $blast\"\n";
}

1;
