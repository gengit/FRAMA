#!/usr/bin/env perl

# ACE to Table

use strict;
use warnings;

my $file = $ARGV[0];

my $current_cluster;
my @members;

open my $fh, "<", $file or die "Can't open file for reading: $file\n";
while(my $line = <$fh>) {

    if ($line =~ /^CO\s(.+?)\s/) {
        if (@members > 0) {
            print "$current_cluster\t".join ",",@members;
            print "\n";
        }
        @members = ();
        $current_cluster = $1;
    }
    if ($line =~ /^AF\s(.+?)\s/) {
        push @members, $1;
    }
}
close $fh;

if (@members > 0) {
    print "$current_cluster\t".join ",",@members;
    print "\n";
}

1;
