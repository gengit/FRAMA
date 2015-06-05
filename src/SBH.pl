#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;
use Bio::SearchIO;

use constant QUERY_ID  => 0;
use constant TARGET_ID => 1;

my $man  = 0;
my $help = 0;
my $file_blast;
my $file_annotation;
GetOptions(
    'blast=s'     => \$file_blast,
    'annotated=s' => \$file_annotation,
    'help|?'      => \$help,
    man           => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

my %known_contigs;
my %known_references;
open my $fh, "<", $file_annotation or die $!;
while (<$fh>) {
    my @e = split "\t";
    $known_contigs{$e[0]}    = 1;
    $known_references{$e[1]} = 1;
}
close $fh;

my $blastfh = \*STDIN;
if ($file_blast) {
    open $blastfh, "<", $file_blast or die "Can't open file: $file_blast\n";
}

while (<$blastfh>) {
    my ($query, $target) = split "\t";

    # already annotated/added reference or contig
    next if (exists $known_references{$target});
    next if (exists $known_contigs{$query});

    # SBH found
    print;

    # adding SBH to known contigs
    $known_references{$target} = 1;
    $known_contigs{$query}   = 1;
}
close $blastfh;

1;

__END__

