#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
03/23/2014 01:08:23 PM

=head1 DESCRIPTION

=head1 OPTIONS

=over 8

=item B<-a>

file 1

=item B<-b>

file 2

=item B<-out>

output

=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;

pod2usage(-message => "\n\tNo arguments\n", -verbose => 1) if (@ARGV == 0);

my $man  = 0;
my $help = 0;
my ($file_a, $file_b, $output);
GetOptions(
	'a=s' => \$file_a,
	'b=s' => \$file_b,
	'out=s' => \$output,
    'help|?' => \$help,
    'man|m'  => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage( -verbose => 2 ) if $man;

unless (-e $file_a && -e $file_b) {
    print STDERR "Input files not found\n";
    die;
}

my %mapping_a;
open my $fh, "<", $file_a or die $!;
while(<$fh>) {
	my ($query, $target, undef, undef, undef, undef, undef, $strand, $start, $end) =  split;
	next if ($mapping_a{$query});

    if ($strand && $start && $end) {
        $mapping_a{$query} = [$target, $strand, $start, $end];
    } else {
        $mapping_a{$query} = [$target];
    }

}
close $fh;

unless (%mapping_a) {
    die "No blast results or incompatible file: $file_a\n";
}

my %mapping_b;
open $fh, "<", $file_b or die $!;
while(<$fh>) {
	my ($query, $target) =  split "\t";
	next if ($mapping_b{$query});
	$mapping_b{$query} = $target;
}
close $fh;

unless (%mapping_b) {
    die "No blast results or incompatible file: $file_b\n";
}

## reverse mapping
my %back_b;
for (sort keys %mapping_b) {
	push @{$back_b{$mapping_b{$_}}}, $_;
}

## find one-to-one
my %result;
for (sort keys %mapping_a) {
    next unless (exists($back_b{$_}));

	my $value_a =  $mapping_a{$_}->[0];
	for my $v (@{$back_b{$_}}) {
		$result{$_} = $v if ($v eq $value_a)
	}
}

my $outfh;
if (defined $output) {
    open $outfh, ">", $output or die $!;
} else {
    $outfh = \*STDOUT;
}
for (sort keys %result) {
    if ($mapping_a{$_}->[2] && $mapping_a{$_}->[3]) {
        print $outfh $_."\t".$result{$_}."\t".$mapping_a{$_}->[1]."\t1\tbbh\t".$mapping_a{$_}->[2]."\t".$mapping_a{$_}->[3]."\n";
    } else {
        print $outfh $_."\t".$result{$_}."\n";
    }
}
close $outfh;

1;

__END__


