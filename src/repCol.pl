#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
10/22/2013 15:03:58

=head1 DESCRIPTION

Replaces content of column in Table1 by assignments in Table2.

Example:

    Table1       Table2      Result
    | A | 1 |    | A | Z |   | Z | 1 |
    | B | 2 |    | B | Y |   | Y | 2 |
    | B | 3 |    | C | X |   | Y | 3 |
    | C | 4 |                | X | 4 |

=head1 OPTIONS

=over 8

=item B<-table1> tab-sep

Table with column to be replaced.

=item B<-table2> tab-sep

Table with replacements.

=item B<-column> num

Select column to be replaced in table1 by number (1-based).

=item B<-tkey> num

Select column in table2 containing the IDs corresponding to IDs in table2 by number (1-based).

=item B<-tvalue> num

Select column in table2 with replacements by number (1-based).

=item B<-add>

Adds column to table1 instead of replacing values.

=item B<-NA>

Use "NA" if no replacement was found (otherwise keeps value).

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
my ($table1, $table2);
my ($col1, $tkey, $tvalue) = (1, 1);
my $add = 0;
my $na  = 0;

GetOptions(
    'table1|t1=s' => \$table1,
    'table2|t2=s' => \$table2,
    'column|c=i'  => \$col1,
    'tkey=i'      => \$tkey,
    'tvalue=i'    => \$tvalue,
    'na'          => \$na,
    'add'         => \$add,
    'help|?'      => \$help,
    man           => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

$col1--; $tkey--; $tvalue--;

my %mapping;
open my $fh, "<", $table2 or die "Can't open file: $table2\n";
while (<$fh>) {
    next if (/^#/);
    chomp;
    my @e = split;
    $mapping{$e[$tkey]} = $e[$tvalue];
}
close $fh;

my $c_block = "";
my $c_replace;
open $fh, "<", $table1 or die "Can't open file: $table1\n";
while (<$fh>) {
    next if (/^#/);
    chomp;
    my @e = split;
    next unless ($#e >= $col1);
    if ($e[$col1] ne $c_block) {
        $c_replace = $mapping{$e[$col1]};
        if ($na && not defined $c_replace) {
            $c_replace = "NA";
        } elsif (not defined $c_replace) {
            $c_replace = $e[$col1];
            print STDERR "Replacement for $e[$col1] not found!\n";
        }

        #$c_replace = $e[$col1] unless ($c_replace);
        $e[$col1] = $c_replace if (!$add);
        push @e, $c_replace if ($add);
    } else {
        $e[$col1] = $c_replace if (!$add);
        push @e, $c_replace if ($add);
    }
    print join("\t", @e) . "\n";
}
close $fh;

1;

__END__

