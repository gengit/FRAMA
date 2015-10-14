#!/usr/bin/env perl

use strict;
use warnings;

=pod
=head1 INFO
Martin Bens, bensmartin@gmail.com
2014-07-24

=head1 DESCRIPTION
Stores table as storable object, which can be retrieved later. 

2 Variants:
Specify -c1 and -c2 => -c1 = key; -c2 = value; [one2many, array reference as value] 
Specify -c1 only => -c1 = key, complete line = value [one2many, array reference as value]

=head1 OPTIONS

=over 8
=item B<-input>

Input table

=item B<-output> 

Output file (default: input.storable)

=item B<-column1|c1> 

column which contains key 

=item B<-column2|c2> 

(optional) column which contains value

=item B<-noarray> 

one2one


=back

=cut

use Getopt::Long;
use Pod::Usage;
use Storable;
use Data::Dumper;

pod2usage(-message => "\n\tNo arguments\n", -verbose => 1) if (@ARGV == 0);

my $man  = 0;
my $help = 0;
my ($input, $col1, $col2, $output);
$col1 = 1;
$col2 = 0;
my $noarray = 0;
GetOptions(
    'input|i=s'    => \$input,
    'output|o=s'   => \$output,
    'column1|c1=i' => \$col1,
    'column2|c2=i' => \$col2,
    'help|?'       => \$help,
    'noarray' => \$noarray,
    'man|m'        => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if (!$input || !-e $input) {
    die "\nImport file not specified or defined\n\n";
}
if (-z $input) {
    die "Input file empty\n";
}

unless ($output) {
    $output = "$input.storable";
}

my %hash;

$col1--;
if ($col2) {

    # store column values in one to many hash
    $col2--;

    open my $fh, "<", $input or die $!;
    while (defined (my $line = <$fh>)) {
        chomp $line;
        my @e = split ' ', $line;
        if ($noarray) {
            $hash{$e[$col1]} = $e[$col2];
        } else {
            push @{$hash{$e[$col1]}}, $e[$col2];
        }

    }
    close $fh;

} else {

    # store complete line based. use column as index.
    my %hash;
    open my $fh, "<", $input or die $!;
    while (defined(my $line = <$fh>)) {
        chomp $line;
        my @e = split ' ', $line;
        if ($noarray) {
            $hash{$e[$col1]} = $line;
        } else {
            push @{$hash{$e[$col1]}}, $line;
        }
    }
    close $fh;
}

store \%hash, $output;

1;

__END__

