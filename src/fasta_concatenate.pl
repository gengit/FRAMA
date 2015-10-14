#!/usr/bin/env perl 

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
06/10/2014 17:04:41

=head1 DESCRIPTION

Concatenates multi-fasta files into a single long sequence. ID + Start + End of
every single sequence is stored in output.index (1based coordinates).

=head1 OPTIONS

=over 8

=item B<-file> 

Multi-Fasta file.

=item B<-output> 

Output stem for output files (index and fa). [default: output] 

=item B<-space> 

Separate sequences by sequence of characters. 

=back

=cut

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use Bio::SeqIO;

pod2usage( -message => "\n\tNo arguments\n", -verbose => 1 ) if ( @ARGV == 0 );

my $man  = 0;
my $help = 0;
my ( $file, $output );
my $space = "";
GetOptions(
    'file|f=s'   => \$file,
    'output|o=s' => \$output,
    'spacer|s=s' => \$space,
    'help|?'     => \$help,
    'man|m'      => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage( -verbose => 2 ) if $man;

if ( -z $file ) {
    die "\n\tFile not found or empty: $file\n\n";
}

unless ($output) {
    die "\n\tSpecify output!\n\n";
}

my $io = Bio::SeqIO->new( -file => $file, -format => 'fasta' );

my @info;
my $start = 1;
my @seqs;
my $length_space = length($space);
while ( my $seq = $io->next_seq ) {
    my $end = $start + length( $seq->seq ) - 1;
    if ( $seq->desc ) {
        push @info, [ $seq->id . " " . $seq->desc, $start, $end + $length_space ];
    } else {
        push @info, [ $seq->id, $start, $end ];
    }

    $start = $end + 1 + $length_space;
    push @seqs, $seq->seq;
}
open my $out_fasta, ">", $output . ".fa" or die "Can't open: $output.fa\n";
print $out_fasta ">combined\n";
toFasta( undef, join( $space, @seqs ), undef, $out_fasta );
close $out_fasta;

open my $out_index, ">", $output . ".index" or die "Can't open: $output.index\n";
for (@info) {
    print $out_index "$_->[0]\t$_->[1]\t$_->[2]\n";
}
close $out_index;

sub toFasta {
    my ( $seqName, $seq, $len, $fh ) = @_;

    $fh  = \*STDOUT unless ($fh);
    $len = 60       unless $len;

    print $fh ">$seqName\n" if ( defined $seqName );
    while ( my $chunk = substr( $seq, 0, $len, "" ) ) {
        print $fh $chunk . "\n";
    }
}

1;

__END__

