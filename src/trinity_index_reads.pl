#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
07/11/2014 11:10:38

=head1 DESCRIPTION

This is a simple, but helpful, script that stores "component - readfile"
relationship in a tab separated file. Recursivly indexing the
"chrysalis/Component_bins" directory. Absolute paths to readfile are stored.

Readfile contains reads which have been mapped by trinity to the inchworm
bundle, containing the largest number of common kmers.

=head1 OPTIONS

=over 8

=item B<-dir>

Base directory of trinity (which must contain intermediate files)

=item B<-output>

Output file of table (default: STDOUT)

=item B<-version>

Trinity changed the name pattern of contig ids and intermediate filename/paths.
Specify "c" for contig ids following name pattern like "c104698.graph_c0_seq1"
and "comp" for "comp27267_c0_seq1". [default: c]

=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;
use File::Glob;
use Storable;
use File::Spec;
use File::Path qw(make_path);

pod2usage(-message => "\n\tNo arguments\n", -verbose => 1) if (@ARGV == 0);

my $man  = 0;
my $help = 0;
my $dir;
my $output;
my $version = "c";
GetOptions(
    'dir=s'     => \$dir,
    'output=s'  => \$output,
    'help|?'    => \$help,
    'version=s' => \$version,
    'man|m'     => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

my ($subdir, $read_file_pattern, $read_path);

# Ya..trinity changed its contig IDs and storage system
if ($version eq "c") {
    $subdir = "Cbin";

    # c104698.graph_c0_seq1
    $read_file_pattern = '.*\.reads.tmp$';
    $read_path         = "/chrysalis/Component_bins";
} else {
    $subdir = "RawComps";

    # comp27267.raw.fasta
    $read_file_pattern = '.*\.raw.fasta$';
    $read_path         = "/chrysalis";
}

unless (defined $dir && -e $dir) {
    die "\nDirectory not found or not specified\n\n";
}

my $component_dir = File::Spec->catdir($dir, $read_path);
unless (-d $component_dir) {
    die "\nERROR:Intermediate output of Trinity not found!\n";
}

my $fh;
unless ($output) {
    $fh = \*STDOUT;
} else {
    my (undef, $dir, undef) = File::Spec->splitpath($output);
    make_path($dir);
    open $fh, ">", $output or die $!;
}

opendir my $dir_fh, $component_dir || die "Error opening: $component_dir\n";
my @subdir = grep {/$subdir/} readdir($dir_fh);
closedir $dir_fh;

for my $cbin (@subdir) {
    my $subdir = File::Spec->catdir($component_dir, $cbin);
    opendir $dir_fh, $subdir || die $!;
    for my $readfile (grep {/$read_file_pattern/} readdir($dir_fh)) {
        my ($comp) = split '\.', $readfile;
        print $fh "$comp\t" . File::Spec->catfile($subdir, $readfile) . "\n";
    }
    closedir $dir_fh;

}
close $fh;

1;

__END__

