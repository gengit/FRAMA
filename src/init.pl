#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
2013-09-19

=head1 DESCRIPTION

Checks if necessary software and modules are installed.

=head1 OPTIONS

=over 8

=item B<-configuration>

Configuration file.

=back

=cut

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use Getopt::Long;
use FindBin;
use Data::Dumper;
use File::Spec;
use Pod::Usage;

my %mandatory = (
    "PATH_BAMTOOLS" => [ ["bamtools"] ],
    "PATH_BOWTIE" => [ ["bowtie2"], ["bowtie2-build"] ],

    # ubuntu binary different name: cdhit-est
    "PATH_EMBOSS" => [ ["needle"], ["getorf"] ],
    "PATH_MAFFT" => [ ["mafft"] ],
    "PATH_SAMTOOLS" => [ ["samtools"] ],

    # mac homebrew different name: Trinity
    "PATH_BLAST" => [ [ "xdformat", "makeblastdb" ], ["blastn"] ],

    # look for wublast or ncbiblast
    "PATH_PERL"    => [ ["perl"] ],
    "PATH_RSCRIPT" => [ ["Rscript"] ],
    "PATH_GENSCAN" => [ ["genscan"] ],
    "PATH_TRINITY" => [ [ "Trinity", "Trinity.pl" ] ],
);

my %optional = (
    "PATH_TGICL"        => [ ["tgicl"] ],
    "PATH_CD_HIT_EST"   => [ [ "cd-hit-est", "cdhit-est" ] ],
    "PATH_REPEATMASKER" => [ ["RepeatMasker"] ],
);

my %executables = ( %mandatory, %optional );
my %dep1 = map { $_ => 1 } keys %mandatory;
my %dep2 = map { $_ => 0 } keys %optional;
my %dependence = ( %dep1, %dep2 );

die "Need configuration file\n" unless ( @ARGV > 0 );
my $file = $ARGV[0];

my $missing = 0;

# check perl modules --> TODO: cpanfile
my @output = `grep "^use " $FindBin::Bin/*.pl`;
my @modules =
  grep { !/constant/ && !/lib/ && !/warnings;/ && !/strict;/ } @output;
for (@modules) { /.+:use\s(.+?)(\s|;)/; $_ = $1; }
my %modules = map { $_ => 1 } @modules;

for ( sort keys %modules ) {
    next if (/^\d/);
    eval "require $_";
    if ($@) {
        print "FAILED: $_\n";
        $missing = 1;
    }
}

# check software
my %visited = map { $_ => 0 } keys %executables;

#open my $fh, "<", $file or die $!;
#while (<$fh>) {
#    chomp;
#
#    next if (/^#/);
#
#    if (/^(.+?)\s*:=\s*(.+?)$/) {
#        my $variable = $1;
#        my $value    = $2;
#
#        if ( defined $value ) {
#            $value =~ s/\s+//g;
#            if ( length($value) == 0 ) {
#                $value = 0;
#            } else {
#                $value = File::Spec->rel2abs($value);
#            }
#        }
#
#        $visited{$variable} = $value;
#    }
#}
#close $fh;

while ( my ( $variable, $value ) = each(%visited) ) {
    if ( $variable =~ /^REF_/ ) {
        $missing += checkFile($value);
    } elsif ( $variable =~ /^PATH_GENSCAN_MAT/ ) {
        $missing += checkFile($value);
    } elsif ( $variable =~ /^ORTHOLOG_/ ) {
        $missing += checkFile($value);
    } elsif ( $variable =~ /^PATH_/ ) {
        $missing += checkExecutable( $variable, $value );
    }
}

die "\nPlease check mandatory required software, modules and files!\n\n"
  if ($missing);

sub checkFile {
    my ($file) = @_;
    $file =~ s/\s+$//g;

    unless ( -e $file ) {

        print "FAILED: $file not found\n";
        return 1;
    }

    return 0;
}

sub checkExecutable {
    my ( $variable, $path ) = @_;

    # ignore unknown variables but give warning
    unless ( exists $executables{$variable} ) {
        print "FAILED: [unknown variable] $variable\n";
        return 0;

        #$missing = 1; next;
    }

    my $programs = $executables{$variable};

    for my $program_options (@$programs) {

        my $failed = 0;
        for my $program (@$program_options) {
            my $exe;
            if (length($program) == 0) {
                $exe = $path
            } else {
                $exe = $program;
                $exe = File::Spec->catfile( $path, $program ) if ($path);
            }

            unless ( executable($exe) ) {
                $failed++;
            }
        }

        if ( $failed == scalar @$program_options ) {
            if ( $dependence{$variable} ) {
                print "FAILED: [mandatory] $variable\n";
            } else {
                print "FAILED: [optional] $variable\n";
            }
            print "\tPATH: $path\n";
            return ( $dependence{$variable} );
        }
    }
    return 0;
}

sub executable {
    my ($command) = @_;
    my $failed = `command -v "$command"`;
    chomp($failed);

    # command -v returns path if found
    if ($failed) {
        return 1;
    }
    return 0;
}

1;

__END__





