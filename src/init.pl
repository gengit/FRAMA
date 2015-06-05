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

use Getopt::Long;
use FindBin;
use Data::Dumper;
use File::Spec;
use Pod::Usage;

use lib "$FindBin::Bin/../lib/perl";

my %executables = (
    "PATH_BAMTOOLS" => "bamtools",
    "PATH_BOWTIE" => ["bowtie2", "bowtie2-build"],
    "PATH_CD_HIT_EST" => "cd-hit-est",
    "PATH_EMBOSS" => ["getorf", "needle"],
    "PATH_GENSCAN" => "genscan",
    "PATH_MAFFT" => "mafft",
    "PATH_PERL" => "perl",
    "PATH_REPEATMASKER" => "RepeatMasker",
    "PATH_RSCRIPT" => "Rscript",
    "PATH_SAMTOOLS" => "samtools",
    "PATH_TGICL" => "tgicl",
    "PATH_TRINITY" => "Trinity.pl",
    "PATH_WUBLAST" => "blastn",
    "PATH_XDFORMAT" => "xdformat",
    "PATH_GENSCAN_MAT" => 0,
    "PATH_BLAST" => "blastn"
);

die "Need configuration file\n" unless (@ARGV > 0);
my $file = $ARGV[0];

my $missing = 0;

# check perl modules
my @output = `grep "^use " $FindBin::Bin/*.pl`;
my @modules =  grep {!/constant/ && !/lib/ && !/warnings;/ && !/strict;/} @output;
for (@modules) { /.+:use\s(.+?)(\s|;)/; $_ = $1; }
my %modules = map { $_ => 1 } @modules;

for ( sort keys %modules ) {
    next if (/^\d/);
    eval "require $_";
    if ($@) { print "FAILED: $_\n"; $missing = 1;
    } else { print "OK: $_\n"; }
}

# check software
my %visited = map { $_ => 0 } keys %executables;
my %paths_to_executables; 

my %paths;
my $stop = 0;
open my $fh, "<", $file or die $!;
while (<$fh>) {
    chomp;

    next if (/^#/);

    if (/^(PATH.+?)\s+?:=\s+?(.+?)$/) {
        my $variable = $1;
        $visited{$variable} = 1;
        unless (exists $executables{$variable}) {
            print "FAILED: Unkown variable: $variable.\n";
            $missing = 1; next;
        }

        my $path = $2;
        $path =~ s/\s+//g;

        my $exe = $executables{$variable};
        if (length($exe) == 1) {
            if (-e $path) {
                print "OK: $variable\n";
            } else {
                print "FAILED: $variable\n";
                $missing = 1; 
            }
            next;
        }
        $path .= "/" unless ($path =~ /\/$/);

        if (ref($exe) eq "ARRAY") {
            for my $e (@$exe) {
                my $executable = File::Spec->catfile($path, $e) if (length($path));
                test($executable);
                $paths_to_executables{$variable} = $executable;
            }
        } else {
            my $executable = File::Spec->catfile($path, $exe) if (length($path));
            test($executable);
            $paths_to_executables{$variable} = $executable;
        }
    }
    if (/^(REF.+)\s*:=\s*(.+?)$/) {
        my $file = $2;
        $file =~ s/\s+//g;
        next unless (length($file) > 0);
        if (-e $file) {
            print "OK: $file\n";
        } else {
            print "FAILED: $file not found\n";
            $missing = 1;
        }
    }
    if (/^ORTHOLOG_TABLE\s*:=\s*(.+?)$/) {
        my $file = $1;
        $file =~ s/\s+//g;
        next unless (length($file) > 0);
        if (-e $file) {
            print "OK: $file\n";
        } else {
            print "FAILED: $file not found\n";
            $missing = 1; 
        }
    }
    if (/^ORTHOLOG_CDS\s*:=\s*(.+?)$/) {
        my $file = $1;
        $file =~ s/\s+//g;
        next unless (length($file) > 0);
        if (-e $file) {
            print "OK: $file\n";
        } else {
            print "FAILED: $file not found\n";
            $missing = 1; 
        }
    }
}
close $fh;

for (keys %visited) {
    next if  ($visited{$_}); 
    if ((ref $executables{$_}) =~ /ARRAY/) {
        for my $executable (@{$executables{$_}}) { test($executable); $paths_to_executables{$_} = $executable}
    } else {
        test($executables{$_});
        $paths_to_executables{$_} = $executables{$_};
    }
}

die "\nPlease install required software and modules!\n\n" if ($missing); 

sub test {
    my $command = shift;
    unless ($command) {
        print "Not going to do that..\n";
        return;
    }
    my $failed = `command -v "$command"`;
    chomp($failed);
    unless ($failed) {
        print "FAILED: $command $failed\n";
        $missing = 1;
    } else {
        print "OK: $command\n";
    }
}

1;

__END__

