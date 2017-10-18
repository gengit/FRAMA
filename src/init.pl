#!/usr/bin/env perl

# --------------------------------------------------------------------
# Check required tools
#
# Martin Bens | bensmartin@gmail.com
# 2017-10-18
# --------------------------------------------------------------------

use strict;

use v5.10;

use YAML;
use Data::Dumper;
use File::Spec;
use File::Copy;

my $execute_dir  = "bin";

sub getFinal {
    my ($data) = @_;

    my @target  = split " ", $data->{executable};
    my @final = ();
    for (my $i = 0; $i < @target; $i++) {
        push @final, File::Spec->catfile($execute_dir, $target[$i]);
    }
    return @final;
}

sub checkProgram {
    my ($data) = @_;

    say "- $data->{description}";

    my @programs = split " ", $data->{executable};
    my $installed = 0;
    if (@programs) {
        for my $e (@programs) {
            next if ($e eq "mafftdir" | $e eq "bin");
            my $command = "which $e 2> /dev/null";
            my ($output) = `$command`;
            if (defined $output) {
                chomp $output;
                say "\t$e installed: $output";
                $installed++;
                symlink($output, File::Spec->catfile($execute_dir, $e));
            } else {
                say "\t$e needs to be installed";
            }
        }
    }
    say "";

    $data->{installed} = 0;
    if ($installed == @programs) {
        $data->{installed} = 1;
    }
}

my $yaml = YAML::LoadFile($ARGV[0]);

# Check installation
foreach my $program (@$yaml) {
    checkProgram($program);
}

1;

# THE END
