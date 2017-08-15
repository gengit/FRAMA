#!/usr/bin/env perl

use strict;
use warnings;

use v5.10;

use YAML;
use Data::Dumper;
use File::Spec;

my $download_dir = "download";
my $execute_dir  = "bin_external";

sub checkProgram {
    my ($data) = @_;

    my @programs = split " ", $data->{executable};
    if (@programs) {
        for my $e (@programs) {
            my $command = "which $e";
            my ($output) = `$command`;
            if (defined $output) {
                chomp $output;
                say "$data->{description} is installed: $output";
                $data->{installed} = $output;
            } else {
                say "$data->{description} needs to be installed";
                $data->{installed} = 0;
            }
        }
    }
}

sub download {
    my ($h) = @_;

    my $filename = "$download_dir/$h->{description}";
    $filename .= "." . $h->{extension} if (exists $h->{extension});

    my $wget = "wget -O $filename -o $filename.log";
    $wget .= " " . $h->{wget} if (exists $h->{wget});

    $wget .= " " . $h->{remote};

    if (-e $filename) {
        say "File already exists: $filename. Skipping.";
        $h->{filename} = $filename;
        return 1;
    }

    say "Running:\n\t$wget";
    my $error = system($wget);
    if (!$error) {
        $h->{filename} = $filename;
    } else {
        say "Download failed: $h->{description} from $h->{remote}";
        return 0;
    }

    return 1;
}

sub unpackMe {
    my ($h) = @_;

    # nothing to do
    if (not exists $h->{dir}) {
        return 1;
    }

    # gz bz2 zip
    my $command = "";
    if ($h->{extension} eq "tar.gz") {
        $command = "tar zxf";
    } elsif ($h->{extension} eq "tar.bz2") {
        $command = "tar jxf";
    } elsif ($h->{extension} eq "zip") {
        $command = "unzip -qq -d $download_dir";
    }
    $command .= " " . $h->{filename};

    # custom directory
    if ($command =~ /^tar/) {
        $command .= " -C $download_dir";
    }
    my $mydir = "$download_dir/$h->{dir}";

    if (-e $mydir) {
        say "File already unpacked: $mydir. Skipping";
        $h->{unpacked} = $mydir;
        return 1;
    }

    say "Running:\n\t$command";

    my $error = system($command);
    if (!$error) {
        $h->{unpacked} = $mydir;
    } else {
        say "Unpacking failed: $h->{description} from $h->{filename}";
    }
}

sub compileMe {
    my ($h) = @_;

    if (exists $h->{compile} && exists $h->{dir}) {
        my $command_file = "$h->{unpacked}/FRAMA_run.sh";
        say "Write commands; $command_file";
        open my $fh, ">", $command_file || die;
        print $fh $h->{compile};
        close $fh;

        say "Running $command_file";
        my $command = "cd $h->{unpacked} && bash FRAMA_run.sh";
        say "Running:\n\t$command\n";
        my $error = system($command);
        if (!$error) {
            $h->{compilation} = 1;
        } else {
            say "Installation failed: $h->{description} from $command_file"
        }
    }
}

sub linkfiles {
    my ($data) = @_;

    my @sources = split " ", $data->{links};
    my @target  = split " ", $data->{executable};
    for (my $i = 0; $i < @sources; $i++) {
        my $source = File::Spec->catfile($download_dir, $sources[$i]);
        $source = File::Spec->catfile($download_dir, $data->{dir}, $sources[$i])
            if (exists $data->{dir});

        $source = File::Spec->rel2abs($source);
        chmod(0755, $source);
        my $target = File::Spec->catfile($execute_dir, $target[$i]);
        symlink($source, $target);

        say "Executable added: $target";
    }
}

sub installRM {
    my ($data) = @_;

    # ncbi-blast
    #compileMe($data->{Blast});
    linkfiles($data->{Blast});

    linkfiles($data->{rmblast});
    linkfiles($data->{trf});

    # RepeatMasker
    my $trf_bin = File::Spec->rel2abs("$execute_dir/trf");
    my $exed = File::Spec->rel2abs("$execute_dir");
    my $command = "cd $download_dir/$data->{repeatmasker}->{dir} && printf \"\\n$^X\\n\\n$trf_bin\\n2\\n$exed\\nY\\n5\\n\" | perl configure";
    my $error = system($command);
    linkfiles($data->{repeatmasker});
}

my $yaml = YAML::LoadFile('programs.yaml');

mkdir $download_dir;
mkdir $execute_dir;

# Check installation
foreach my $program (@$yaml) {
    checkProgram($program);
}
print Dumper $yaml;
exit;

# Download files
foreach my $program (@$yaml) {
    download($program) unless ($program->{installed});
}
# Unpack files
foreach my $program (@$yaml) {
    unpackMe($program) if ($program->{filename});
}

# Install RepeatMasker
# repeatmasker depence on blast, rmblast, trf
my $repeatmasker;
foreach my $program (@$yaml) {
    if ($program->{description} eq "repeatmasker") {
        $repeatmasker->{repeatmasker} = $program
    }
    if ($program->{description} eq "trf") {
        $repeatmasker->{trf} = $program
    }
    if ($program->{description} eq "rmblast") {
        $repeatmasker->{rmblast} = $program
    }
    if ($program->{description} eq "Blast") {
        $repeatmasker->{Blast} = $program
    }
}
installRM($repeatmasker);


# Install remaining packages
foreach my $program (@$yaml) {
    if ($program->{description} !~ "(repeatmasker|trf|rmblast|Blast)") {
        compileMe($program) unless ($program->{installed});
    }
}

foreach my $program (@$yaml) {
    if ($program->{description} !~ "(repeatmasker|trf|rmblast|Blast)") {
        linkfiles($program) if ($program->{filename});
    }
}

