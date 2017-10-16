#!/usr/bin/env perl

# --------------------------------------------------------------------
# Install required tools
#
# Martin Bens | bensmartin@gmail.com
# 2017-08-27
# --------------------------------------------------------------------

use strict;
use warnings;

use v5.10;

use YAML;
use Data::Dumper;
use File::Spec;
use File::Copy;

my $download_dir = "download";
my $unpack_dir = "sources";
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

# Determine software to install
sub checkProgram {
    my ($data) = @_;

    say "- $data->{description}";

    my @programs = split " ", $data->{executable};
    my $installed = 0;
    if (@programs) {
        for my $e (@programs) {
            my $command = "which $e";
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

# Download required software
sub download {
    my ($h) = @_;

    my $filename = "$download_dir/$h->{description}";
    $filename .= "." . $h->{extension} if (exists $h->{extension});

    my $wget = "wget -O $filename -o $filename.log";
    $wget .= " " . $h->{wget} if (exists $h->{wget});

    $wget .= " " . $h->{remote};

    if (-e $filename) {
        say "File already exists: $filename. Skipping.";
        $h->{download} = $filename;
        return 1;
    }

    if (defined $h->{licence}) {
        say "Licence";
        say $h->{licence};
        say "";
    }

    say "Running:\n\t$wget";
    my $error = system($wget);
    if (!$error) {
        $h->{download} = $filename;
    } else {
        say "Download failed: $h->{description} from $h->{remote}";
        return 0;
    }

    return 1;
}

# Unpack everything accordingly
# tar.gz, tar.bz2, zip
sub unpackMe {
    my ($h) = @_;

    # nothing to do
    if (not exists $h->{dir}) {
        my $s = File::Spec->rel2abs("$download_dir/$h->{links}");
        copy($s, "$unpack_dir/$h->{links}");
        return 1;
    }

    # gz bz2 zip
    my $command = "";
    if ($h->{extension} eq "tar.gz") {
        $command = "tar zxf";
    } elsif ($h->{extension} eq "tar.bz2") {
        $command = "tar jxf";
    } elsif ($h->{extension} eq "zip") {
        $command = "unzip -qq -d $unpack_dir";
    }
    $command .= " " . $h->{download};

    # custom directory
    if ($command =~ /^tar/) {
        $command .= " -C $unpack_dir";
    }
    my $mydir = "$unpack_dir/$h->{dir}";

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
        say "Unpacking failed: $h->{description} from $h->{download}";
    }
}

# start default compilation in /sources ..
sub compileMe {
    my ($h) = @_;

    if ((not exists $h->{compile}) || (not exists $h->{dir})) {
        say "Nothing to compile: $h->{description}";
        return 1;
    }

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

# create link in /bin
sub linkfiles {
    my ($data, $type) = @_;

    my @sources = split " ", $data->{links};
    my @target  = getFinal($data);

    for (my $i = 0; $i < @sources; $i++) {
        if (-e $target[$i]) {
		print "$target[$i] exists.\n";
		next;
	}

        my $source = File::Spec->catfile($unpack_dir, $sources[$i]);
        $source = File::Spec->catfile($unpack_dir, $data->{dir}, $sources[$i])
            if (exists $data->{dir});

        $source = File::Spec->rel2abs($source);
        chmod(0755, $source);
        #my $target = File::Spec->catfile($execute_dir, $target[$i]);
        symlink($source, $target[$i]);

	print "$source --> $target[$i]\n";

        say "Executable added: $target[$i]";
    }
}

# Repeatmasker is special..
sub installRM {
    my ($data) = @_;

    # ncbi-blast
    if (not $data->{Blast}->{installed}) {
        compileMe($data->{Blast});
    }
    linkfiles($data->{Blast});

    my $repmask = File::Spec->catfile($execute_dir, $data->{repeatmasker}->{executable});
    if (-e $repmask) {
	    say "RepeatMasker installed. Skipping.\n";
	    return 1;
    }

    linkfiles($data->{rmblast});
    linkfiles($data->{trf});

    # required files..
    my $blast_bin = File::Spec->rel2abs("$execute_dir");
    my $trf_bin = File::Spec->rel2abs("$execute_dir/trf");

    # install repeatmasker
    my $command = "cd $unpack_dir/$data->{repeatmasker}->{dir} && printf \"\\n$^X\\n\\n$trf_bin\\n2\\n$blast_bin\\nY\\n5\\n\" | perl configure";
    my $error = system($command);
    linkfiles($data->{repeatmasker});
}

my $yaml = YAML::LoadFile('programs.yaml');

mkdir $download_dir;
mkdir $unpack_dir;
mkdir $execute_dir;

# Check installation
foreach my $program (@$yaml) {
    checkProgram($program);
}

# Download files
foreach my $program (@$yaml) {
    download($program) unless ($program->{installed});
}

# Unpack files
foreach my $program (@$yaml) {
    unpackMe($program) if ($program->{download});
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
    next if ($program->{description} =~ /^(repeatmasker|trf|rmblast|Blast)$/);

    my ($name) = getFinal($program);
    if (-e $name || $program->{installed}) {
        say "Skipping: $name";
        next;
    }
    say "Compiling: $name";
    compileMe($program);
    if ($program->{download}) {
	    say "Linking: $name";
	    linkfiles($program)
    }
}

print "\nPlease install \"GENSCAN\" manually: http://genes.mit.edu/license.html\n\n";

print << 'END';
    FRAMA_DIR=$(readlink -f .)
    wget http://genes.mit.edu/XXX
    mkdir genscan && tar xvf genscanlinux.tar -C genscan
    mv genscan $FRAMA_DIR/sources/.
    ln -f -s $(readlink -f $FRAMA_DIR/sources/genscan/genscan) $FRAMA_DIR/bin/genscan
    ln -f -s $(readlink -f $FRAMA_DIR/sources/genscan/) $FRAMA_DIR/bin/genscan.dir
END

1;

# THE END
