#!/usr/bin/env perl

use warnings;
use strict;

=head1 AUTOR

Martin Bens (bensmartin@gmail.com)

=head1 DESCRIPTION

Uses SGE to execute BLAST|BLAT|... in parallel.

=head1 SYNOPSIS

perl parallel_blast.pl -db name -query file.fa -out myresult -jobs 2 [-format xml -flavor
ncbi -type blastn -parameter "-num_alignments 1" -debug]

=head1 OPTIONS

=over 8

=item B<-db>

Database

=item B<-query>

Query (in fasta format)

=item B<-out>

Output file

=item B<-type>

Only if using blast:

[optional] blastn|blastp|blastx|tblastx [default: blastn]

=item B<-jobs>

Number of jobs.

=item B<-flavor>

blat|megablast|ncbi|wu|genblast [default: wu]

=item B<-format>

text|xml|tabular [default: tabular]

=item B<-parameter>

Optional parameters passed to blast|blat.

=item B<-bin>

Directory of executable (blastn, blastall, megablast, ...)

=back

=cut

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use File::Basename;
use File::Path qw(make_path remove_tree);
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use Bio::SeqIO;

if ( @ARGV == 0 ) {
    pod2usage( -message => "\n\tNot enough parameter\n", -verbose => 1 );
    die;
}

$SIG{'INT'} = 'CLEANUP';

my %formats = (
    "wu-xml"            => "-mformat 7",
    "wu-tabular"        => "-mformat 2",
    "wu-text"           => "",
    "ncbi-xml"          => "-m 7",
    "ncbi-tabular"      => "-m 8",
    "ncbi-text"         => "",
    "megablast-text"    => "-D 2",
    "megablast-tabular" => "-D 3",
    "blat-tabular"      => "-noHead",
    "blat-text"         => "-out=wublast",
    "custom"            => "",
);

my %opt = (
    dir       => File::Spec->curdir(),
    out       => File::Spec->catfile( File::Spec->curdir(), "output.txt" ),
    flavor    => "wu",
    format    => "tabular",
    type      => "blastn",
    parameter => "",
    debug     => 0,
    sge_param => "-l mem_free=3000M -m n -q smallram -S /usr/local/bin/bash",
    user      => `whoami`,
    engine    => "sge",
    help      => 0,
    jobs      => 10,
    path      => {
        "bin"          => "",
        "bash"         => "/usr/local/bin/bash",
        "qsub"         => "qsub",
        "qstat"        => "qstat",
    }
);


chomp( $opt{user} );

GetOptions(
    \%opt,         "query|q=s",     "db|d=s",   "out|o=s",
    "jobs|j=i",    "parameter|p=s", "type|t=s", "format|fo=s",
    "flavor|fl=s", "debug|?",       'help|?',   "engine=s",
    "sge_param=s", "path=s%"
) or pod2usage(1);

pod2usage( -verbose => 2, -exitval => 2 ) if $opt{help};

unless ( -e $opt{query} ) {
    die "Couldn't find query: $opt{query}\n";
}

unless ( -e $opt{db} ) {
    die "Couldn't find database: $opt{db}\n";
} else {
    if ($opt{flavor} eq "wu") {
        my $command = $opt{path}->{bin} . "xdformat -n $opt{db} &> /dev/null" if ($opt{type} =~ /n$/);
        $command = $opt{path}->{bin} . "xdformat -p $opt{db} &> /dev/null" if ($opt{type} =~ /[px]$/);
        system($command);
    } elsif ($opt{flavor} eq "ncbi") {
        system( $opt{path}->{bin} . "makeblastdb -dbtype nucl -in $opt{db} &> /dev/null") if ($opt{type} =~ /n$/);
        system( $opt{path}->{bin} . "makeblastdb -dbtype prot -in $opt{db} &> /dev/null") if ($opt{type} =~ /[px]$/);
    }
}

my $jobid = 0;             # sge job
my @chars = ( 0 .. 9 );    # alphabet for random identitfier

my $rand_string = "SGE_";
$rand_string .= $chars[ rand @chars ] for 1 .. 5;

checkParameter();

( undef, my $outdir, my $outfile ) = File::Spec->splitpath( $opt{out} );
$outdir = "." if ( $outdir eq "" );
make_path($outdir) unless ( -d $outdir );
$opt{tmp} = File::Spec->catdir( $outdir, "_tmp_$rand_string" );
make_path( $opt{tmp} );

my $status_files_combined = 1;
if ( $opt{engine} eq "parallel" || $opt{engine} eq "sge" ) {

    my @input_files;

    my $filename_query = basename( $opt{query} );
    $filename_query =~ s/(\.[^\.]+$)//;

    # get number of sequences
    my $num_seqs = get_SeqNumber( $opt{query} );

    if ( $num_seqs == 0 ) {
        die "\nNo sequences found!\n\n";
    }

    # create prefix with path and filename
    my $prefix = File::Spec->catfile( $opt{tmp}, $rand_string . "_input" );

    # start splitting
    $opt{jobs} = $num_seqs if ( $num_seqs < $opt{jobs} );

    my $seq_per_file = int( ( $num_seqs / int( $opt{jobs} ) ) + 0.5 );
    my $seq_count    = 0;
    my $file_count   = 0;

    # first file
    my $out_file = $prefix . "." . $file_count;
    push @input_files, $out_file;

    open my $out_seq, ">", $out_file
      or die "Can't open file for writing: $out_file\n";

    my $in = Bio::SeqIO->new( -file => $opt{query} );
    while ( my $seq = $in->next_seq ) {
        $seq_count++;
        if ( $seq_count <= $seq_per_file ) {
            toFasta( $seq->id, $seq->seq, undef, $out_seq );
        } else {
            $file_count++;
            $out_file = "$prefix.$file_count";
            push @input_files, $out_file;

            close $out_seq;
            open $out_seq, ">", $out_file
              or die "Can't open file for writing: $out_file\n";

            toFasta( $seq->id, $seq->seq, undef, $out_seq );
            $seq_count = 1;
        }
    }
    close $out_seq;
    $opt{jobs} = $file_count + 1;

    my @output_files;
    if ( $opt{engine} eq "parallel" ) {
        #my $pm = Parallel::ForkManager->new( $opt{jobs} );
	my $qsub_filename = File::Spec->catfile( $opt{tmp}, $rand_string . "_commands.sh" );
	open my $fh, ">", $qsub_filename;
	for (@input_files) {
		my ( $command, $output_file ) = getCommand($_);
		push @output_files, $output_file;
		print $fh $command."\n";
	}
	close $fh;
	system("cat $qsub_filename | parallel");

    } elsif ( $opt{engine} eq "sge" ) {
        my ( $blast_command, $output ) = getCommand( $prefix . ".\$i" );

        my $sge_out = $opt{out} . ".sgeout.log";
        my $sge_err = $opt{out} . ".sgeerr.log";

        my @options = (
            "#!" . $opt{path}->{bash},
            "#\$ -t 1-$opt{jobs}",
            "#\$ -o $sge_out",
            "#\$ -e $sge_err"
        );
        my @commands =
          ( "i=\$(expr \$SGE_TASK_ID - 1)", $blast_command, "exit 0" );

        my $qsub_filename =
          File::Spec->catfile( $opt{tmp}, $rand_string . "_commands.sh" );

        open( FH, ">$qsub_filename" );
        print FH join( "\n", ( @options, @commands ) );
        close(FH);

        # send array job to sge
        my $command = $opt{path}->{qsub} . " $opt{sge_param} $qsub_filename";
        my $output_line = `$command`;
        unless ($output_line) {
            die "Couldn't send job to SGE\n";
        }

        if ( $output_line
            =~ /Your job-array (\d+)\.(.+?) \((.+?)\) has been submitted/ )
        {
            $jobid = $1;
        } else {
            removeFiles(1);
        }
        my $qstat_command = $opt{path}->{qstat} . " -u \$(whoami)";
        my $running       = 1;
        while ($running) {
            sleep(30);
            $running = 0;
            my $output = `$qstat_command`;
            $running = 1 if ( $output =~ /$jobid/ );
        }

        if ( $opt{flavor} eq "repeatmasker" ) {
            @output_files =
              map {
                File::Spec->catfile( $opt{tmp},
                    $rand_string . "_input.$_.masked" )
              } 0 .. ( $opt{jobs} - 1 );
        } else {
            @output_files = map {
                File::Spec->catfile( $opt{tmp}, $rand_string . "_output.$_" )
            } 0 .. ( $opt{jobs} - 1 );
        }
    }

	sleep(1);
    my $status = combineResults( \@output_files );
    removeFiles($status);
} else {
    my ( $command, $output ) = getCommand( $opt{query}, $opt{out} );

    my $status = system($command);
    removeFiles($status);
}

if (-z $opt{out}) {
    die "Failed";
}

sub getCommand {
    my $input  = shift;
    my $output = shift;

    ( $output = $input ) =~ s/input/output/ unless ($output);

    my $command;
    if ( $opt{flavor} eq 'wu' ) {
        $command = $opt{path}->{bin}
          . "$opt{type} $opt{db} $input -o=$output $formats{$opt{format}} $opt{parameter} > /dev/null";
    } elsif ( $opt{flavor} eq 'megablast' ) {
        $command = $opt{path}->{bin}
          . "megablast -d $opt{db} -i $input -o $output $formats{$opt{format}} $opt{parameter}";
    } elsif ( $opt{flavor} eq 'blat' ) {
        $command = $opt{path}->{bin}
          . "blat $opt{db} $input $output $formats{$opt{format}} $opt{parameter} &> /dev/null";
    } elsif ( $opt{flavor} eq 'ncbi' ) {
        $command = $opt{path}->{bin}
          . "blastall -p $opt{type} -d $opt{db} -i $input -o $output $formats{$opt{format}} $opt{parameter}";
    } elsif ( $opt{flavor} eq 'genblasta' ) {
        $command = $opt{path}->{bin}
          . "genblasta -P blast -pg $opt{type} -q $input -t $opt{db} -o $output $opt{parameter} &> /dev/null";
    } elsif ( $opt{flavor} eq 'repeatmasker' ) {
        $command =
            "cd $opt{tmp} && "
          . $opt{path}->{bin}
          . "RepeatMasker $input $opt{parameter}";
        $output = $input . ".masked";
    }
    return ( $command, $output );
}

sub removeFiles {
    my $status = shift || 0;
    return if ( $status != 0 );

    unless ( $opt{debug} ) {
        File::Path::rmtree( $opt{tmp} ) or die $!;
    }
}

sub get_SeqNumber {
    my $file  = shift;
    my $i_seq = 0;
    open FILE, "<", $opt{query} or die "Can't open $opt{query}\n";
    while (<FILE>) { $i_seq++ if (/^>/); }
    close FILE;
    return $i_seq;
}

sub checkParameter {

    unless ( defined $opt{format} ) {
        $opt{format} = "custom";
    } elsif ( $opt{format} eq "table" ) {
        $opt{format} = "tabular";
    }

    $opt{format} = $opt{flavor} . "-" . $opt{format};
    $opt{db}     = File::Spec->rel2abs( $opt{db} );
    $opt{query}  = File::Spec->rel2abs( $opt{query} );
    $opt{out}    = File::Spec->rel2abs( $opt{out} );

}

sub combineResults {
    my @files = @{ shift() };
    my $status;
    my $command;

    # adding repeatmasker (no output if no repeats)
    if ( $opt{flavor} eq "repeatmasker" ) {
        for (@files) {
            unless ( -e $_ ) {
                $_ =~ s/\.masked//;
            }
        }
        $command = "cat " . join( " ", @files ) . " > $opt{out}";

    } elsif ( $opt{flavor} eq "genblasta" ) {
        $command = "cat " . join( " ", @files ) . " > $opt{out};";
        my @blast_files = map { s/SGE_(\d+)_output\./SGE_$1_input./; $_ . ".blast" } @files;
        $command .= "cat " . join( " ", @blast_files ) . " > $opt{out}.blast;";
    } else {
        if ( $opt{format} =~ /tabular/ ) {
            $command = "cat " . join( " ", @files ) . " > $opt{out}";
        } elsif ( $opt{format} =~ /text/ ) {
            $command = "cat " . join( " ", @files ) . " > $opt{out}";
        }
    }

    $status = system($command);
    return $status;
}

sub CLEANUP {

    #$logger->info("Caught Interrupt (^C), Aborting");

    print STDERR "Aborting.\n";
    system("qdel $jobid") if ($jobid);

    print STDERR "Removing files..\n";
    File::Path::rmtree( $opt{tmp} );

    exit(1);
}

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

