#!/usr/bin/env perl

use warnings;
use strict;

=pod
=head1 INFO
Martin Bens, bensmartin@gmail.com
2014-07-30

=head1 DESCRIPTION

Infers gene boundaries of single or multi-gene contigs.

Example output:

    #strand  prev_cds_end  centerL  centerR  prom_start  prom_end  clip  clipscore  start  end   transcript  gene
    1        1             1        2968     -1          -1        3296  1.250      1      3296  ACC_GENA   geneA
    -1       2529          2968     4471     4511        4550      3273  1.019      3273   4471  ACC_GENB   geneB
    1        4084          4471     5494     4813        4852      5481  3.642      4852   5494  ACC_GENC   geneC

Last 4 columns are important. Includes start and end of gene. Other columns
describe possible cut positions. CenterL and centerR represent center of
intergenic region to neighbooring coding sequences. Prom_start/end refers to
promotor information. Clip and Clipscore result of 3' UTR clipping (based on
read data, ortholog and PAS).

=head1 OPTIONS

=over 8

=item B<-utr-length>

Specify length of reference UTR (length applies to UTR after trimming of
'N' and 'A').

=item B<-prefix>

Prefix of output files (default: mRNA)

=item B<-output>

Output directory (default: out).

=item B<-ortholog>

Fasta file with full length sequence of orthologs.

=item B<-contig>

Fasta file with sequence to clip.

=item B<-readfile>

All reads in fasta format (poly(A) read step). Separate multiple files with
",".

=item B<-input>

Tab delimited table with annotated CDS regions (and Promoters if known).
Otherwise, only clipping 3' end based on best score.

    #feature	symbol	accession	start	end	  strand
    CDS	        geneA	ACC_GENA	1261	2529   1
    CDS	        geneB	ACC_GENB	3406	4084  -1
    PRO	        geneB	ACC_GENB	4511	4550  -1
    CDS	        geneC	ACC_GENC	4856	5344   1
    PRO	        geneC	ACC_GENC	4813	4852   1


=back

=cut

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use File::Temp qw(tempfile tempdir);
use File::Path qw(make_path remove_tree);
use File::Spec::Functions qw(catfile catdir splitpath);

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SeqIO;
use Pod::Usage;
use Data::Dumper;

use Getopt::Long;
pod2usage(-message => "\n\tNo arguments. See -help.\n", -verbose => 1) if (@ARGV == 0);

use constant NEEDLE_OPEN => 10;
use constant NEEDLE_EXT  => 10;

use constant BOWTIE_GBAR => 1000;
use constant BOWTIE_DPAD => 0;

use constant UTR_LENGTH => 50;

use constant WEIGHT_PAS  => 1.18;
use constant WEIGHT_READ => 0.335;
use constant WEIGHT_DROP => 0.285;
use constant WEIGHT_UTR  => 0.99;
use constant MIN_SCORE   => 1.0;

use constant MAX_INT => 9**9**9;

my %opt = (
    PATH_POLYA    => "$FindBin::Bin/polyamisc/",
    'utr-length' => 50,
    'prefix'     => "mRNA",
    'output'     => 'out',
    'graphics'   => 0,
    'debug'      => 0,
);
GetOptions(
    \%opt,        'help|?', 'output_dir=s', 'prefix=s', 'contig=s',
    'ortholog=s', 'readfile=s', 'utr-length=i', 'graphics', 'debug',
    'input=s',    'output=s'
);

GetOptions() or pod2usage(1);

pod2usage(-verbose => 2) if $opt{help};

# mandatory
unless ($opt{readfile} && $opt{ortholog} && $opt{contig} && $opt{input}) {
    pod2usage(-message => "\n\tNot enough arguments. See -help.\n");
}

my %cds_range;
my %promoters;

# read cds ranges
if (exists $opt{input}) {
    unless (-e $opt{input}) {
        die "\nInput file not found: $opt{input}\n\n";
    }

    open my $fh, "<", $opt{input} or die $!;
    while (defined(my $line = <$fh>)) {
        next if ($line =~ /^#/);
        chomp $line;
        my @e = split "\t", $line;
        unless (@e >= 6) {
            print STDERR "Could not read line: $.\n";
            next;
        }

        # ignore results not supported by blast
        next if ($e[2] eq "UNKNOWN");

        $cds_range{$e[2]} = {start => $e[3], end => $e[4], strand => $e[5], sym => $e[1], seq => undef, transcript => $e[2]} if ($e[0] eq "CDS");
        $promoters{$e[2]} = {start => $e[3], end => $e[4], strand => $e[5], sym => $e[1], seq => undef, transcript => $e[2]} if ($e[0] eq "PRO");
    }
    close $fh;
}

# add ortholog sequence
my $io = Bio::SeqIO->new(-file => $opt{ortholog}, -format => "fasta");
if (%cds_range) {
    while (my $seq = $io->next_seq) {
        if ($cds_range{$seq->id}) {
            $cds_range{$seq->id}->{seq} = $seq;
        }
    }
} else {
    my $seq = $io->next_seq;
    $cds_range{UNKNOWN} = {start => 1, end => MAX_INT, strand => 1, sym => "UNKNOWN", seq => $seq, transcript => "UNKNOWN"};
}

# read contig
$io = Bio::SeqIO->new(-file => $opt{contig}, -format => "fasta");
my $contig_seq = $io->next_seq if ($io);

unless (defined $contig_seq) {
    die "\n\tContig sequence not found\n";
}

# create output dir
unless (defined $opt{output_dir}) {
    $opt{output_dir} = tempdir(CLEANUP => 1);
} elsif (not -d $opt{output_dir}) {
    make_path($opt{output_dir});
}

# check directions of coding regions
my $both_directions = 0;
for my $value (values %cds_range) {
    if ($value->{strand} eq "-1") {
        $both_directions = 1;
        # TODO: do not realign reads!
        last;
    }
}

# trim contig
($contig_seq) = maskEnds($contig_seq);

my $polyA = 0;
if ($contig_seq->seq =~ /(N+)$/) {
    $polyA = length($1);
}
my $length = $contig_seq->length - $polyA;

my ($contig_fh, $contig_file) = tempfile(UNLINK =>1);

my $out = Bio::SeqIO->new(-file => ">".$contig_file, -format => "fasta");
$out->write_seq($contig_seq);

# compute all features
my %files = (
    contig => $contig_file,
    pas    => catfile($opt{output_dir}, "$opt{prefix}_PAS.txt"),
    cov    => catfile($opt{output_dir}, "$opt{prefix}_coverage.txt"),
    drop   => catfile($opt{output_dir}, "$opt{prefix}_dropcoverage.txt"),
    polyA  => catfile($opt{output_dir}, "$opt{prefix}_polyA.txt"),
    bowtie => catfile($opt{output_dir}, "$opt{prefix}_bowtie_alignment"),
);
run_orth_ind(\%files, $contig_seq);

# compute all feature for opposite strand
my %files_revcom;
my ($contig_seq_rev);
if ($both_directions) {
    ($contig_seq_rev) = ($contig_seq->revcom);

    %files_revcom = (
        contig => catfile($opt{output_dir}, "$opt{prefix}_contig_rev.fa"),
        pas    => catfile($opt{output_dir}, "$opt{prefix}_PAS_rev.txt"),
        cov    => catfile($opt{output_dir}, "$opt{prefix}_coverage_rev.txt"),
        drop   => catfile($opt{output_dir}, "$opt{prefix}_dropcoverage_rev.txt"),
        polyA  => catfile($opt{output_dir}, "$opt{prefix}_polyA_rev.txt"),
        bowtie => catfile($opt{output_dir}, "$opt{prefix}_bowtie_alignment_rev"),
    );

    my $io = Bio::SeqIO->new(-file => ">".$files_revcom{contig}, -format => "fasta");
    $io->write_seq($contig_seq_rev);

    run_orth_ind(\%files_revcom,$contig_seq_rev);
}

# align UTR from corresponding ortholog and score features
for my $key (keys %cds_range) {

    my $pair = $cds_range{$key};

    # output files
    my $utr_aln_file  = catfile($opt{output_dir}, "$opt{prefix}_" . $pair->{sym} . "_UTR.aln");
    my $utr_id_file   = catfile($opt{output_dir}, "$opt{prefix}_" . $pair->{sym} . "_UTR.txt");

    my %parameter = (
        ortholog     => $pair->{seq},
        output_aln   => $utr_aln_file,
        output_txt   => $utr_id_file,
    );
    if ($pair->{strand} eq "-1") {
        $parameter{contig} = $contig_seq_rev;
    } else {
        $parameter{contig} = $contig_seq;
    }

    UTRaln(\%parameter) == 0 || die "\nError: UTR alignment failed\n\n";

    my $sumscore = catfile($opt{output_dir}, "$opt{prefix}_" . $pair->{sym} . "_sumscore.txt");
    my $sumscore_graph = catfile($opt{output_dir}, "$opt{prefix}_" . $pair->{sym} . "_sumscore.png") if ($opt{graphics});

    %parameter = (
        UTR          => $utr_id_file,
        output       => $sumscore,
        output_graph => $sumscore_graph,
        length => $contig_seq->length,
    );
    if ($pair->{strand} eq "-1") {
        $parameter{polyA}        = $files_revcom{polyA};
        $parameter{PAS}          = $files_revcom{pas};
        $parameter{COV}          = $files_revcom{drop};
        $parameter{revcom} = 1;
    } else {
        $parameter{polyA}        = $files{polyA};
        $parameter{PAS}          = $files{pas};
        $parameter{COV}          = $files{drop};
        $parameter{revcom} = 0;
    }

    mRNA_featureSelection(\%parameter) == 0 || die "\nError: Couldn't compute sumscore.\n\n";

    $pair->{scorefile} = $sumscore;
}

my @utr_regions;

my @sorted_keys = sort { $cds_range{$a}->{start} <=> $cds_range{$b}->{start} } keys %cds_range;

print join "\t", "#strand", "prev_cds_end", "centerL", "centerR", "prom_start", "prom_end", "clip", "clipscore", "start", "end", "transcript", "gene";
print "\n";

if (@sorted_keys == 1) {
    my ($key_A) = @sorted_keys;
    my $A = $cds_range{$key_A};
    if ($A->{strand} == 1) {

        # use promotor information if available
        my $start = 1;
        my $prom = $promoters{$A->{transcript}};
        my ($prom_start, $prom_end ) = (-1, -1);
        if ($prom) {
            ($prom_start, $prom_end) = ($prom->{start}, $prom->{end});
            $start = $prom_end;
        }

        my $end = $contig_seq->length;
        my ($pos, $score) = bestClip($A->{scorefile}, $A->{end});
        if ($pos > 0 && $pos < ($length-30)) {
            $end = $pos
        }

        print "1\tNA\tNA\tNA\t$prom_start\t$prom_end\t$pos\t$score\t$start\t$end\t$key_A\t".$A->{sym}."\n";
    } else {
        # annotated CDS should be in transcript region.
        die "Something went totally wrong. Transcript wasn't orientated.\n";
    }
} else {
    my @centers = (1);
    for (my $i = 0; $i < @sorted_keys-1; $i++) {
        my $key_A = $sorted_keys[$i];
        my $key_B = $sorted_keys[$i+1];
        push @centers, ($cds_range{$key_A}->{end} +(int ( ( ($cds_range{$key_B}->{start} - $cds_range{$key_A}->{end} + 1) / 2 ) + 0.5) ));
    }
    push @centers, $contig_seq->length;

    my $i = 0;
    my $previous_cds = 1;
    for (; $i < @sorted_keys; $i++) {
        my $key_A = $sorted_keys[$i];
        my $key_B = $sorted_keys[$i+1];

        my $A = $cds_range{$key_A};
        my $B = $cds_range{$key_B} if ($key_B);

        if ($A->{strand} eq -1) {

            # start
            my ($pos, $score) = bestClip($A->{scorefile}, $previous_cds, $A->{start});
            my $start = $pos if ($pos > 30);
            if ($pos > 30) {
                $start = $pos;
            } else {
                $start = $centers[$i];
            }

            # end
            my $prom = $promoters{$A->{transcript}};
            my ($prom_pos, $prom_start, $prom_end) = (MAX_INT, -1, -1);
            if ($prom) {
                ($prom_start, $prom_end) = ($prom->{start}, $prom->{end});
                $prom_pos = $prom_start;
            }
            my $end = ($centers[$i+1] < $prom_pos) ? $centers[$i+1]  : $prom_pos;

            print "-1\t$previous_cds\t$centers[$i]\t$centers[$i+1]\t$prom_start\t$prom_end\t$pos\t$score\t$start\t$end\t$key_A\t".$A->{sym}."\n";

        } else {

            my $next_cds = $contig_seq->length;
            if ($B) {
                $next_cds = $B->{start};
            }

            # extract first gene based on which feature produces smaller
            # transcripts (promotor end, center between CDSs)
            my ($prom_pos, $prom_start, $prom_end) = (-1, -1, -1);

            # start
            my $prom = $promoters{$A->{transcript}};
            if ($prom) {
                ($prom_start, $prom_end) = ($prom->{start}, $prom->{end});
                $prom_pos = $prom->{end};
            }
            my $start = ($centers[$i] > $prom_pos) ? $centers[$i] : $prom_pos;

            # end
            my ($pos, $score) = bestClip($A->{scorefile}, $A->{end}, $next_cds);
            my $end = $centers[$i+1];
            $end = $pos if ($pos > 0 && $pos < $length-30);

            print "1\t$previous_cds\t$centers[$i]\t$centers[$i+1]\t$prom_start\t$prom_end\t$pos\t$score\t$start\t$end\t$key_A\t".$A->{sym}."\n";
        }
        $previous_cds = $A->{end};
    }
}

sub bestClip {
    my ($file, $start, $end) = @_;

    $start = 0 unless ($start);
    $end = MAX_INT unless ($end);

    my $highscore = -1;
    my $maxpos    = -1;
    open my $fh, "<", $file or die $!;
    while (<$fh>) {
        next if (/^#/);
        chomp;

        my @e = split;

        if ($e[2] > MIN_SCORE && $e[1] > $start && $e[1] < $end && $e[2] > $highscore ) {
            $highscore = $e[2];
            $maxpos    = $e[1];
        }
    }
    close $fh;

    return ($maxpos, $highscore);
}

sub run_orth_ind {
    my ($files, $contig_seq) = @_;

    polyASignal({contig => $files->{contig}, output => $files->{pas}}) == 0
    || die "\nError: Identification of PAS failed\n\n";

    alignBowtie({contig => $files->{contig}, output => $files->{bowtie}}) == 0
    || die "\nError: read alignment failed\n\n";

    covDrop({input => "$files->{bowtie}.bam", cov => $files->{cov}, drop => $files->{drop}}) == 0
    || die "\nError: coverage\n\n";

    polyAreads({input => $files->{bowtie}, output => $files->{polyA}, cov => $files->{cov}}) == 0
    || die "\nError: filtering of poly(A)reads failed\n\n";

}

sub mRNA_featureSelection {
    my ($args) = @_;

    my @command = (
        "$opt{PATH_POLYA}polyaclip_sumscore.pl",
        $args->{PAS}, $args->{polyA}, $args->{COV}, $args->{UTR},
        WEIGHT_PAS, WEIGHT_READ, WEIGHT_DROP, WEIGHT_UTR, ">", $args->{output}
    );
    my $failed = system(join " ", @command);

    if ($args->{revcom}) {
        open my $fh, "<", $args->{output} or die $!;
        my @content = <$fh>;
        close $fh;

        my @header = @content[0..2];
        @content = reverse @content[3..$#content];

        open my $out, ">", $args->{output} or die $!;

        for (@header) {
            print $out $_;
        }
        my $counter = 0;
        for (@content) {
            $counter++;
            s/^(.+?)\t(.+?)\t/$1\t$counter\t/g;
            print $out $_;
        }
        close $out;
    }

    return $failed;
}

sub covDrop {
    my ($args) = @_;

    my $failed = 0;
    unless (-z $args->{input} ) {

        my $command = "bamtools coverage -in $args->{input} > $args->{cov}";

        $failed  = system($command);
        unless ($failed) {
            $command =
            "cat $args->{cov} | $opt{PATH_POLYA}polyacov_dropcov_simple.pl - > $args->{drop} ";
            $failed = system($command);
        }
    } else {
        # subsequent scripts handle empty files
        system("touch $args->{cov}");
        system("touch $args->{drop}");
        system("touch $args->{input}");
    }

    return $failed;
}

# Description
#
#   Returns end of 3p UTR
#
# Input
#
#   {seq} object
#   {length} of UTR
sub getUTRend {
    my ($args) = @_;

    my $utr = $args->{seq}->seq;
    $utr =~ s/[AN]+$//gi;
    $utr = substr $utr, -$args->{length};
    return Bio::Seq->new(
        -id  => $args->{seq}->id . " _ " . $args->{length} . " UTR ",
        -seq => $utr
    );
}

sub maskEnds {
    my ($seq) = @_;

    my $seq_string = $seq->seq;
    $seq_string =~ s/([AN]+)$/"N" x length($1)/egi;
    $seq_string =~ s/^([TN]+)/"N" x length($1)/egi;

    my $clone = $seq->clone();
    $clone->seq($seq_string);
    return $clone;
}

#
# Input
#
#   {utr} seq
#   {contig} seq
#   {output} file
#
# Returns
#
#   status
#
sub alignNeedle {
    my ($args) = @_;

    for (qw/output contig/) {
        unless (defined $args->{$_}) {
            print STDERR " Not defined : $_ \n ";
            return;
        }
    }

    return 0 if (-e $args->{output});

    my ($fh_a, $file_a) = tempfile(UNLINK => 1);
    my ($fh_b, $file_b) = tempfile(UNLINK => 1);

    print $fh_a ">" . $args->{utr}->id . "\n ";
    print $fh_a $args->{utr}->seq . "\n";
    close $fh_a;

    (my $seq_string = $args->{contig}->seq) =~ s/N+$//gi;
    print $fh_b ">" . $args->{contig}->id . "\n";
    print $fh_b $seq_string. "\n";
    close $fh_b;

    my $command =
        "needle -gapopen "
      . NEEDLE_OPEN
      . " -gapextend "
      . NEEDLE_EXT
      . " -aformat fasta -asequence $file_a -bsequence $file_b -outfile $args->{output} &> /dev/null";
    my $failed = system($command);

    close $fh_a;
    close $fh_b;

    return $failed;
}

sub polyASignal {
    my ($args) = @_;

    unless (defined $args->{contig} && defined $args->{output}) {
        print STDERR "Poly(A)-Signal: not enough arguments!";
        return;
    }

    # polyA scoring
    my @commands = (
        "cat $args->{contig} |",
        #"$opt{PATH_POLYA}/fasta_trimpolya.pl - |",
        "$opt{PATH_POLYA}polyaspwm_predict.pl - |",
        "$opt{PATH_POLYA}polyaspwm_score.pl - >",
        "$args->{output}"
    );
    my $failed = system(join " ", @commands);

    return $failed;
}

# Input
#   {contig} Bio::Seq
#   {readfile} string
#   {output} string
#
# Returns
#   1 if failes.
#
sub alignBowtie {
    my ($args) = @_;

    for (qw/output contig/) {
        unless (defined $args->{$_}) {
            print STDERR "Not defined: $_ \n ";
            return;
        }
    }


    my $tempdir = tempdir (CLEANUP => 1);
    my $bowtie_index = catfile($tempdir, "bowtie_index");
    my $sam = catfile($tempdir, "output");

    my $failed = system("bowtie2-build -q $args->{contig} $bowtie_index 2> /dev/null");

    if ($failed) {
        print STDERR "Failed to build bowtie index. Skipping.\n";
        return 1;
    }

    my $reads = catfile($tempdir, "read_input");

    my @files = split ",", $opt{readfile};
    for (@files) {
        $failed =  system("perl -pe 's/^>@/>/' $_ >> $reads");
        if ($failed) {
            print STDERR "Read preprocessing went wrong\n\n";
            return 1;
        }
    }
    my $command = "bowtie2 --gbar "
    . BOWTIE_GBAR
    . " --dpad "
    . BOWTIE_DPAD
    . " --no-unal --n-ceil C,0 --quiet -f --very-sensitive-local -p 1 -x $bowtie_index -U $reads -S $sam";

    $failed = system($command) unless ($failed);

    if ($failed) {
        print STDERR "Failed to align reads\n";
        return 1;
    }

    my $aligned = 0 ;
    open my $fh, "<", $sam or die "Can't open file for reading: $sam\n";
    while(<$fh>) {
        next if (/^@/);
        $aligned++;
        last;
    }
    close $fh;

    if ($aligned == 0) {
        print STDERR "No reads could be aligned\n";
        system("rm $sam");
    } else {

		$failed = system("samtools view -bS $sam | samtools sort -o $args->{output}.bam -");

        if ($failed) {
            print STDERR "Failed to convert SAM to BAM.\n";
            return $failed;
        }

    }

    return $failed;
}

sub UTRaln {
    my ($args) = @_;

    for (qw/contig ortholog output_aln output_txt/) {
        unless (defined $args->{$_}) {
            print STDERR " Not defined : $_ \n ";
            return;
        }
    }

    my $UTR = getUTRend(
        {
            seq    => $args->{ortholog},
            length => $opt{'utr-length'}
        }
    );
    my $failed = alignNeedle(
        {
            utr    => $UTR,
            contig => $args->{contig},
            output => $args->{output_aln}
        }
    );
    return 1 if ($failed);

    my @command = (
        "$^X", "$FindBin::Bin/polyamisc/UTR.pl", "-end-gap 0.5",
        "-i $args->{output_aln} |",
        "$^X $opt{PATH_POLYA}polyahomol_score.pl - >",
        $args->{output_txt}
    );
    $failed = system(join " ", @command);

    return $failed;
}

sub polyAreads {
    my ($args) = @_;

    for (qw/input output/) {
        unless (defined $args->{$_}) {
            print STDERR " Not defined : $_ \n ";
            return;
        }
    }

    my @commands = (
        "$^X $FindBin::Bin/polyamisc/bam_polyA.pl",
        "-i $args->{input}.bam |",
        "$^X $opt{PATH_POLYA}/polyaread_score.pl -d $args->{cov} - >",
        $args->{output}
    );
    return system(join " ", @commands);
}

1;

__END__
