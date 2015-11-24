#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens
bensmartin@gmail.com
12/07/13 18:31:17

=head1 DESCRIPTION

Wrapper, that calls all steps necessary to scaffold, trim and predict CDS based
on assigned orthologs. 

=head1 OPTIONS

=over 8

=item B<-reindex>

Reindex fasta and genbank files.

=item B<-taxon-assembly>

Taxonomy ID of species to assemble.

=item B<-taxon-ortholog>

Taxonomy ID of reference species.

=item B<-species-order>

We keep a note in GenBank output about the ortholog, used to annotated CDS.
Please specify prefered order of species, multiple equal regions originating
from different species have been found. Number (0-based!) refers column in
-ortholog-table. If not specied, column order of -ortholog-table will be used.

=item B<-ortholog-cds> 

Fasta file containing all coding sequences for accessions used in
-ortholog-table.

=item B<-trinity> Trinity.fa

Trinity output.

=item B<-reference> reference.gbk

Reference Transcriptome

=item B<-assignment> annotation.csv

Annotation: Trinity_annotation.csv,

=item B<-ortholog> 

Only write genbank entries for specific transcripts. Comma separated list of
accession numbers.

=item B<-ortholog-file> 

=item B<-contig> 

Only write genbank entries for specific contigs. Comma separated list of contig
ids.

=item B<-contig-file> 

=item B<-output-file> 

Output GenBank file.

=item B<-new> 

Specify if you want to rebuild a transcript (computes all steps again). 

=item B<-cpus> 

Number of CPUs

=item B<-scaffolding> 

Blast result with fragments for scaffolding.

=item B<-blast> 

Blast result to annotated CDS predicted by GENSCAN.

=item B<-genscan-matrix> 

Path to matrix to use for GENSCAN.

=back

=cut

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use GenbankHelper;

use Pod::Usage;
use Benchmark;
use Storable;
use Parallel::ForkManager;
use POSIX qw/strftime/;
use Data::Dumper;
use Getopt::Long;

use File::Path qw(make_path remove_tree);
use File::Spec::Functions qw(catfile catdir splitpath);

use Bio::AlignIO;
use Bio::DB::Fasta;
use Bio::DB::Taxonomy;
use Bio::Location::Simple;
use Bio::Location::Split;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::Tools::Genscan;
use Bio::Taxon;




use constant TDATE => uc( strftime "%d-%b-%Y", localtime );

# Blast
use constant QUERY_ID     => 0;
use constant TARGET_ID    => 1;
use constant STRAND_BLAST => 7;
use constant TARGET_START => 10;
use constant TARGET_END   => 11;
use constant IDENTITY     => 16;
use constant QUERY_COV    => 15;
use constant SYMBOL       => 17;

# Trinity_annotation.csv
use constant TRIN_ID  => 0;
use constant ORTH_ID  => 1;
use constant STRAND   => 2;
use constant HIT_TYPE => 4;
use constant STARTS   => 5;
use constant ENDS     => 6;

# mRNA clipping
use constant UTR_LENGTH => 50;

use constant MAX_INT => 9**9**9;

use constant ORDER_FEATURES => qw/source mRNA CDS promoter misc_feature/;

my %DIVISION = (
    'Primates'      => "PRI",
    'Rodents'       => "ROD",
    'Mammals'       => "MAM",
    'Vertebrates'   => "VRT",
    'Invertebrates' => "INV",
    'Viruses'       => "VRL",
    'Plants'        => "PLN",
    'Phages'        => "PHG",
    'Unassigned'    => "UNA",
    'Bacteria'      => "BCT",
    'UNKNOWN'       => "UNK"
);

if ( @ARGV == 0 ) {
    pod2usage( -message => "\n\tNo arguments. See -help.\n" );
    exit;
}

my $debug = 0;
my %opt   = (
    reindex  => 0,
    debug    => 0,
    reindex  => 0,
    new      => 0,
    cpus     => 0,
    clipping => 1,
);

GetOptions(
    \%opt,              'trinity=s',
    'reference=s',      'assignment=s',
    'output-dir=s',     'output-file=s',
    'reindex',          'scaffolding=s',
    'contig-file=s',    'ortholog-file=s',
    'ortholog=s',       'contig=s',
    'read-index=s',     'debug',
    'reindex',          'new',
    'blast=s',          'ortholog-table=s',
    'ortholog-cds=s',   'species-order=s',
    'cpus=i',           'help',
    'taxon-assemble=s', 'taxon-ortholog=s',
    "genscan-matrix=s", 'predictions!',
    'combineonly',      'tgicl=s',
    'clipping!'
);

pod2usage( -verbose => 2 ) if $opt{help};

$debug = 1 if ( $opt{debug} );

if ( $opt{combineonly} && $opt{"output-file"} && $opt{"output-dir"} ) {
    &combine;
    exit;
}

# MANDATORY PARAMETER

my @mandatory_parameter =
  ( "assignment", "reference", "trinity", "blast", "read-index", "output-dir" );
for (@mandatory_parameter) {
    if ( not exists $opt{$_} ) {
        die "\nMissing parameter: $_\n\n";
    }
}

my @missing_files;
push @missing_files, "Assignment table not found."
  unless ( -e $opt{assignment} );
push @missing_files, "Reference transcriptome not found."
  unless ( -e $opt{reference} );
push @missing_files, "Trinity output not found." 
  unless ( -e $opt{trinity} );
push @missing_files, "Blast results for CDS annotation not found."
  unless ( -e $opt{blast} );
push @missing_files, "GENSCAN matrix not found."
  unless ( -e $opt{'genscan-matrix'} );

if ( $opt{'clipping'} ) {
    push @missing_files, "Hash containing path to readfiles not found."
      unless ( -e $opt{'read-index'} );
}

# OPTIONAL PARAMETER
if ( exists $opt{'ortholog-table'} ) {
    push @missing_files, "Ortholog table not found."
      unless ( -e $opt{'ortholog-table'} );
    push @missing_files, "Ortholog CDS sequences not found."
      unless ( -e $opt{'ortholog-cds'} );
}

if (@missing_files) {
    print STDERR "ERROR: Missing files:\n";
    for (@missing_files) {
        print STDERR "\t$_\n";
    }
    die;
}

# OPTIONAL PARAMETER

unless ( $opt{scaffolding} ) {
    print "Scaffolding disabled\n";
} else {
    push @missing_files, "Blast result for scaffolding not found."
      unless ( -e $opt{scaffolding} );
}

my @species_order = split ",", $opt{'species-order'}
  if ( $opt{'species-order'} );

my $quiet = 0;

# TAXON OBJECTS
my ( $taxon, $taxon_id, $taxon_name, $division ) =
  getTaxon( $opt{'taxon-assemble'} );
my ( undef, undef, $orth_taxon_name, undef ) =
  getTaxon( $opt{'taxon-ortholog'} );
$orth_taxon_name =~ s/\s/_/g;

# READING/INDEXING INPUT FILES
print "Reading reference genbank\n" if ( !$quiet );
my $db_orth =
  GenbankHelper::getIndex( $opt{reference}, -reindex => $opt{reindex} );
die("\nError during indexing of: $opt{reference}\n\n") unless ($db_orth);

print "Reading trinity contigs\n" if ( !$quiet );
my $db_trinity =
  Bio::DB::Fasta->new( $opt{trinity}, -reindex => $opt{reindex} );
die("\nError during indexing of: $opt{trinity}\n\n") unless ($db_trinity);

print "Reading blast result for CDS annotation\n" if ( !$quiet );
my $blast_all = getBlastHash( $opt{blast} );
die("\nError reading blast results for CDS annotation.\n\n")
  unless ($blast_all);

print "Reading blast result for scaffolding fragments\n" if ( !$quiet );
my $blast = getScaffoldFragments( $opt{scaffolding} );
die("\nError reading blast results for scaffolding.\n\n") unless ( !$quiet );

my ( $ortholog_table, $taxon_names );
my $ortholog_table_db;
if ( exists $opt{'ortholog-table'} ) {
    print "Reading ortholog table\n" if ( !$quiet );
    ( $ortholog_table, $taxon_names ) = indexTable( $opt{'ortholog-table'} );
    die("\nError reading ortholog table.\n\n")
      unless ( $ortholog_table && $taxon_names );

    print "Indexing ortholog CDS fasta\n" if ( !$quiet );
    $ortholog_table_db =
      Bio::DB::Fasta->new( $opt{'ortholog-cds'}, -reindex => $opt{reindex} );
    die("\nError reading ortholog table.\n\n") unless ($ortholog_table_db);
}

my $component2readfile;
if ( $opt{'clipping'} ) {
    print "Retrieving hash: component => readfile\n" if ( !$quiet );
    eval { $component2readfile = retrieve $opt{'read-index'}; };
    if ($@) { die "\nCould not retrieve read file hash. Rebuild!\n\n"; }
}

my $tgicl;
if ( exists $opt{'tgicl'} ) {
    print "Indexing cluster members (TGICL)\n" if ( !$quiet );
    $tgicl = indexTGICL( $opt{'tgicl'} );
}

# SUBSET OF ORTHOLOGS/TRANSCRIPTS ONLY

# provided by argument
my %contigs_only = map { $_ => 1 } split ",", $opt{contig} if ( $opt{contig} );
my %orthologs_only = map { $_ => 1 } split ",", $opt{ortholog}
  if ( $opt{ortholog} );

# provided by file
if ( $opt{"contig-file"} && -e $opt{'contig-file'} ) {
    print "Reading target contigs\n" if ( !$quiet );
    open my $fh, "<", $opt{"contig-file"} or die $!;
    while (<$fh>) { chomp; $contigs_only{$_} = 1; }
    close $fh;
}
if ( $opt{"ortholog-file"} && -e $opt{'ortholog-file'} ) {
    print "Reading target orthologs\n" if ( !$quiet );
    open my $fh, "<", $opt{"ortholog-file"} or die $!;
    while (<$fh>) { chomp; $orthologs_only{$_} = 1; }
    close $fh;
}

print "Starting...\n";

# RUN EACH ASSIGNMENT IN PARALLEL
my $pm = Parallel::ForkManager->new( $opt{cpus} );

# START
my $t0    = Benchmark->new();
my $count = 0;
open FH, "<", $opt{assignment} or die $!;
while ( defined( my $line = <FH> ) ) {
    next if ( $line =~ /^#/ );

    my @e = split "\t", $line;

    unless ( @e >= 4 ) {
        print STDERR "ERROR: $opt{assignment}: line $. invalid. Skipping.\n";
        next;
    }

    next if ( %contigs_only   && not exists $contigs_only{ $e[TRIN_ID] } );
    next if ( %orthologs_only && not exists $orthologs_only{ $e[ORTH_ID] } );

    $count++;

    $pm->start and next;
    process(
        (
            $e[TRIN_ID],  $e[ORTH_ID], $e[STRAND],
            $e[HIT_TYPE], $e[STARTS],  $e[ENDS]
        )
    );
    $pm->finish;
}
close FH;
$pm->wait_all_children;
my $t1 = Benchmark->new();

print "Processing $count transcripts took  "
  . timestr( timediff( $t1, $t0 ) ) . "\n";
&combine;

sub combine {

    # COMBINE RESULTS TO SINGLE GENBANK
    if ( $opt{'output-file'} ) {
        remove_tree( $opt{'output-file'} ) if ( -e $opt{'output-file'} );
        $t0 = Benchmark->new();
        print "Combining results...\n" if ($debug);
        my @command = (
            "$FindBin::Bin/combine_fast.sh",
            $opt{'output-dir'}, "_final.gbk", $opt{'output-file'},
            "&> /dev/null"
        );
        system( join " ", @command );
        $t1 = Benchmark->new();
        print "Combining files took " . timestr( timediff( $t1, $t0 ) ) . "\n";
    }
}

sub process {
    my ( $trin_id, $orth_id, $strand, $hit_type ) = @_;

    my $t0 = Benchmark->new();
    print "Processing: $trin_id [$orth_id,$hit_type]..\n" if ($debug);

    my ( $type, $fusion ) = map { uc } split ",", $hit_type;

    my %obj = (
        contig_id      => $trin_id,
        orth_id        => $orth_id,
        contig_seq     => undef,
        orth_seq       => undef,
        orth_file      => undef,
        transcript_dir => undef,
        fragments      => undef,
        cds_start      => undef,
        cds_end        => undef,
        revcom         => 0,
    );

    # get contig
    my $db_contig_seq = $db_trinity->get_Seq_by_acc($trin_id);
    if (
        !(
            defined $db_contig_seq
            && ref $db_contig_seq eq "Bio::PrimarySeq::Fasta"
        )
      )
    {
        print STDERR "ERROR: Contig not found in database! Skipping: "
          . $trin_id
          . " [$orth_id]\n";
        return;
    }

    # get ortholog
    # MSG: Unrecognized DBSOURCE data: BioProject: PRJNA175699 ...
    $obj{orth_seq} = $db_orth->get_Seq_by_acc($orth_id);
    unless (
        !(
            defined $obj{orth_seq}
            && ref $obj{orth_seq} eq "Bio::PrimarySeq::Fasta"
        )
      )
    {
        print STDERR "ERROR: Transcript not found in database! Skipping: "
          . $trin_id
          . " [$orth_id]\n";
        return;
    }

    # get symbol
    my ($orth_gene_feature) = $obj{orth_seq}->get_SeqFeatures('gene');
    ( $obj{orth_sym} ) = $orth_gene_feature->get_tag_values('gene')
      if ($orth_gene_feature);

    unless ( defined $obj{orth_sym} ) {
        $obj{orth_sym} = "NOT_DEFINED";
        print STDERR "ERROR: Transcript has no gene-symbol assigned.\n";
    }

    # transcript folder
    $obj{transcript_dir} = getTranscriptDir( $obj{orth_seq} );
    createTranscriptDir( $obj{transcript_dir} );

    # fasta file with ortholog
    $obj{orth_file} = catfile( $obj{transcript_dir}, "ortholog.fa" );
    my $io =
      Bio::SeqIO->new( -file => ">" . $obj{orth_file}, -format => "fasta" );
    $io->write_seq( $obj{orth_seq} );

    # orientated! "clone" of contig
    if ( $strand eq "-1" ) {
        my $db_contig_seq = $db_contig_seq->revcom;
        $obj{contig_seq} = Bio::Seq->new(
            -seq => $db_contig_seq->seq,
            -id  => $db_contig_seq->id
        );
        $obj{revcom} = 1;
    } else {
        $obj{contig_seq} = Bio::Seq->new(
            -seq => $db_contig_seq->seq,
            -id  => $db_contig_seq->id
        );
    }

    # collect all features
    my %features;

    # SCAFFOLDING
    if ($blast) {
        if ( $blast->{ $obj{orth_id} } ) {

            my ( $scaffold_seq, $scaffold_features, $fragment_ids ) =
              performScaffolding( \%obj );

            if ( defined $scaffold_seq && defined $scaffold_features ) {
                $obj{fragments}  = $fragment_ids;
                $obj{contig_seq} = Bio::Seq->new(
                    -id  => $obj{contig_id},
                    -seq => $scaffold_seq->seq
                );
                push @{ $features{misc_feature} }, @{$scaffold_features};
            }
        }
    }

    # WRITE CONTIG
    $io = Bio::SeqIO->new(
        -file   => ">" . catfile( $obj{transcript_dir}, "contig.fa" ),
        -format => "fasta"
    );
    $io->write_seq( $obj{contig_seq} );

    # GENSCAN
    $obj{genscan_file} = catfile( $obj{transcript_dir}, "CDS_genscan.txt" );
    unless ( -e $obj{genscan_file} ) {
        my @command = (
            "genscan", $opt{'genscan-matrix'},
            catfile( $obj{transcript_dir}, "contig.fa" ),
            ">", $obj{genscan_file},
            "2>", "$obj{genscan_file}.err"
        );
        my $failed = system( join " ", @command );
        if ($failed) {
            print STDERR "ERROR: Genscan failed: "
              . $obj{contig_seq}->id
              . " [$obj{orth_id}]. Skipping assignment.\n";
            print STDERR join " ", @command;
            print STDERR "\n";
            return;
        }
    }

    # CDS PREDICTION
    my ( $cds_feature, @seleno_features ) = performCDSprediction( \%obj );

    if ( defined $cds_feature ) {
        $features{CDS}  = $cds_feature;
        $obj{cds_start} = $cds_feature->start;
        $obj{cds_end}   = $cds_feature->end;
        push @{ $features{misc_feature} }, @seleno_features
          if (@seleno_features);
    } else {
        print STDERR "ERROR: No CDS predicted. Skipping: "
          . $obj{contig_seq}->id . " ["
          . $obj{orth_id} . "]\n";
        return;
    }

    # mRNA CLIPPING
    my $clipped = "clipped:0";
    my ( $start, $stop ) = ( 1, $obj{contig_seq}->length );
    my $promoter = undef;
    if ( $opt{'clipping'} ) {
        my ( $gene, $clippedB, $promoter, $polyA_signal, $polyA_site ) =
          performClipping( \%obj );
        if ($gene) {

            ( $start, $stop ) = ( $gene->[0], $gene->[1] );

            # trust CDS prediction and drop 5' clipping prediction
            if ( $start > $cds_feature->start ) {
                $start    = 1;
                $promoter = undef;
            }

            if ( $start == 1 && $obj{"contig_seq"}->length == $stop ) {
                $clippedB = 0;
            }

            $clipped = "clipped:1" if ($clippedB);
        }
    }

    # Feature: SOURCE
    $features{source} = Bio::SeqFeature::Generic->new(
        -primary_tag => "source",
        -start       => "1",
        -end         => $obj{contig_seq}->length,
        -tag         => {
            organism => $taxon_name,
            db_xref  => "taxon:" . $taxon_id,
            mol_type => $obj{orth_seq}->molecule,
            note     => $type . ":ortholog:$orth_id",
        }
    );

    # Feature: PROMOTER
    $features{promoter} = Bio::SeqFeature::Generic->new(
        -primary_tag => "promoter",
        -start       => $promoter->[0],
        -end         => $promoter->[1],
        -tag         => { gene => $obj{orth_sym}, note => "genscan" }
    ) if ( defined $promoter );

    # Feature: GENE
    if ( ref $features{CDS}->location eq "Bio::Location::Split" ) {
        my $tmp_split_location = Bio::Location::Split->new();

        my @sub_locations =
          sort { $a->start <=> $b->start }
          $features{CDS}->location->sub_Location();

        # adjust CDS first and last exons start and end coordinates, respectively
        my $first_exon = ( shift @sub_locations )->clone;
        my $last_exon  = ( pop @sub_locations )->clone;

        $first_exon->start($start);
        $last_exon->end($stop);

        $tmp_split_location->add_sub_Location($first_exon);
        for (@sub_locations) { $tmp_split_location->add_sub_Location($_); }
        $tmp_split_location->add_sub_Location($last_exon);

        $features{mRNA} = Bio::SeqFeature::Generic->new(
            -primary_tag => "mRNA",
            -location    => $tmp_split_location,
            -tag         => { gene => $obj{orth_sym}, note => $clipped }
        );

    } else {
        $features{mRNA} = Bio::SeqFeature::Generic->new(
            -primary_tag => "mRNA",
            -start       => $start,
            -end         => $stop,
            -tag         => { gene => $obj{orth_sym}, note => $clipped }
        );
    }

    my $accession = $obj{contig_id};
    $accession .= "_$fusion" if ($fusion);

    #
    # DEFINITION - A concise description of the sequence. Mandatory keyword/one or more records.
    #
    # Example:
    #   DEFINITION  Heterocephalus glaber (BHMT) mRNA
    my $description =
      "$taxon_name ($obj{orth_sym}) " . $obj{orth_seq}->molecule;

    my $annotated_contig;
    if ($taxon) {
        $annotated_contig = Bio::Seq::RichSeq->new(
            -seq              => $obj{contig_seq}->seq,
            -id               => $accession,
            -accession_number => $accession,
            -division         => $division,
            -molecule         => $obj{orth_seq}->molecule,
            -dates            => TDATE,
            -desc             => $description
        );
    } else {
        $annotated_contig = Bio::Seq::RichSeq->new(
            -seq              => $obj{contig_seq}->seq,
            -id               => $accession,
            -accession_number => $accession,
            -division         => $division,
            -molecule         => $obj{orth_seq}->molecule,
            -dates            => TDATE,
            -desc             => $description
        );
    }

    print "Adding features...\n" if ($debug);

    for (ORDER_FEATURES) {
        my $value = $features{$_};
        if ($value) {
            if ( ref $value eq "ARRAY" ) {
                for my $f ( @{ $features{$_} } ) {
                    $annotated_contig->add_SeqFeature($f);
                }
            } else {
                $annotated_contig->add_SeqFeature($value);
            }
        }
    }

    print "Writing GenBank\n" if ($debug);
    $io = Bio::SeqIO->new(
        -file   => ">" . catfile( $obj{transcript_dir}, "_final.gbk" ),
        -format => "GenBank"
    );
    $io->write_seq($annotated_contig);

    my $te = Benchmark->new();
    print "$obj{contig_id}\t$obj{orth_id}\t$obj{orth_sym}\t"
      . timestr( timediff( $te, $t0 ) ) . "\n";
}

sub indexTGICL {
    my ($file) = @_;

    my %tgicl;
    if ( -e $file ) {
        open my $fh, "<", $file or die "Can't open file for reading: $file\n";
        while (<$fh>) {
            chomp;
            if (/^(CL\d+Contig\d+)\t(.+)/) {
                my @elements = split ",", $2;
                my $clusterid = $1;
                my %components;
                for (@elements) {
                    $components{$1} = 1 if (/^(.+?)[\._]/);
                }
                $tgicl{$clusterid} = [ keys %components ];
            }
        }
        close $fh;
    }
    return \%tgicl;
}

sub getReadFiles {
    my @ids = @_;

    my %components;
    for my $frag_id (@ids) {
        if ( $frag_id =~ /^(.+?)[\._]/ ) {

            # Example: c121735.graph_c0_seq1
            $components{$1} = 1;
        } elsif ( defined $tgicl && $frag_id =~ /^CL(\d+)Contig(\d+)/ ) {

            # TGICL Clusters
            for my $member ( @{ $tgicl->{$frag_id} } ) {
                $components{$member} = 1;
            }
        }
    }

    my @readfiles;
    for ( keys %components ) {
        my $readfile = $component2readfile->{$_};
        unless ($readfile) {
            print STDERR "ERROR: Readfile not found for component: $_\n";
        } else {
            push @readfiles, $readfile;
        }
    }

    return join ",", @readfiles;
}

sub getTranscriptDir {
    my ($orth_seq) = @_;

    # get symbol
    my ($orth_gene_feature) = $orth_seq->get_SeqFeatures('gene');
    my ($orth_sym)          = $orth_gene_feature->get_tag_values('gene');
    if ( $orth_sym =~ /(.+?)\s+/ ) {
        $orth_sym = $1;
    }

    return catdir( $opt{'output-dir'}, $orth_sym . "_" . $orth_seq->id );
}

sub createTranscriptDir {
    my ($transcript_dir) = @_;

    print "Output path: $transcript_dir\n" if ($debug);
    if ( -d $transcript_dir && $opt{new} ) {
        print "Removing previously computed results\n" if ($debug);
        remove_tree($transcript_dir);
        make_path($transcript_dir);
    } elsif ( !-d $transcript_dir ) {
        make_path($transcript_dir);
    }

    return 1;
}

sub performClipping {
    my ($obj) = @_;

    my $cds_file   = catfile( $obj->{transcript_dir}, "3UTR_known_CDS.csv" );
    my $blast_file = catfile( $obj->{transcript_dir}, "3UTR_blast.csv" );
    my $output_genscan =
      catfile( $obj->{transcript_dir}, "3UTR_annotation.csv" );
    my $output_clip = catfile( $obj->{transcript_dir}, "3UTR_sumscore.txt" );
    my $output_clip_err = catfile( $obj->{transcript_dir}, "3UTR.err" );

    my $orth_file   = catfile( $obj->{transcript_dir}, "3UTR_ortholog.fa" );
    my $contig_file = catfile( $obj->{transcript_dir}, "3UTR_contig.fa" );

    # assign symbols to CDS regions (genscan_annotation.pl)
    if ( !-e $output_genscan ) {
        print "Annotated genscan predictions..\n" if ($debug);

        # prepare CDS annotation table
        open my $fh, ">", $cds_file or die $!;
        print $fh join "\t", $obj->{cds_start}, $obj->{cds_end}, 1,
          $obj->{orth_sym}, $obj->{orth_id};
        print $fh "\n";
        close $fh;

        # prepare hit table
        open $fh, ">", $blast_file or die $!;
        for ( @{ $blast_all->{ $obj->{contig_id} } } ) {
            print $fh join "\t", @$_;
            print $fh "\n";
        }
        close $fh;

        # BUILD COMMAND
        my @command = (
            "$^X $FindBin::Bin/genscan_annotation.pl",
            "-cds $cds_file",
            "-hits $blast_file",
            "-genscan",
            $obj->{genscan_file},
            "> $output_genscan"
        );

        # RUN
        my $failed = system( join " ", @command );
        if ($failed) {
            print STDERR
              "ERROR: genscan_annotate.pl failed. Ignoring further CDS regions if any. ["
              . $obj->{contig_id} . ",",
              . $obj->{orth_sym} . "]\n";
        }
        remove_tree($cds_file);
    }

    # calculate clipping positions for each CDS region
    if ( !-e $output_clip ) {
        print "3 UTR clipping..\n" if ($debug);

        # get accessions of orthologous transcripts
        my @transcripts;
        open my $fh, "<", $output_genscan or die $!;
        while (<$fh>) {
            next if (/^#/);
            ( undef, undef, my $transcript ) = split;
            push @transcripts, $transcript;
        }
        close $fh;

        # prepare fasta with sequences of orthologs
        my $io =
          Bio::SeqIO->new( -file => ">" . $orth_file, -format => "fasta" );
        for (@transcripts) {
            next if ( $_ eq "UNKNOWN" );
            my $tmp_orth_seq = $db_orth->get_Seq_by_acc($_);
            $io->write_seq($tmp_orth_seq) if ($tmp_orth_seq);
            die "Accession not found in database: $_\n" unless ($tmp_orth_seq);
        }

        # fasta with contig
        $io =
          Bio::SeqIO->new( -file => ">" . $contig_file, -format => "fasta" );
        $io->write_seq( $obj->{contig_seq} );

        # get reads from scaffolding fragments
        my $readfiles =
          getReadFiles( $obj->{contig_id}, @{ $obj->{fragments} } );
        unless ($readfiles) {
            print STDERR "ERROR: No read file found: "
              . $obj->{contig_id} . "\n";
            return;
        }

        # BUILD COMMAND
        my @command = (
            "$^X $FindBin::Bin/clipping.pl",
            "-contig $contig_file",
            "-ortholog $orth_file",
            "-input $output_genscan",
            "-readfile $readfiles",
            "-utr-length ",
            UTR_LENGTH,
            "-output_dir",
            catdir( $obj->{transcript_dir}, "3UTR" ),
            "-prefix 3UTR"
        );
        push @command, "-debug" if ($debug);
        push @command, "1> $output_clip";
        push @command, "2> $output_clip_err";

        my $failed = system( join " ", @command );
        if ($failed) {
            print STDERR "ERROR: mRNA clipping failed.\n";
            return;
        }
    } else {
        print "Clip scores already computed. Skipping.\n" if ($debug);
    }

    # PARSE OUTPUT
    my ( $gene, $promotor, $clipped );
    open my $fh, "<", $output_clip or die $!;
    while (<$fh>) {
        next if (/^#/);
        my @e = split;
        next unless ( @e == 12 );
        if ( $e[10] eq $obj->{orth_id} ) {
            $gene = [ $e[8], $e[9] ];
            $clipped = 1 if ( $e[8] > 1 || $e[9] < $obj->{contig_seq}->length );
            if ( $e[4] != -1 && $e[4] ne "inf" ) {
                $promotor = [ $e[4], $e[5] ];
            }
        }
    }
    close $fh;

    return ( $gene, $clipped, $promotor );
}

sub performCDSprediction {
    my ($obj) = @_;

    my $out_aln_file = catfile( $obj->{transcript_dir}, "CDS_alignment.aln" );
    my $out_txt_file = catfile( $obj->{transcript_dir}, "CDS_result.txt" );
    my $out_genscan_file = $obj->{genscan_file};

    if ( -e $out_aln_file && -e $out_txt_file ) {
        print "CDS prediction already performed. Skipping.\n" if ($debug);
    } else {
        print "CDS prediction..\n" if ($debug);

        my $ortholog_file =
          catfile( $obj->{transcript_dir}, "CDS_ortholog_cds.fa" );
        my $io =
          Bio::SeqIO->new( -file => ">" . $ortholog_file, -format => "fasta" );

        # adding full length ortholog (better 5' alignments)
        $io->write_seq(
            Bio::Seq->new(
                -id  => "full:" . $orth_taxon_name . ":" . $obj->{orth_id},
                -seq => $obj->{orth_seq}->seq
            )
        );

        my $orthologs = $ortholog_table->{ $obj->{orth_id} }
          if ( defined $ortholog_table );
        if ($orthologs) {

            # fasta with CDS of orthologs
            for ( my $i = 0; $i < @$orthologs; $i++ ) {
                next if ( $orthologs->[$i] eq "NA" );
                my $seq =
                  $ortholog_table_db->get_Seq_by_acc( $orthologs->[$i] );
                unless ($seq) {
                    print STDERR "Sequence not found: "
                      . $orthologs->[$i]
                      . ". Ortholog of "
                      . $obj->{orth_id} . "\n";
                    next;
                }

                my $last_codon = uc( substr $seq->seq, -3 );
                my $termination =
                  (      $last_codon eq "TAG"
                      || $last_codon eq "TAA"
                      || $last_codon eq "TGA" );

                $io->write_seq(
                    Bio::Seq->new(
                        -id  => "cds:" . $taxon_names->[$i] . ":" . $seq->id,
                        -seq => $termination
                        ? substr $seq->seq,
                        0, -3
                        : $seq->seq
                    )
                );
            }
        } else {
            print "Using assigned ortholog for CDS prediction only\n"
              if ($debug);
            my ($cds_feature) = $obj->{orth_seq}->get_SeqFeatures("CDS");
            my $last_codon = uc( substr $cds_feature->spliced_seq->seq, -3 );
            my $termination =
              (      $last_codon eq "TAG"
                  || $last_codon eq "TAA"
                  || $last_codon eq "TGA" );

            $io->write_seq(
                Bio::Seq->new(
                    -id  => "cds:" . $orth_taxon_name . ":" . $obj->{orth_id},
                    -seq => $termination
                    ? substr $cds_feature->spliced_seq->seq,
                    0, -3
                    : $cds_feature->spliced_seq->seq
                )
            );
        }

        # contig file
        my $contig_file = catfile( $obj->{transcript_dir}, "CDS_contig.fa" );
        $io =
          Bio::SeqIO->new( -file => ">" . $contig_file, -format => 'fasta' );
        $io->write_seq( $obj->{contig_seq} );

        # get number of selenocysteines
        my $selenopos = 0;
        for ( $obj->{orth_seq}->get_SeqFeatures('misc_feature') ) {
            if ( $_->has_tag('note') ) {
                for my $note ( $_->get_tag_values('note') ) {
                    $selenopos++ if ( $note =~ /Selenocysteine/ );
                }
            }
        }

        my @command = (
            "$^X $FindBin::Bin/CDS_predict.pl", "-contig $contig_file",
            "-compareTo",                       $obj->{orth_id},
            "-seleno $selenopos",               "-out-aln $out_aln_file",
            "-out-genscan $out_genscan_file",   "-unaligned",
            "-msa",                             $ortholog_file,
            "-msa-format",                      "fasta",
            "-genscan-matrix",                  $opt{'genscan-matrix'}
        );
        push @command, "-predictions" if ( $opt{'predictions'} );

        push @command, "> $out_txt_file 2> $out_txt_file.err";
        print STDERR ( join " ", @command ) if ($debug);
        my $failed = system( ( join " ", @command ) );
        if ($failed) {
            print STDERR "ERROR: CDS prediction failed! ["
              . $obj->{contig_id} . "; "
              . $obj->{orth_id} . "; "
              . $obj->{orth_sym} . "]\n";
            print STDERR "$!";
            return undef;
        }
    }

    unless ( -e $out_txt_file ) {
        print STDERR "ERROR: CDS prediction failed! ["
          . $obj->{contig_id} . "; "
          . $obj->{orth_id} . "; "
          . $obj->{orth_sym} . "]\n";
        return undef;
    }

    my @features = parseCDSpredict($out_txt_file);

    # add infered CDS to alignment
    my $io    = Bio::AlignIO->new( -file => $out_aln_file, -format => "fasta" );
    my $aln   = $io->next_aln;
    my $seq   = $aln->get_seq_by_id( $obj->{contig_id} );
    my $left  = $seq->column_from_residue_number( $features[0]->start );
    my $right = $seq->column_from_residue_number( $features[0]->end );
    my $cds_slice = $seq->subseq( $left, $right );

    my $new_cds_seq =
        ( "-" x ( $left - 1 ) )
      . $cds_slice
      . ( "-" x ( $aln->length - $right ) );
    $aln->add_seq(
        Bio::LocatableSeq->new(
            -seq => $new_cds_seq,
            -id  => "cds:" . $seq->id
        )
    );

    $io = Bio::AlignIO->new( -file => ">" . $out_aln_file, -format => "fasta" );
    $io->write_aln($aln);

    return @features;
}

sub parseCDSpredict {
    my ($out_txt_file) = @_;

    my @features;

    my ( $method, $location, $completeness, $orthologs, $prot, $seleno );
    if ( -e $out_txt_file ) {

        open my $fh, "<", $out_txt_file or die $!;
        while (<$fh>) {
            chomp;
            next if (/^#/);
            my @e = split;
            next unless ( @e == 6 );

            ( $method, $location, $completeness, $orthologs, $prot, $seleno ) =
              @e;
            last;
        }
        close $fh;

        unless ( defined $method ) {
            print STDERR "ERROR: Could not find any CDS!\n";
            return;
        }
    }

    my $inference_string;
    if ( $method eq "genscan" ) {
        $inference_string = "genscan";
    } elsif ( $method eq "alignment" ) {
        if ( $orthologs =~ /(.+?:.+?)(,|$)/ ) {
            $inference_string = "alignment:$1";
        }
    } elsif ( $method eq "longest_orf" ) {
        $inference_string = "longest_orf";
    } else {
        $inference_string = "UNKNOWN";
    }

    # CDS feature
    my $cds_feature = Bio::SeqFeature::Generic->new(
        -primary_tag => 'CDS',
        -location    => getLocationObject( $completeness, $location ),
        -tag         => {
            translation => $prot,
            inference   => $inference_string,
            note        => $completeness,
        }
    );
    push @features, $cds_feature;

    if ( $seleno ne "NA" ) {
        my @starts = split ",", $seleno;
        for my $seleno_start (@starts) {
            my $seleno_end = $seleno_start + 2;

            push @features,
              Bio::SeqFeature::Generic->new(
                -primary_tag => 'misc_feature',
                -start       => $seleno_start,
                -end         => $seleno_end,
                -tag         => {
                    note => "Selenocysteine"
                }
              );
            $cds_feature->add_tag_value( 'transl_except',
                "(pos:$seleno_start..$seleno_end,aa:Sec)" );
        }

    }
    return @features;
}

sub getLocationObject {
    my ( $completeness, $location ) = @_;

    my @locations = split ",", $location;

    my $loc;
    if ( @locations == 1 ) {
        $loc = getLocation( $completeness, "5end3end", $locations[0] );
    } else {
        $loc = Bio::Location::Split->new();
        my $start_loc = shift @locations;
        my $end_loc   = pop @locations;
        $loc->add_sub_Location(
            getLocation( $completeness, "5end", $start_loc ) );
        for (@locations) {
            $loc->add_sub_Location( getLocation( $completeness, "exact", $_ ) );
        }
        $loc->add_sub_Location(
            getLocation( $completeness, "3end", $end_loc ) );
    }
    return $loc;
}

sub getLocation {
    my ( $completeness, $end, $string ) = @_;

    my ( $start, $stop ) = split "-", $string;
    my $loc = Bio::Location::Fuzzy->new( -start => $start, -end => $stop );

    if ( $end eq "exact" || !defined $completeness ) {
        $loc->start_pos_type('EXACT');
        $loc->end_pos_type('EXACT');
        return $loc;
    }

    if ( $completeness =~ /partial/ ) {
        $loc->start_pos_type('BEFORE') if ( $end =~ /5end/ );
        $loc->end_pos_type('AFTER')    if ( $end =~ /3end/ );
        if ( $completeness =~ /5prime/ && $end =~ /5end/ ) {
            $loc->end_pos_type('EXACT');
        }
        if ( $completeness =~ /3prime/ && $end =~ /3end/ ) {
            $loc->start_pos_type('EXACT');
        }
    }

    return $loc;
}

sub indexTable {
    my ($file) = @_;

    my %ortholog_table;
    my @names;

    my @column_order;

    open my $fh, "<", $file or die $!;
    while (<$fh>) {
        if (/^#(.+)/) {
            my @e = split "\t", $1;
            for (@e) {
                ( undef, undef, my $name ) = getTaxon($_);
                $name = "UNKNOWN" unless ($name);
                $name =~ s/\s/_/g;
                push @names, $name;
            }

            # apply required order of orthologs
            my %numbers = map { $_ => 1 } grep { $_ < @e } @species_order;
            my @missing;
            for ( 0 .. $#e ) {
                push @missing, $_ unless ( exists $numbers{$_} );
            }
            @column_order = ( @species_order, @missing );
        } else {
            my @e = split;
            $ortholog_table{ $e[0] } = [ @e[@column_order] ];
        }
    }
    close $fh;
    @names = @names[@column_order];

    if ($debug) {
        print "Order of species: " . join ", ", @names;
        print "\n";
    }

    return \%ortholog_table, \@names;
}

sub performScaffolding {
    my ($obj) = @_;

    my $output_contig_file =
      catfile( $obj->{transcript_dir}, "SCAFFOLD_scaffold.fa" );
    my $output_txt_file =
      catfile( $obj->{transcript_dir}, "SCAFFOLD_result.txt" );

    if ( -e $output_txt_file ) {
        print "Scaffolded already performed. Skipping.\n" if ($debug);
    } else {
        print "Scaffolding..\n" if ($debug);
        my $ortholog_file = catfile( $obj->{transcript_dir}, "ortholog.fa" );
        my $contig_file =
          catfile( $obj->{transcript_dir}, "SCAFFOLD_contig.fa" );
        my $fragment_file =
          catfile( $obj->{transcript_dir}, "SCAFFOLD_fragments.fa" );
        my $output_aln_file =
          catfile( $obj->{transcript_dir}, "SCAFFOLD_alignment.aln" );
        my $combined_file =
          catfile( $obj->{transcript_dir}, "SCAFFOLD_temp.aln" );

        #my $out_dir = catdir($obj->{transcript_dir}, "tgicl");
        #make_path($out_dir);

        my $io =
          Bio::SeqIO->new( -file => ">" . $contig_file, -format => "fasta" );
        $io->write_seq( $obj->{contig_seq} );

        $io =
          Bio::SeqIO->new( -file => ">" . $fragment_file, -format => "fasta" );
        for ( @{ $blast->{ $obj->{orth_id} } } ) {
            my $fragment_seq = $db_trinity->get_Seq_by_acc( $_->[0] );

            unless ($fragment_seq) {
                warn("Fragment not found: $_->[0]");
                next;
            }
            $fragment_seq = $fragment_seq->revcom if ( $_->[1] );
            $io->write_seq($fragment_seq);
        }

        my @command = (
            "$^X $FindBin::Bin/scaffolding.pl",
            "-out-combined $combined_file",
            "-ortholog $ortholog_file",
            "-contig $contig_file",
            "-fragments $fragment_file",
            "-out-final $output_aln_file",
            "-output $output_contig_file",
            "> $output_txt_file",
            "2> $output_txt_file.err"
        );

        #print join " ", @command;
        my $failed = system( join " ", @command );
        if ($failed) {
            warn(   "Scaffolding failed: "
                  . $obj->{contig_id} . "/"
                  . $obj->{orth_id}
                  . "\n" );
            return undef;
        }
    }

    my @features;
    my @fragments;
    my $new_seq;
    if ( -e $output_txt_file ) {

        # scaffolding worked, but no resulting scaffold produced
        if ( -e $output_contig_file ) {
            my $io = Bio::SeqIO->new(
                -file   => $output_contig_file,
                -format => "fasta"
            );
            $new_seq = $io->next_seq if ($io);
        }
        return undef unless ($new_seq);

        # scaffolding worked, and we have a new contig!
        open my $fh, "<", $output_txt_file or die $!;
        while (<$fh>) {
            next if ( /^#/ || /^$/ );
            chomp;
            my @e = split "\t";

            if ( $e[0] eq "best" ) {
                push @features,
                  Bio::SeqFeature::Generic->new(
                    -primary_tag => 'misc_feature',
                    -start       => $e[2],
                    -end         => $e[3],
                    -tag         => {
                        note => "best contig:$e[1]"
                    }
                  );
            } elsif ( $e[0] eq "fragment" ) {
                push @features,
                  Bio::SeqFeature::Generic->new(
                    -primary_tag => 'misc_feature',
                    -start       => $e[2],
                    -end         => $e[3],
                    -tag         => {
                        note => "sequence inferred during scaffolding:$e[1]"
                    }
                  );
                push @fragments, $e[1];

            } elsif ( $e[0] eq "gap" ) {
                push @features,
                  Bio::SeqFeature::Generic->new(
                    -primary_tag => 'assembly_gap',
                    -start       => $e[2],
                    -end         => $e[3],
                    -tag         => {
                        gap_type         => "within scaffold",
                        estimated_length => ( $e[3] - $e[2] + 1 ),
                    }
                  )

            }
        }
        close $fh;
    }

    unless (@features) {
        warn( "Something went wrong during feature creation for scaffolding: "
              . $obj->{contig_id} );
        return undef;
    }

    return $new_seq, \@features, \@fragments;

}

sub getTaxon {
    my $id = shift;

    die "Specify TaxonID!" unless ($id);

    #
    # Taxonomy currently requires internet connection.
    #

    my $dbh;
    my ($taxon, $taxon_id, $taxon_name, $division) = (undef, $id, "UNKOWN", "UNKNOWN");
    eval {
        $dbh = Bio::DB::Taxonomy->new(-source => 'entrez');
        $taxon = $dbh->get_taxon( -taxonid => $id );
    };
    if ($@) {
        print STDERR "WARNING: Failed to connect to taxonomy database. No internet access?\n";
    } 
    unless ($taxon) {
        print STDERR "WARNING: Could not retrieven taxon information for ID: $id\n";
    } else {
        $taxon_id   = $taxon->id;
        $taxon_name = $taxon->scientific_name;
        $division   = $taxon->division;
    }
    return ( $taxon, $taxon_id, $taxon_name, $DIVISION{$division} );
}

sub getBlastHash {
    my ($file) = @_;

    my $hash;

    if ( -e "$file.index" && !$opt{reindex} ) {
        $hash = retrieve "$file.index";
    } else {
        open my $fh, "<", $file or die $!;
        while (<$fh>) {
            next if (/^#/);
            my @e = split;
            push @{ $hash->{ $e[TARGET_ID] } },
              [
                $e[QUERY_ID], $e[TARGET_START], $e[TARGET_END],
                $e[IDENTITY], $e[SYMBOL]
              ];
        }
        close $fh;
        store $hash, "$file.index";
    }

    return $hash;
}

sub getScaffoldFragments {
    my ($file) = @_;

    unless ( defined $file && -e $file ) {
        print STDERR
          "Could not find file with BLAST results for scaffolding! No scaffolding!\n";
        return undef;
    }

    my $current_hit = 0;
    my @contigs;
    my $hash = {};

    open FH, "<", $file or die $!;
    while (<FH>) {
        my @e = split;

        $current_hit = $e[TARGET_ID] unless ($current_hit);
        if ( $current_hit ne $e[TARGET_ID] ) {
            $hash->{$current_hit} = [@contigs];
            @contigs              = ();
            $current_hit          = $e[TARGET_ID];
        }

        push @contigs, [
            $e[QUERY_ID],
            $e[STRAND_BLAST] == -1 ? 1 : 0,    # reverse complement?
        ];
    }
    $hash->{$current_hit} = [@contigs];
    close FH;

    return $hash;
}


1;

__END__

