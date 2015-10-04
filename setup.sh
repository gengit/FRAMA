#!/bin/bash

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do 
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" 
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

install_perl() {

printf "\n\nInstalling R-packages in personal library..."

}

install_R() {

printf "\n\nInstalling R-packages in personal library..."

R --vanilla --slave <<END

personal.lib = Sys.getenv("R_LIBS_USER")
if (!file.exists(personal.lib)) { dir.create(personal.lib, rec = T) }

# CRAN packages
mypackages = c("plyr", "ggplot2", "reshape", "gridExtra", "grid")
installme = mypackages[!(mypackages %in% installed.packages()[,"Package"])]
if(length(installme) > 0) {
    install.packages(installme, repos =  "http://cran.us.r-project.org", lib = personal.lib)
} 

# Bioconductor
source("http://bioconductor.org/biocLite.R")
mypackages = c("GO.db", "org.Hs.eg.db", "KEGG.db")
installme = mypackages[!(mypackages %in% installed.packages()[,"Package"])]
if(length(installme) > 0) {
    biocLite(installme, suppressAutoUpdate = T, suppressUpdates = T, lib = personal.lib)
} 

END

}

install_software() {

DIR=$DIR/bin_external

mkdir -p $DIR && cd $DIR

BLAST="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.31+-x64-linux.tar.gz"
RMBLAST="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/LATEST/ncbi-rmblastn-2.2.28-x64-linux.tar.gz"
TRF="http://tandem.bu.edu/trf/downloads/trf407b.linux64"
REPEATMASKER="http://www.repeatmasker.org/RepeatMasker-open-4-0-5.tar.gz"
TGICL="ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/tgicl/tgicl_linux.tar.gz"
BOWTIE1="http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip/download"
BOWTIE2="http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip/download"
MAFFT="http://mafft.cbrc.jp/alignment/software/mafft-7.245-linux.tgz"
GENBLASTA="http://genome.sfu.ca/genblast/latest/genblasta_v1.0.4_linux_x86_64"
BAMTOOLS="https://github.com/pezmaster31/bamtools/archive/v2.4.0.zip"
TRINITY="https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.0.6.tar.gz"
SAMTOOLS="http://sourceforge.net/projects/samtools/files/samtools/1.2/samtools-1.2.tar.bz2/download"
PARALLEL="http://ftp.gnu.org/gnu/parallel/parallel-20150922.tar.bz2"
EMBOSS="http://science-annex.org/pub/emboss/EMBOSS-6.6.0.tar.gz"
CDHIT="https://github.com/weizhongli/cdhit/releases/download/V4.6.4/cd-hit-v4.6.4-2015-0603.tar.gz"

printf "\n\nInstalling software in $DIR\n"

BLAST_BIN=""
if [[ ! $(which blastn 2> /dev/null) ]]; then
    echo -e "--> Downloading NCBI-BLAST\n\n\t$BLAST\n"
    wget --ftp-user=anonymous --ftp-password=anonymous $BLAST
    tar zxf ncbi-blast-2.2.31+-x64-linux.tar.gz 
    rm ncbi-blast-2.2.31+-x64-linux.tar.gz
    BLAST_BIN=$(readlink -f ncbi-blast-2.2.31+/bin)

    ln -s $BLAST_BIN/blastn .
    ln -s $BLAST_BIN/makeblastdb .
else
    BLAST_BIN=$(dirname $(which blastn))
fi

if [[ ! $(which RepeatMasker 2> /dev/null) ]]; then
    echo -e "--> Downloading RMBLAST\n\n\t$RMBLAST\n"
    wget --ftp-user=anonymous --ftp-password=anonymous $RMBLAST
    tar zxf ncbi-rmblastn-2.2.28-x64-linux.tar.gz 
    rm ncbi-rmblastn-2.2.28-x64-linux.tar.gz

    RMBLAST_BIN=$(readlink -f ncbi-rmblastn-2.2.28/bin)

    ln -s $BLAST_BIN/blastx $RMBLAST_BIN/blastx
    ln -s $BLAST_BIN/makeblastdb $RMBLAST_BIN/makeblastdb

    echo -e "--> Downloading Tandem Repeat Finder\n\n\t$TRF\n"
    cat <<EOF
    Tandem Repeats Finder License Terms

    The author of this software grants to any individual or organization the
    right to use and to make an unlimited number of copies of this software.
    You may not de-compile, disassemble, reverse engineer, or modify the
    software. This software cannot be sold, incorporated into commercial
    software or redistributed. The author of this software accepts no
    responsibility for damages resulting from the use of this software and
    makes no warranty or representation, either express or implied, including
    but not limited to, any implied warranty of merchantability or fitness for
    a particular purpose. This software is provided as is, and the user assumes
    all risks when using it.

EOF

    wget $TRF
    mv trf407b.linux64 trf && chmod +x trf
    TRF_BIN=$(readlink -f trf)

    echo -e "--> Downloading RepeatMasker\n\n\t$REPEATMASKER\n"
    wget $REPEATMASKER
    tar xzf RepeatMasker-open-4-0-5.tar.gz 
    rm RepeatMasker-open-4-0-5.tar.gz
    mv RepeatMasker RepeatMasker_4-0-5

    PERL_BIN=$(which perl)
    cd RepeatMasker_4-0-5 
    echo -e "\n$PERL_BIN\n\n$TRF_BIN\n2\n$RMBLAST_BIN\nY\n5\n" | perl configure
    cd ..

    ln -s RepeatMasker_4-0-5/RepeatMasker .
fi

if [[ ! $(which tgicl 2> /dev/null) ]]; then
    echo -e "--> Downloading TGICL\n\n\t$TGICL\n"
    wget --ftp-user=anonymous --ftp-password=anonymous $TGICL
    tar zxf tgicl_linux.tar.gz 
    rm tgicl_linux.tar.gz
fi

if [[ ! $(which bowtie 2> /dev/null) ]]; then
    echo -e "--> Downloading Bowtie1\n\n\t$BOWTIE1\n"
    wget $BOWTIE1
    mv download bowtie1.zip
    unzip -qq bowtie1.zip 
    rm bowtie1.zip

    ln -s bowtie-1.1.2/bowtie .
    ln -s bowtie-1.1.2/bowtie-build .
fi

if [[ ! $(which bowtie2 2> /dev/null) ]]; then
    echo -e "--> Downloading Bowtie2\n\n\t$BOWTIE2\n"
    wget $BOWTIE2
    mv download bowtie2.zip
    unzip -qq bowtie2.zip
    rm bowtie2.zip

    ln -s bowtie2-2.2.6/bowtie2 .
    ln -s bowtie2-2.2.6/bowtie2-build .
fi

if [[ ! $(which mafft 2> /dev/null) ]]; then
    echo -e "--> Downloading MAFFT\n\n\t$MAFFT\n"
    wget $MAFFT
    tar zxf mafft-7.245-linux.tgz 
    rm mafft-7.245-linux.tgz

    cd mafft-linux64
    ln -s mafft.bat mafft
    cd ..
fi

if [[ ! $(which genblasta 2> /dev/null) ]]; then
    echo -e "--> Downloading GENBLASTA\n\n\t$GENBLASTA\n"
    wget $GENBLASTA
    chmod +x genblasta_v1.0.4_linux_x86_64
    mv genblasta_v1.0.4_linux_x86_64 genblasta
fi

# --------------------------------------------------------------------
# SOURCES 
# --------------------------------------------------------------------

if [[ ! $(which bamtools 2> /dev/null) ]]; then
    echo -e "--> Downloading BamTools\n\n\t$BAMTOOLS\n"
    wget $BAMTOOLS
    mv v2.4.0.zip bamtools.zip
    unzip -qq bamtools.zip 
    rm bamtools.zip
    cd bamtools-2.4.0 
    mkdir build
    cd build 
    cmake ..
    make 
    cd ../..

    ln -s bamtools-2.4.0/bin/bamtools .
fi

if [[ ! $(which cd-hit-est 2> /dev/null) ]]; then
    echo -e "--> Downloading CD-HIT-EST\n\n\t$CDHIT\n"
    wget $CDHIT
    tar zxf cd-hit-v4.6.4-2015-0603.tar.gz
    rm cd-hit-v4.6.4-2015-0603.tar.gz
    cd cd-hit-v4.6.4-2015-0603 
    make 
    cd ..

    ln -s cd-hit-v4.6.4-2015-0603/cd-hit-est .
fi

if [[ ! $(which needle 2> /dev/null) ]]; then
    echo -e "--> Downloading EMBOSS\n\n\t$EMBOSS\n"
    wget $EMBOSS
    tar zxf EMBOSS-6.6.0.tar.gz 
    rm -rf EMBOSS-6.6.0.tar.gz
    cd EMBOSS-6.6.0 && ./configure --without-x && make 
    cd ..

    ln -s EMBOSS-6.6.0/emboss/needle .
    ln -s EMBOSS-6.6.0/emboss/getorf .
fi

if [[ ! $(which parallel 2> /dev/null) ]]; then
    echo -e "--> Downloading Parallel\n\n\t$PARALLEL\n"
    wget $PARALLEL 
    tar xfj parallel-20150922.tar.bz2 
    rm parallel-20150922.tar.bz2
    cd parallel-20150922 
    ./configure --prefix $(pwd) 
    make 
    make install 
    cd ..

    ln -s parallel-20150922/bin/parallel .
fi

if [[ ! $(which samtools 2> /dev/null) ]]; then
    echo -e "--> Downloading SAMTOOLS\n\n\t$SAMTOOLS\n"
    wget $SAMTOOLS
    mv download samtools.tar.bz2
    tar xfj samtools.tar.bz2 
    rm samtools.tar.bz2
    cd samtools-1.2 
    make prefix=$(pwd) 
    cd ..

    ln -s samtools-1.2/samtools .
fi

if [[ ! $(which Trinity.pl 2> /dev/null) ]]; then
    echo -e "--> Downloading Trinity\n\n\t$TRINITY\n"
    wget $TRINITY
    tar zxf v2.0.6.tar.gz 
    rm v2.0.6.tar.gz
    cd trinityrnaseq-2.0.6 
    make 
    ln -s Trinity Trinity.pl
    cd ..
fi
}

clear

printf "This script will install external software that requires (among others)
cmake, zlib, ncurses, R and java. Please install these prerequesites beforehand
(consult README.md).\n\n"
printf "Are you ready to continue? [Y/n]\n"
read -n 1 continue
printf "\n";

if [ "$continue" == "n" ]; then
    printf "Stopped\n"
    exit;
fi

printf "Do you want to install external software? [y/n]\n"
read -n 1 external_software
printf "\nDo you want to install R packages? [y/n]\n"
read -n 1 R_packages
printf "\nDo you want to install Perl modules? [y/n]\n"
read -n 1 perl_modules
printf "\n"

if [ "$external_software" == "y" ]; then
    install_software
fi

if [ "$R_packages" == "y" ]; then
    install_R
fi

if [ "$perl_modules" == "y" ]; then
    install_perl
fi

