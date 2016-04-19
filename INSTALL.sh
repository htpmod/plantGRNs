### INSTALL required software packages to run this pipeline 

## change if needed 
BINDIR=`pwd` ## $HOME/flowerGRNs/bin 

cd $BINDIR 
## FastQC   
set -x
cd $BINDIR
echo "Installing FastQC ..."
wget -q http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -O fastqc_v0.11.5.zip
unzip -o fastqc_v0.11.5.zip && rm fastqc_v0.11.5.zip
chmod +x FastQC/fastqc
set +x

cd $BINDIR 
## picard   
set -x
cd $BINDIR
echo "Installing picard ..."
wget -q https://github.com/broadinstitute/picard/releases/download/2.2.1/picard-tools-2.2.1.zip -O picard-tools-2.2.1.zip 
unzip -o picard-tools-2.2.1.zip && mv picard-tools-2.2.1 picard && rm picard-tools-2.2.1.zip
set +x


## Fastx-Toolkit (root?)
set -x
cd $BINDIR
wget -q https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz -O libgtextutils-0.7.tar.gz
tar -zxf libgtextutils-0.7.tar.gz && cd libgtextutils-0.7 && ./configure && make && sudo make install
set +x
rm -rf $BINDIR/libgtextutils-0.7* 

set -x
cd $BINDIR
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH
wget -q https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2 -O fastx_toolkit-0.0.14.tar.bz2
tar -xjf fastx_toolkit-0.0.14.tar.bz2 && cd fastx_toolkit-0.0.14
./configure --prefix=$BINDIR/Fastx && make && make install
set +x
rm -rf $BINDIR/fastx_toolkit-0.0.14* 

## samtools 
echo "Checking samtools installation ..." 
if type samtools 2>/dev/null; then
    echo "samtools was already installed. " 
else 
    set -x
    cd $BINDIR
    echo "Installing samtools ..."
    wget -q https://github.com/samtools/samtools/archive/0.1.19.tar.gz -O samtools-0.1.19.tar.gz
    tar -zxf samtools-0.1.19.tar.gz && rm -f samtools-0.1.19.tar.gz && cd samtools-0.1.19
    make && cp samtools $BINDIR
    rm -rf $BINDIR/samtools-0.1.19
    set +x
fi 

## bedtools 
echo "Checking bedtools installation ..." 
if type bedtools 2>/dev/null; then
    echo "bedtools was already installed. " 
else 
    set -x
    cd $BINDIR
    echo "Installing bedtools ..."
    wget -q wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz -O bedtools-2.25.0.tar.gz
    tar -zxf bedtools-2.25.0.tar.gz && rm -f bedtools-2.25.0.tar.gz && cd bedtools2 && make
    cp bin/* $BINDIR 
    rm -rf $BINDIR/bedtools2 
    set +x
fi 

## UCSC toolkits  
set -x
cd $BINDIR
echo "Installing UCSC toolkits ..."
## choose x86_64 platform here 
wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig -O bedGraphToBigWig 
wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed -O bedToBigBed 
chmod +x bedGraphToBigWig bedToBigBed 
set +x

## IDR 
set -x
cd $BINDIR
echo "Installing IDR ..." 
wget -q https://sites.google.com/site/anshulkundaje/projects/idr/idrCode.tar.gz?attredirects=0 -O idrCode.tar.gz
tar -xzf idrCode.tar.gz && rm idrCode.tar.gz
set +x


## phantompeakqualtools 
set -x
cd $BINDIR
echo "Installing phantompeakqualtools ..." 
wget -q https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/phantompeakqualtools/ccQualityControl.v.1.1.tar.gz -O ccQualityControl.v.1.1.tar.gz
tar -xzf ccQualityControl.v.1.1.tar.gz && rm ccQualityControl.v.1.1.tar.gz
set +x

## MACS2
echo "Checking MACS2 installation ..." 
if type macs2 2>/dev/null; then
    echo "MACS2 was already installed. " 
else 
    echo "Installing MACS2 ..."
    set -x
    pip install MACS2 --user
    set +x
fi 

## BiocParallel
echo "Checking BiocParallel installation ..." 
Rscript - <<EOS
    getwd()
    if(!suppressMessages(require(BiocParallel,quietly=T,warn.conflicts=F))) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("BiocParallel", ask=FALSE)
        library(BiocParallel)
    }
EOS

## ChIPseeker
echo "Checking ChIPseeker installation ..." 
Rscript - <<EOS
    getwd()
    if(!suppressMessages(require(ChIPseeker,quietly=T,warn.conflicts=F))) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("ChIPseeker", ask=FALSE)
        library(ChIPseeker)
    }
EOS
