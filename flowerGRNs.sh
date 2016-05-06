## Note: be aware that this pipeline only works for sing-end read files 
##       if you are working on paired-end data, some parts of code should 
##       have minor changes (e.g, using different parameters)

# **************** 1. workspace     ****************
# create a workspace under $HOME/plantGRNs
mkdir -p $HOME/plantGRNs && myWorkDIR=$HOME/plantGRNs

# **************** 2. get source code   ****************
# change directory to workspace
cd $myWorkDIR 
# clone the plantGRNs repo from https://github.com/PlantENCODE/
rm -rf pipeline && git clone --recursive https://github.com/PlantENCODE/plantGRNs.git pipeline
# parallelRun is used for parallel analysis 
chmod +x pipeline/parallelRun 
# metadata files were prepared according to Table 1
mv $myWorkDIR/pipeline/meta $myWorkDIR/ 

# **************** 3. installation  **************** 
# create a directory bin store executable programs 
mkdir -p $myWorkDIR/bin 
# install third-party software packages into this directory 
# NOTE: follow the INSTALL.sh for instruction on how to install third-party software packages 
cd $myWorkDIR/bin && sh $myWorkDIR/pipeline/INSTALL.sh
# make executable programs available through the PATH environmental variable
export PATH=$myWorkDIR/bin:$myWorkDIR/pipeline:$PATH

# **************** 4. genome data   **************** 
# create a sub-directory db to save reference genome data 
mkdir -p $myWorkDIR/db && cd $myWorkDIR/db 
# prepare genome sequence FASTA files 
for chrom in $(seq 1 5); do 
    wget -O TAIR10_chr$chrom.fa -q ftp://ftp.arabidopsis.org/Sequences/whole_chromosomes/TAIR10_chr$chrom.fas
done 
cat TAIR10_chr*.fa > tair10.fa && rm -f TAIR10_chr*.fa 
# make bowtie2 index files 
bowtie2-build tair10.fa tair10 
# prepare genome annotation files 
wget -qO- tair10Genes.gff3 ftp://ftp.arabidopsis.org/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff | grep -v "Chr[CM]" > tair10Genes.gff3 
wget -O gene_aliases.txt -q ftp://ftp.arabidopsis.org/Genes/gene_aliases_20140331.txt
awk -vFS="\t" -vOFS="\t" -vTF="$myWorkDIR/meta/TF.list" 'BEGIN{
    while(getline<"gene_aliases.txt"){alias[$1]=$2} while(getline<TF){alias[$1]=$2}} \
    ($3=="gene" && $9~/Note=protein_coding_gene/) \
        {match($9,"ID=([^;$]+)",ID); \
        print $1,$4,$5,ID[1],alias[ID[1]],$7}' \
    tair10Genes.gff3 | sed -e "s/'//g" \
    | sort -k1,1 -k2,2n > tair10Genes.bed 
# miRNA genes 
wget -O tair10miRNAs.gff3 -q ftp://mirbase.org/pub/mirbase/CURRENT/genomes/ath.gff3 
grep -v "#" tair10miRNAs.gff3 | sed -e 's/^chr/Chr/' \
    | awk -vFS="\t" -vOFS="\t" '($3=="miRNA_primary_transcript"){
    match($9,"Name=([^;$]+)",NM);match($9,"ID=([^;$]+)",ID);
    print $1,$4,$5,ID[1],NM[1],$7}' | sed -e 's/ath-//g' \
    | sort -k1,1 -k2,2n > tair10miRNAs.bed 
    
# **************** 5. ChIP-seq data **************** 
# get experiment meta files 
while read SRP; do 
    perl -MLWP::Simple -e "getprint 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=${SRP}&result=read_run&fields=secondary_study_accession,experiment_accession,run_accession,fastq_ftp,sample_alias,scientific_name,library_strategy&download=txt'" | awk '(NR>1)'
done < $myWorkDIR/meta/study.list \
    > $myWorkDIR/meta/experiment.meta 

cd $myWorkDIR/meta
# get metadata for ChIP-seq data sets of interest 
awk -vOFS="\t" 'BEGIN{while(getline < "chip.list"){list[$3]++}} \
    (list[$2]){print $2,$3,"ftp://"$4,$1}' \
    $myWorkDIR/meta/experiment.meta \
    > $myWorkDIR/meta/sample.meta 

# create a sub-directory fq to save ChIP-seq FASTQ files 
mkdir -p $myWorkDIR/fq && cd $myWorkDIR/fq 
# download FASTQ files in parallel (maximum 10 runs)
cut -f 3 $myWorkDIR/meta/sample.meta | xargs -P 10 \
    -I % sh -c 'O=$(basename %); wget -O $O -q %'

# merge multiple SRR files (if any) into a single SRX file 
awk -vOFS="\t" '{a[$1]=a[$1]" "$2".fastq.gz"} \
    END{for(x in a){print x,a[x]}}' \
    $myWorkDIR/meta/sample.meta \
    | while read SRX SRRs; do 
        cat $SRRs > $SRX.fastq.gz
    done 
rm -f SRR*.fastq.gz 

# **************** 6. FastQC    **************** 
cut -f 3 $myWorkDIR/meta/chip.list | sort | uniq > $myWorkDIR/meta/SRX.list 
# save QC report files under fastqc directory  
mkdir -p $myWorkDIR/fastqc && cd $myWorkDIR/fastqc
# be sure that FastQC is available from PATH 
export PATH=$myWorkDIR/bin/FastQC:$PATH
# run FastQC in parallel 
cat $myWorkDIR/meta/SRX.list | while read SRX; do 
        FQ="$SRX.fastq.gz"
        [[ ! -d $myWorkDIR/fastqc/$SRX ]] && \
        echo "mkdir -p $myWorkDIR/fastqc/$SRX; \
        fastqc --quiet --outdir $myWorkDIR/fastqc/$SRX \
        $myWorkDIR/fq/${FQ} "
    done | parallelRun 
## cat $myWorkDIR/meta/SRX.list | xargs -P 10 -I % sh -c 'FQ=%".fastq.gz"; mkdir %; fastqc --quiet --outdir % ../fq/${FQ}'

# **************** 7. Trimming  **************** 
# samples from SRP002328, SRP003928, SRP026163, SRP009053 and SRP020612: quality score >= 25 
awk '(/SRP0(02328|03928|09053|20612|26163)/)' $myWorkDIR/meta/sample.meta \
    | cut -f 1 | sort | uniq \
    | while read SRX; do 
        echo "zcat $myWorkDIR/fq/$SRX.fastq.gz \
            | fastq_quality_filter -Q33 -q 25 -p 80 \
            | fastx_trimmer -z -Q33 -l 35 \
            -o $myWorkDIR/fq/$SRX.clean.fastq.gz "
    done | parallelRun 
# samples from SRP037581: only keep the first 60 bp 
awk '($4=="SRP037581")' $myWorkDIR/meta/sample.meta \
    | cut -f 1 | sort | uniq \
    | while read SRX; do 
        echo "zcat $myWorkDIR/fq/$SRX.fastq.gz \
            | fastx_trimmer -z -Q33 -f 13 -l 60 \
            -o $myWorkDIR/fq/$SRX.clean.fastq.gz "
    done | parallelRun 

# **************** 8. Alignment **************** 
# save alignment result under align directory  
mkdir -p $myWorkDIR/align && cd $myWorkDIR/align
# FASTQ encoding: Sanger/Illumina_1_8 format (--phred33); 
# use 2 threads for alignment (--threads 2); 
# suppress SAM records for unaligned reads (--no-unal)
NTHREADS=2 ## number of threads
cat $myWorkDIR/meta/SRX.list | while read SRX; do 
        FQ="$myWorkDIR/fq/$SRX.fastq.gz"
        [[ -f "$myWorkDIR/fq/$SRX.clean.fastq.gz" ]] && \
        FQ="$myWorkDIR/fq/$SRX.clean.fastq.gz"
        [[ ! -f $SRX.bam ]] && echo "bowtie2 --phred33 \
            --threads $NTHREADS --no-unal --sensitive -k 1 -q \
            -x $myWorkDIR/db/tair10 -U ${FQ} 2> $SRX.bt2.log \
            | samtools view -@ $NTHREADS -bS - \
            | samtools sort -@ $NTHREADS - $SRX && \
            samtools index $SRX.bam "
    done | parallelRun

# **************** 9. Filtering **************** 
#===================================
# 1. Remove unmapped, mate unmapped not primary alignment, 
#    reads failing platform
# 2. Remove low MAPQ reads
# 4. Remove duplicate reads
#===================================
NTHREADS=2 ## number of threads
PICARD=$myWorkDIR/bin/picard ## picard DIR 
cd $myWorkDIR/align
cat $myWorkDIR/meta/SRX.list | while read SRX; do 
        # Filter bam file, based on FLAG 1804: segment unmapped (4) + next segment in the template unmapped (8) + secondary alignments (256) + not passing filters, such as platform/vendor quality controls (512) + PCR or optical duplicates (1024); -q 30 exclude MAPQ < 30 
        [[ ! -f $SRX.final.bam ]] && \
        echo "
            java -jar $PICARD/picard.jar MarkDuplicates \
            I=$SRX.filter.bam O=$SRX.dupmark.bam M=$SRX.dup.qc \
            VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false \
            ASSUME_SORTED=true; samtools view -@ ${NTHREADS} \
            -F 1804 -b $SRX.dupmark.bam > $SRX.final.bam; \
            samtools index $SRX.final.bam "
    done | parallelRun
# merge replicates (if any) for peak calling and visualization 
cd $myWorkDIR/align
awk -vFS="\t" -vOFS="\t" '{bam=$3".final.bam"; \
    if(($2~/_control/)){c=$1"_control"}else{c=$1"_chip"} a[c]=a[c]?a[c]";"bam:bam;} \
    END{for(i in a){print i,a[i]}}' $myWorkDIR/meta/chip.list \
    | while read TF bams; do 
        out=$TF".bam"
        if [[ $bams =~ ";" ]]; then 
            inbams=$(echo $bams | sed -e 's/;/ /g')
            echo "samtools merge $out $inbams; samtools index $out"
        else 
            ln -sf $bams $out;
            ln -sf $bams".bai" $out".bai";
        fi
    done | parallelRun 

# alignment summary 
echo -e "Sample\t#Raw\t#UniqAlign\t%Align\t#Filtered\tFinal" \
    > $myWorkDIR/meta/align.table.txt 
cat $myWorkDIR/meta/SRX.list | while read SRX; do 
        if [[ -f $SRX.filter.bam ]]; then 
            TR=$(grep "reads; of these" $SRX.bt2.log | awk '{print $1}')
            UA=$(grep "aligned exactly 1 time" $SRX.bt2.log | awk '{print $1}')
            PA=$(grep "overall alignment rate" $SRX.bt2.log | awk '{print $1}')
            FA=$(samtools idxstats $SRX.filter.bam | awk '{S+=$3}END{print S}')
            LA=$(samtools idxstats $SRX.final.bam | awk '{S+=$3}END{print S}')
            echo -e "$SRX\t$TR\t$UA\t$PA\t$FA\t$LA"
        fi 
    done >> $myWorkDIR/meta/align.table.txt 

# **************** 10. Quality Metrics **************** 
# save quality metrics under qc2 directory  
mkdir -p $myWorkDIR/qc2 && cd $myWorkDIR/qc2
# calculate various quality metrics
cat $myWorkDIR/meta/SRX.list | while read SRX; do 
        # PBC metrics: TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
        [[ ! -f $SRX.pbc.qc ]] && \
        echo "bamToBed -i $myWorkDIR/align/$SRX.filter.bam \
            | awk 'BEGIN{OFS=\"\t\"} {print \$1,\$2,\$3,\$6}' \
            | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} \
            (\$1==1){m1=m1+1}(\$1==2){m2=m2+1}{m0=m0+1}{mt=mt+\$1} \
            END{printf \"%d\t%d\t%d\t%d\t%f\t%f\t%f\n\",\
            mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > $SRX.pbc.qc "
        # Cross-correlation analysis
        CMD="Rscript $myWorkDIR/bin/phantompeakqualtools/run_spp.R"
        [[ ! -f $SRX.cc.qc ]] && \
            echo "$CMD -c=$myWorkDIR/align/$SRX.filter.bam \
            -savp=$SRX.cc.plot.pdf -out=$SRX.cc.qc -odir=$myWorkDIR/qc2 \
            -tmpdir=$myWorkDIR/qc2/ > $myWorkDIR/qc2/$SRX.cc.log 2>&1 "
    done | parallelRun 
# SPOT: 
# IDR: https://sites.google.com/site/anshulkundaje/projects/idr

# extract quality metrics  
# QualityTag: Quality tag based on thresholded RSC 
#   (codes: -2:veryLow,-1:Low,0:Medium,1:High,2:veryHigh)
echo -e "Sample\t#Reads\tFragLen\tNSC\tRSC\tQualityTag\tNRF\tPBC1\tPBC2" \
    > $myWorkDIR/meta/qc.table.txt 
cat $myWorkDIR/meta/SRX.list | while read SRX; do 
        if [[ -f $SRX.pbc.qc && -f $SRX.cc.qc ]]; then 
            CC=$(awk -vFS="\t" '{split($3,L,",");FL=L[1]>0?L[1]:(L[2]>0?L[2]:L[3]); \
            printf "%d\t%d\t%f\t%f\t%d",$2,FL,$9,$10,$11}' $SRX.cc.qc)
            PBC=$(awk '{printf "%f\t%f\t%f",$5,$6,$7}' $SRX.pbc.qc)
            echo -e "$SRX\t$CC\t$PBC"
        fi 
    done >> $myWorkDIR/meta/qc.table.txt 

# **************** 11. MACS2 peak and signal **************** 
# prepare a metadata file for MACS2 
awk -vFS="\t" -vOFS="\t" '{a[$2]=$3;if(!($2~/_control/)){b[$2]=$1;}} \
    END{for(i in b){print b[i],i,a[i],a[i"_control"]}}' \
    $myWorkDIR/meta/chip.list | sort -k 1 -k 2 \
    > $myWorkDIR/meta/macs2.meta 

#=====================================
# 1. Generate narrow peaks (FDR<0.001)
# 2. Fold enrichment signal tracks
# 3. -log10(p-value) signal tracks
#=====================================
mkdir -p $myWorkDIR/macs2/signal && cd $myWorkDIR/macs2
# peak calling on each chip-control pairs 
cat $myWorkDIR/meta/macs2.meta | while read TF Rep ChIP Ctrl; do 
        ChIP=$myWorkDIR/align/$ChIP.final.bam 
        Ctrl=$myWorkDIR/align/$Ctrl.final.bam 
        OUT=${Rep}_sc 
        [[ ! -f chrom.sizes ]] && samtools view -H $ChIP | \
        awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); \
            match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' \
            > chrom.sizes 
        GS=$(awk '{GS+=$2}END{print GS}' chrom.sizes)
        [[ ! -f signal/${OUT}.peaks.bb ]] && \
        echo "macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} \
            -g $GS -p 1e-2 --mfold 2 20; \
            maxS=\$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak | head -n 1 | cut -f 5); \
            minS=\$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak | tail -n 1 | cut -f 5); \
            awk -vOFS='\t' -vm=\$minS -vM=\$maxS '{\$5=int(((\$5-m)*(1000-10)/(M-m))+10); print}' ${OUT}_peaks.narrowPeak > ${OUT}.peaks.bed; \
            bedToBigBed -type=bed6+4 ${OUT}.peaks.bed chrom.sizes signal/${OUT}.peaks.bb; \
            rm ${OUT}_model.r ${OUT}_peaks.narrowPeak ${OUT}_summits.bed " 
    done | parallelRun 
# peak calling on merged data
# 1. use merged control data as control 
mkdir -p $myWorkDIR/macs2/signal && cd $myWorkDIR/macs2
cat $myWorkDIR/meta/macs2.meta | while read TF Rep ChIP Ctrl; do 
        ChIP=$myWorkDIR/align/$ChIP.final.bam 
        Ctrl=$myWorkDIR/align/${TF}_control.bam 
        OUT=$Rep 
        [[ ! -f chrom.sizes ]] && samtools view -H $ChIP | \
        awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); \
            match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' \
            > chrom.sizes 
        GS=$(awk '{GS+=$2}END{print GS}' chrom.sizes)
        [[ ! -f signal/${OUT}.peaks.bb ]] && \
        echo "macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} \
            -g $GS -p 1e-2 --mfold 2 20; \
            maxS=\$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak | head -n 1 | cut -f 5); \
            minS=\$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak | tail -n 1 | cut -f 5); \
            awk -vOFS='\t' -vm=\$minS -vM=\$maxS '{\$5=int(((\$5-m)*(1000-10)/(M-m))+10); print}' ${OUT}_peaks.narrowPeak > ${OUT}.peaks.bed; \
            bedToBigBed -type=bed6+4 ${OUT}.peaks.bed chrom.sizes signal/${OUT}.peaks.bb; \
            rm ${OUT}_model.r ${OUT}_peaks.narrowPeak ${OUT}_summits.bed " 
    done | parallelRun 
# 2. use both merged chip and merged control data 
mkdir -p $myWorkDIR/macs2/signal && cd $myWorkDIR/macs2
cut -f 1 $myWorkDIR/meta/macs2.meta | uniq | while read TF; do 
        ChIP=$myWorkDIR/align/${TF}_chip.bam 
        Ctrl=$myWorkDIR/align/${TF}_control.bam 
        OUT=${TF}_rep0 
        [[ ! -f chrom.sizes ]] && samtools view -H $ChIP | \
        awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); \
            match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' \
            > chrom.sizes 
        cSF=$(samtools idxstats $ChIP | awk '!(/*/){TL+=$3} END{print TL/1000000}')
        tSF=$(samtools idxstats $Ctrl | awk '!(/*/){TL+=$3} END{print TL/1000000}')
        SF=$(echo "$cSF $tSF" | awk '($1>$2){print $2} ($1<=$2){print $1}')
        GS=$(awk '{GS+=$2}END{print GS}' chrom.sizes)
        [[ ! -f signal/${OUT}.peaks.bb ]] && \
        echo "macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} -g $GS -p 1e-2 --mfold 2 20 -B --SPMR; \
            maxS=\$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak | head -n 1 | cut -f 5); \
            minS=\$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak | tail -n 1 | cut -f 5); \
            awk -vOFS='\t' -vm=\$minS -vM=\$maxS '{\$5=int(((\$5-m)*(1000-10)/(M-m))+10); print}' ${OUT}_peaks.narrowPeak > ${OUT}.peaks.bed; \
            bedToBigBed -type=bed6+4 ${OUT}.peaks.bed chrom.sizes signal/${OUT}.peaks.bb; \
            macs2 bdgcmp -m FE -t ${OUT}_treat_pileup.bdg -c ${OUT}_control_lambda.bdg -o ${OUT}_FE.bdg; \
            sortBed -i ${OUT}_FE.bdg | slopBed -i - -g chrom.sizes -b 0 | bedClip stdin chrom.sizes ${OUT}.fc.signal.bdg; \
            bedGraphToBigWig ${OUT}.fc.signal.bdg chrom.sizes signal/${OUT}.fc.signal.bw; \
            macs2 bdgcmp -m ppois -t ${OUT}_treat_pileup.bdg -c ${OUT}_control_lambda.bdg -o ${OUT}_ppois.bdg -S $SF; \
            sortBed -i ${OUT}_ppois.bdg | slopBed -i - -g chrom.sizes -b 0 | bedClip stdin chrom.sizes ${OUT}.pval.signal.bdg; \
            bedGraphToBigWig ${OUT}.pval.signal.bdg chrom.sizes signal/${OUT}.pval.signal.bw; \
            rm ${OUT}*.bdg ${OUT}_model.r ${OUT}_summits.bed " 
    done | parallelRun 

# **************** 12. IDR **************** 
# prepare a metadata file for IDR  
mkdir -p $myWorkDIR/macs2/idr && cd $myWorkDIR/bin/idrCode/ 
pDIR=$myWorkDIR/macs2
awk -vFS="\t" -vOFS="\t" '{a[$1]=a[$1]?a[$1]"\t"$2:$2;} \
    END{for(i in a){print i,a[i]}}' \
    $myWorkDIR/meta/macs2.meta | sort -k 1 -k 2 \
    | tee $myWorkDIR/meta/idr.meta | awk -vFS="\t" -vDIR=$pDIR \
    '(NF>2){for(i=2;i<=NF;i++){for(j=i+1;j<=NF;j++){print "Rscript batch-consistency-analysis.r "DIR"/"$i".peaks.bed "DIR"/"$j".peaks.bed -1 "DIR"/idr/"$1"_"$i"_VS_"$j" 0 F p.value"}}}' \
    | parallelRun 
# examples: 
## Rscript batch-consistency-analysis.r $myWorkDIR/macs2/SVP_rep1.peaks.bed $myWorkDIR/macs2/SVP_rep2.peaks.bed -1 $myWorkDIR/macs2/idr/SVP_rep1_VS_SVP_rep2 0 F p.value 
## Rscript batch-consistency-analysis.r $myWorkDIR/macs2/SVP_rep1.peaks.bed $myWorkDIR/macs2/SVP_rep3.peaks.bed -1 $myWorkDIR/macs2/idr/SVP_rep1_VS_SVP_rep3 0 F p.value 
## Rscript batch-consistency-analysis.r $myWorkDIR/macs2/SVP_rep2.peaks.bed $myWorkDIR/macs2/SVP_rep3.peaks.bed -1 $myWorkDIR/macs2/idr/SVP_rep2_VS_SVP_rep3 0 F p.value 
## Rscript batch-consistency-plot.r 1 $myWorkDIR/macs2/idr/SVP_repAll $myWorkDIR/macs2/idr/SVP_rep1_VS_SVP_rep3
# 
# Get thresholds to truncate peak lists and final set of peak calls 
mkdir -p $myWorkDIR/macs2/final 
while read TF Rep MORE; do 
    if [[ ! -z $MORE ]]; then
        echo "*** $TF with replicates"
        numPeaks=
        for cmp in $myWorkDIR/macs2/idr/${TF}_*-overlapped-peaks.txt; do 
            numPeaks="$numPeaks "$( awk '$11 <= 0.1 {print $0}' $cmp | wc -l )
        done 
        max_numPeak=$(echo $numPeaks | awk '{for(i=1;i<=NF;i++){if(!m || $i>m){m=$i}}}END{print m}')
        echo "    numPeaks: $numPeaks [$max_numPeak]"
        ## peaks with FDR < 0.001 and peaks passing a threshold of IDR < 10%
        awk '($9>3)' $myWorkDIR/macs2/${TF}_rep0.peaks.bed \
            | tee $myWorkDIR/macs2/final/$TF.all_peaks.bed \
            | sort -k8nr,8nr | head -n $max_numPeak \
            | bedtools sort -i > $myWorkDIR/macs2/final/$TF.peaks.bed 
    else
        echo "*** $TF without replicates"
        ## peaks with FDR < 0.001 and top 1000 peaks
        awk '($9>3)' $myWorkDIR/macs2/${TF}_rep0.peaks.bed \
        | tee $myWorkDIR/macs2/final/$TF.all_peaks.bed \
        | sort -k8nr,8nr | head -n 1000 \
        | bedtools sort -i > $myWorkDIR/macs2/final/$TF.peaks.bed 
        echo "    total peaks: "$(cat $myWorkDIR/macs2/final/$TF.all_peaks.bed | wc -l)" ["$(cat $myWorkDIR/macs2/final/$TF.peaks.bed | wc -l)"]"
    fi
done < $myWorkDIR/meta/idr.meta 

# **************** 13. Peak annotation **************** 
# peak annotation using ChIPseeker: 
# http://bioconductor.org/packages/release/bioc/html/ChIPseeker.html  
cut -f 1 $myWorkDIR/meta/idr.meta \
    | awk '{a=$1":"$1".peaks.bed";s=s?s"\t"a:a}END{print "FlowerTF\t"s}' \
    | peakAnnotation -i $myWorkDIR/macs2/final -o $myWorkDIR/annotation \
    -g $myWorkDIR/db/tair10Genes.gff3 
# compare peaks from the original publication

cd $myWorkDIR/macs2/final 
echo `pwd` && rm -f SVP.peaks.bed 
for peak in *.peaks.bed; do 
    TF=$(echo $peak | sed -e 's/.peaks.bed//')
    A=$(cat $peak | wc -l)
    B=$(cat /home/data/pub/studies/2016_Yan_COPB/v3/top1000/$TF.bed | wc -l)
    C=$(bedtools intersect -u -a $peak -b /home/data/pub/studies/2016_Yan_COPB/v3/top1000/$TF.bed | wc -l)
    D=$(bedtools intersect -u -b $peak -a /home/data/pub/studies/2016_Yan_COPB/v3/top1000/$TF.bed | wc -l)
    echo -e "$TF\t$A\t$C\t$B\t$D"
done > peak_overlap_with_publication.txt 

# **************** 14. Target genes **************** 
# Genomic regions were associated with genes if located 3 kb upstream of the start of the gene up to 1 kb downstream of the end of the gene
mkdir -p $myWorkDIR/target && cd $myWorkDIR/target 
awk -vFS="\t" -vOFS="\t" '{if($6=="+"){$2=$2-3000;$3=$3+1000;if($2<0){$2=0;}}else{$2=$2-1000;$3=$3+3000;if($2<0){$2=0;}} print}' $myWorkDIR/db/tair10Genes.bed > gene.3u1d.bed 
awk -vFS="\t" -vOFS="\t" '{if($6=="+"){$2=$2-1500;$3=$3+500;if($2<0){$2=0;}}else{$2=$2-500;$3=$3+1500;if($2<0){$2=0;}} print}' $myWorkDIR/db/tair10miRNAs.bed >> gene.3u1d.bed 

cut -f 1 $myWorkDIR/meta/idr.meta | while read TF; do 
    bedtools intersect -u -a gene.3u1d.bed -b $myWorkDIR/macs2/final/$TF.peaks.bed \
    | cut -f4,5 | awk -vTF=$TF '{print TF"\t"$0}' > $TF.target.list 
    bedtools intersect -wo -a gene.3u1d.bed -b $myWorkDIR/macs2/final/$TF.peaks.bed \
    | cut -f4,5,7-11 | awk -vTF=$TF '{print TF"\t"$0}' > $TF.target.score 
done 
rm -f SVP.target.list 
cat *.target.list | awk -vFS="\t" '{id=$2"\t"$3;n[id]++;a[id]=a[id]?a[id]","$1:$1}END{for(i in a) print i"\t"n[i]"\t"a[i]}' \
    | tee target.table.txt | awk -vFS="\t" 'BEGIN{while(getline<"TF.txt"){a[$1]++}}(a[$1])' > TFtarget.table.txt
awk -vMIR=$myWorkDIR/db/tair10miRNAs.bed 'BEGIN{while(getline<MIR){a[$4]++}}(a[$1])' target.table.txt > MIRtarget.table.txt
cat *.target.list | tee target.list | awk -vFS="\t" 'BEGIN{while(getline<"TF.txt"){a[$1]++}}(a[$2])' > TFtarget.list
cat *.target.list | tee target.list | awk -vFS="\t" -vMIR=$myWorkDIR/db/tair10miRNAs.bed 'BEGIN{while(getline<MIR){a[$4]++}}(a[$2])' > MIRtarget.list

# **************** 15. Summary of target genes **************** 
mkdir -p $myWorkDIR/target/stat && cd $myWorkDIR/target/stat
for list in $myWorkDIR/target/*.target.list; do 
    TF=$(basename $list | sed s/.target.list//)
    cut -f 2 $list > $TF.gene 
    cat $myWorkDIR/macs2/final/$TF.peaks.bed > $TF.peak
done 

file_labels=`ls *.gene | sed -e 's/.gene//g'`
echo name" "$file_labels > gene_pairwise_overlapping_rate.txt
echo name" "$file_labels > gene_pairwise_overlapping_number.txt
for ds1 in *.gene; do
    p1=$(echo $ds1 | sed s/.gene//)
    n1=$(cat $ds1 | wc -l)
    echo -n $p1","$n1 >> gene_pairwise_overlapping_rate.txt
    echo -n $p1","$n1 >> gene_pairwise_overlapping_number.txt
    for ds2 in *.gene; do
        p2=$(echo $ds2 | sed s/.gene//)
        n2=$(cat $ds2 | wc -l)
        a=$(sort $ds1 $ds2 | uniq -d | wc -l)
        b=$(cat $ds1 | wc -l)
        c=$(cat $ds2 | wc -l)
        ## rate=$(bc -l <<< "$a/($b/2+$c/2)")
        rate=$(bc -l <<< "$a/$b")
        echo -n " "$rate >> gene_pairwise_overlapping_rate.txt
        echo -n " "$a >> gene_pairwise_overlapping_number.txt
    done
    echo >> gene_pairwise_overlapping_rate.txt
    echo >> gene_pairwise_overlapping_number.txt
done


file_labels=`ls *.peak | sed -e 's/.peak//g'`
echo name" "$file_labels > peak_pairwise_overlapping_rate.txt
echo name" "$file_labels > peak_pairwise_overlapping_number.txt
for ds1 in *.peak; do
    p1=$(basename $ds1 | sed s/.peak//)
    n1=$(cat $ds1 | wc -l)
    echo -n $p1","$n1 >> peak_pairwise_overlapping_rate.txt
    echo -n $p1","$n1 >> peak_pairwise_overlapping_number.txt
    for ds2 in *.peak; do
        a=$(intersectBed -a $ds1 -b $ds2 -u | wc -l )
        b=$(cat $ds1 | wc -l)
        c=$(cat $ds2 | wc -l)
        ## rate=$(bc -l <<< "$a/($b/2+$c/2)")
        rate=$(bc -l <<< "$a/$b")
        echo -n " "$rate >> peak_pairwise_overlapping_rate.txt
        echo -n " "$a >> peak_pairwise_overlapping_number.txt
    done
    echo >> peak_pairwise_overlapping_rate.txt
    echo >> peak_pairwise_overlapping_number.txt
done

Rscript corrplot.R


# R code 
R 
pdf("Fig. 3c.pdf", width=8.27, heigh=8.27, pointsize=10)
op <- par(mfrow=c(1, 2))
stat <- read.table("peak_overlap_with_publication.txt", row.names=1, head=FALSE, sep="\t")
bplt1 <- t(cbind(stat[,2]/stat[,1]*100, 100-stat[,2]/stat[,1]*100))
bplt2 <- t(cbind(stat[,4]/stat[,3]*100, 100-stat[,4]/stat[,3]*100))
colnames(bplt1) <- colnames(bplt2) <- rownames(stat)
bp1 <- barplot(bplt1, horiz=T, border=NA, col=c("#fdae61", "#878787"), xlab="Percentage of overlapping", main="Re-analysis")
text(100, bp1, labels=stat[,1], pos=4, xpd=T)
bp2 <- barplot(bplt2, horiz=T, border=NA, col=c("#fdae61", "#878787"), xlab="Percentage of overlapping", main="Publication")
text(100, bp2, labels=stat[,3], pos=4, xpd=T)
par(op)
dev.off()

#### summary 
options(stringsAsFactors=FALSE)
TFstarstat <- read.table("TFtarget.table.txt", head=FALSE, sep="\t")
MIRtarstat <- read.table("MIRtarget.table.txt", head=FALSE, sep="\t")
Alltarstat <- read.table("target.table.txt", head=FALSE, sep="\t")
OTHtarstat <- Alltarstat[!(Alltarstat[,1] %in% c(TFstarstat[,1], MIRtarstat[,1])), ]
TFstarstat[TFstarstat[, 2]=="", 2] <- TFstarstat[TFstarstat[, 2]=="", 1]
rownames(TFstarstat) <- TFstarstat[, 2]
rownames(MIRtarstat) <- MIRtarstat[, 2]
rownames(Alltarstat) <- Alltarstat[, 1]

pdf("Fig. 4.pdf", width=8.27, heigh=8.27, pointsize=10)
op <- par(mfrow=c(2, 2))
TFstat <- table(TFstarstat[,3])
NN <- names(TFstat)
MIRstat <- table(MIRtarstat[,3])[NN]
OTHstat <- table(OTHtarstat[,3])[NN]
MIRstat[is.na(MIRstat)] <- 0
OTHstat[is.na(OTHstat)] <- 0

bardata <- cbind("TF targets"=TFstat, "MIR targets"=MIRstat, "Others"=OTHstat)
plotdata <- t(bardata) / colSums(bardata)
barplot(plotdata, ylab="% targets", xlab="# regulators", beside=T, border=NA, col=c('#377eb8', '#ff7f00', '#999999'), args.legend=list(box.col=NA, bg=NA, ncol=1, x="topright", border=NA), legend.text=paste0(rownames(plotdata), " (", c(sum(TFstat), sum(MIRstat), sum(OTHstat)), ")"))
par(op)
dev.off()

#### network  
library(igraph)
library(RColorBrewer)
## display.brewer.all()

topTFs <- TFstarstat[with(TFstarstat, order(-V3)), ][1:10, 2]
topMIR <- MIRtarstat[with(MIRtarstat, order(-V3)), ][1:10, 2]
TFs <- TFstarstat[,2]
MIRs <- MIRtarstat[,2]

TFtargets <- read.table("TFtarget.list", head=FALSE, sep="\t")
TFtargets[TFtargets[, 3]=="", 3] <- TFtargets[TFtargets[, 3]=="", 2]
MIRtargets <- read.table("MIRtarget.list", head=FALSE, sep="\t")

net.data <- rbind(TFtargets[, c(1,3)], MIRtargets[, c(1,3)])
regulators <- unique(c(TFtargets[,1], MIRtargets[,1]))
vertices <- unique(c(TFs, MIRs, regulators))
net <- graph.data.frame(net.data, directed=T, vertices=vertices)
node.color <- rep(adjustcolor("grey80",alpha=.4), vcount(net))
node.size <- rep(3, vcount(net))
node.label <- rep(NA, vcount(net))
node.shape <- rep("circle", vcount(net))
names(node.color) <- names(node.size) <- names(node.label) <- names(node.shape) <- vertices

## node.color[c(topTFs, topMIR)] <- "#ff9d00"
node.degree <- rbind(TFstarstat, MIRtarstat)[c(MIRs, TFs), 3]
node.color[c(MIRs, TFs)] <- colorRampPalette(brewer.pal(n=7,name="Greys"))(max(node.degree))[node.degree]
node.color[regulators] <- "#ed6d45" 
node.label[c(regulators, topTFs, topMIR)] <- c(regulators, topTFs, topMIR)
node.size[c(MIRs, TFs)] <- 2+node.degree
node.size[regulators] <- 12
## node.shape[MIRs] <- "square"
frame.color <- node.color 
frame.color[] <- 'grey50'
frame.color[MIRs] <- '#3aa935'
frame.color[regulators] <- "#be1621"
label.cex <- node.size
label.cex[] <- .5
label.cex[regulators] <- 1

V(net)$color <- node.color
V(net)$frame.color <- frame.color
V(net)$label <- node.label
V(net)$size <- node.size
V(net)$shape <- node.shape
V(net)$label.cex <- label.cex

E(net)$color <- ifelse(net.data[,2] %in% regulators, adjustcolor("#f46d43",alpha=.5), adjustcolor("grey80",alpha=.5))
## E(net)$color <- adjustcolor("grey80",alpha=.5)
E(net)$width <- 2 
E(net)$arrow.size <- .4
E(net)$arrow.width <- .4

g1 <- c("AP1", "AG", "SEP3", "PI", "AP3")
g2 <- c("FLC", "FLM", "SVP", "SOC1")
mkgs <- list(c(which(vertices %in% g1)), c(which(vertices %in% g2)))

pdf("Fig. 5.pdf", width=8.27, heigh=8.27, pointsize=10)
set.seed(1234)
plot(net, layout=layout.auto, mark.groups=mkgs, mark.col=c("#36a9e0","#f9b233"), mark.border=NA, edge.arrow.size=.8, vertex.label.font=3)
set.seed(1234)
plot(net, layout=layout.kamada.kawai, mark.groups=mkgs, mark.col=c("#36a9e0","#f9b233"), mark.border=NA, edge.arrow.size=.8, vertex.label.font=3)
set.seed(1234)
plot(net, layout=layout.fruchterman.reingold, mark.groups=mkgs, mark.col=c("#36a9e0","#f9b233"), mark.border=NA, edge.arrow.size=.8, vertex.label.font=3)
set.seed(1234)
plot(net, layout=layout.fruchterman.reingold, mark.groups=mkgs, mark.col=c("#36a9e0","#f9b233"), mark.border=NA, edge.arrow.size=.8, vertex.label.font=3)
set.seed(1234)
plot(net, layout=layout.drl, mark.groups=mkgs, mark.col=c("#36a9e0","#f9b233"), mark.border=NA, edge.arrow.size=.8, vertex.label.font=3)
## plot(net, layout=layout.svd, mark.groups=mkgs, mark.col=c("#36a9e0","#f9b233"), mark.border=NA, edge.arrow.size=.8, vertex.label.font=3)
## plot(net, layout=layout.spring, mark.groups=mkgs, mark.col=c("#36a9e0","#f9b233"), mark.border=NA, edge.arrow.size=.8, vertex.label.font=3)
dev.off()

pdf("Fig. 5-1.pdf", width=8.27, heigh=8.27, pointsize=10)
layouts <- grep("^layout\\.", ls("package:igraph"), value=TRUE) 
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama", layouts)]

## op <- par(mfrow=c(2,2))
for (layout in layouts) {
    print(layout)
    l <- do.call(layout, list(net)) 
    set.seed(1234)
    plot(net, edge.arrow.mode=0, layout=l, main=layout)
}
## par(op)
dev.off()

#### expression 
library(preprocessCore)
library(pheatmap)
library(MASS)

exprs <- read.table("expression.table.txt", row.names=1, head=TRUE, sep="\t", check.names=FALSE)
annot <- read.table("expression.table.header.txt", head=TRUE, sep="\t", check.names=FALSE)
dim(exprs)
dim(annot)
rid <- grep("\\|", rownames(exprs))
if(length(rid) > 0){
    srows <- rownames(exprs)[rid]
    drows <- sapply(srows, strsplit, split="\\|")
    nrows <- irows <- c()
    for(r in names(drows)){
        nrows <- c(nrows, drows[[r]])
        irows <- c(irows, rep(r, length(drows[[r]])))
    }
    arows <- exprs[irows, ]
    rownames(arows) <- nrows
    exprs <- exprs[-rid, ]
    exprs <- rbind(exprs, arows)
}

cnames <- colnames(exprs)
rnames <- rownames(exprs)
## exprs <- normalize.quantiles(as.matrix(exprs)) 
## colnames(exprs) <- cnames
## rownames(exprs) <- rnames

overlapgenes <- intersect(rownames(Alltarstat), rownames(exprs))

# columns <- intersect(which(annot$control=="No"), grep("flower|inflorescence|Shoot Apex", annot$tissue, ignore.case=TRUE))
columns <- intersect(which(annot$control=="No"), grep("wt", annot$mutant, ignore.case=TRUE))
## columns <- intersect(columns, grep("flower|inflorescence|Shoot Apex", annot$tissue, ignore.case=TRUE))
columns <- intersect(columns, grep("Developmental", annot$category, ignore.case=TRUE))
expdata <- exprs[overlapgenes, columns]
colnames(expdata) <- paste0(annot[columns, 'age'], ": ", annot[columns, 'tissue']) # , "; ", annot[columns, 'mutant']
logdata <- log10(expdata)
CV <- function(x){
    mean <- mean(x, na.rm=TRUE)
    sd <- sd(x, na.rm=TRUE)
    (sd/mean)*100
}
cvdata <- sort(abs(apply(logdata, 1, CV)), decreasing=T)
tops <- c(names(head(cvdata, n=100)), names(tail(cvdata, n=100)))
tops <- overlapgenes
logdata <- logdata[tops, ]

qtl <- quantile(as.matrix(logdata), probs=c(0.025,0.975))
logdata[logdata < qtl[1]] <- qtl[1]
logdata[logdata > qtl[2]] <- qtl[2]

target2TFs <- Alltarstat[tops, 4] 
target2TFs <- sapply(target2TFs, strsplit, split=",") 
names(target2TFs) <- tops
annotation_col <- data.frame(Tissue=rep("Other", ncol(logdata)), row.names=colnames(logdata))
annotation_col[grep("Inflorescence|flower|Transition", rownames(annotation_col), ignore.case=TRUE), "Tissue"] <- "Floral"
annotation_row <- as.data.frame(matrix("N", ncol=length(regulators), nrow=length(tops)))
rownames(annotation_row) <- tops
colnames(annotation_row) <- regulators
for(gene in tops){
    annotation_row[gene, target2TFs[[gene]]] <- "Y"
}
ann_colors <- list(Tissue=c(Floral="#c51b7d", Other="#1b7837"))
for(regulator in regulators){
    ann_colors[[regulator]] <- c(Y="#1b7837", N="#ffffff")
}


pdf("Fig. 6.pdf", width=8.27, heigh=8.27, pointsize=10)
phtm <- pheatmap(logdata, annotation_row=annotation_row, annotation_col=annotation_col, annotation_colors=ann_colors, fontsize=6, border_color=NA, show_rownames=F, useRaster=T) ## ,cutree_rows=8,  cluster_cols=F
geneord <- phtm$tree_row$labels[phtm$tree_row$order]
TFBSs <- Alltarstat[geneord, 3]
plot(TFBSs, type="h", col=colorRampPalette(brewer.pal(n=7,name="Reds"))(max(TFBSs))[TFBSs])
dev.off()


exit()

############
WTcid <- intersect(grep("Shoot Apex, Inflorescence", annot$tissue, ignore.case=TRUE), grep("wt", annot$mutant, ignore.case=TRUE))
lfycid <- intersect(grep("Shoot Apex, Inflorescence", annot$tissue, ignore.case=TRUE), grep("lfy-12", annot$mutant, ignore.case=TRUE))
ap1cid <- intersect(grep("Shoot Apex, Inflorescence", annot$tissue, ignore.case=TRUE), grep("ap1-15", annot$mutant, ignore.case=TRUE))
ap2cid <- intersect(grep("Shoot Apex, Inflorescence", annot$tissue, ignore.case=TRUE), grep("ap2-6", annot$mutant, ignore.case=TRUE))
ap3cid <- intersect(grep("Shoot Apex, Inflorescence", annot$tissue, ignore.case=TRUE), grep("ap3-6", annot$mutant, ignore.case=TRUE))
agcid <- intersect(grep("Shoot Apex, Inflorescence", annot$tissue, ignore.case=TRUE), grep("ag-12", annot$mutant, ignore.case=TRUE))

AP1score <- read.table("AP1.target.score", head=FALSE, sep="\t", check.names=FALSE)
AP1score <- AP1score[AP1score[,2] %in% overlapgenes, ]
px <- as.numeric(AP1score[, 8])
py <- as.numeric(exprs_norm[AP1score[,2], WTcid]-exprs_norm[AP1score[,2], ap1cid])
f1 <- kde2d(log2(px), log2(py), n=400)
plot(log2(px), log2(py))
abline(lm(y ~ x, data=data.frame(x=log2(px), y=log2(py))))
dev.off()
