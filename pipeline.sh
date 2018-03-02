## Feel free to contact me via chendijun2012@gmail.com for bugs, comments and suggestions. 

## Note: be aware that this pipeline only works for sing-end read files 
##       if you are working on paired-end data, some parts of code should 
##       have minor changes (e.g, using different parameters)

# **************** 1. workspace     ****************
# create a workspace under $HOME/flowerGRNs
mkdir -p $HOME/works/flowerGRNs && myWorkDIR=$HOME/works/flowerGRNs

# **************** 2. get source code   ****************
# change directory to workspace
cd $myWorkDIR 
# clone the flowerGRNs repo from https://github.com/PlantENCODE/
rm -rf pipeline && git clone --recursive https://github.com/PlantENCODE/plantGRNs.git pipeline
# make all scripts under pipeline directory executable 
chmod +x pipeline/* 
# put metadata files under meta folder
mv $myWorkDIR/pipeline/meta $myWorkDIR/ 

# **************** 3. installation  **************** 
# create a directory bin store executable programs 
mkdir -p $myWorkDIR/bin 
# install third-party software packages into this directory 
# NOTE: follow the INSTALL.sh for instruction on how to install third-party software packages 
cd $myWorkDIR/bin && sh $myWorkDIR/pipeline/INSTALL.sh
# make executable programs available through the PATH environmental variable
export PATH=$myWorkDIR/bin:$myWorkDIR/pipeline:$PATH

# **************** 4. ChIP-seq data **************** 
# create a sub-directory fq to save ChIP-seq FASTQ files 
mkdir -p $myWorkDIR/fq && cd $myWorkDIR/fq 

# get experiment metadata file 
perl -MLWP::Simple -e "getprint 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=SRP022770&result=read_run&fields=experiment_accession,run_accession,fastq_ftp,sample_alias,library_strategy&download=txt'" | awk '(NR>1)' > chipseq.meta  
# download FASTQ files in parallel (maximum 10 runs)
cut -f 3 chipseq.meta | xargs -P 10 \
    -I % sh -c 'O=$(basename %); wget -O $O -q %'
# merge multiple SRR files (if any) into a single SRX file 
awk -vOFS="\t" '{a[$1]=a[$1]" "$2".fastq.gz"} \
    END{for(x in a){print x,a[x]}}' chipseq.meta \
    | while read SRX SRRs; do 
        cat $SRRs > $SRX.fastq.gz
    done 
rm -f SRR*.fastq.gz 

# **************** 5. genome data   **************** 
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
awk -vFS="\t" -vOFS="\t" 'BEGIN{
    while(getline<"gene_aliases.txt"){alias[$1]=$2}} \
    ($3=="gene" && $9~/Note=protein_coding_gene/) \
        {match($9,"ID=([^;$]+)",ID); \
        print $1,$4,$5,ID[1],alias[ID[1]],$7}' \
    tair10Genes.gff3 | sed -e "s/'//g" \
    | sort -k1,1 -k2,2n > tair10Genes.bed 

# **************** 6. FastQC    **************** 
cut -f 3 $myWorkDIR/meta/chip-seq.txt | sort | uniq > $myWorkDIR/meta/SRX.list 
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
# using command: fastq_quality_filter or fastx_trimmer
# no need here 

# **************** 8. Alignment **************** 
# save alignment result under align directory  
mkdir -p $myWorkDIR/align && cd $myWorkDIR/align
# FASTQ encoding: Sanger/Illumina_1_8 format (--phred33); 
# use 2 threads for alignment (--threads 2); 
# suppress SAM records for unaligned reads (--no-unal)
NTHREADS=2 ## number of threads
cat $myWorkDIR/meta/SRX.list | while read SRX; do 
        FQ="$myWorkDIR/fq/$SRX.fastq.gz"
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
miniQ=1
PICARD=$myWorkDIR/bin/picard ## path for Picard tools
cd $myWorkDIR/align
cat $myWorkDIR/meta/SRX.list | while read SRX; do 
        # Filter bam file, based on FLAG 1804: segment unmapped (4) + next segment in the template unmapped (8) + secondary alignments (256) + not passing filters, such as platform/vendor quality controls (512) + PCR or optical duplicates (1024); -q 30 exclude MAPQ < 30 
        [[ ! -f $SRX.filter.bam || ! -f $SRX.final.bam ]] && \
        echo "samtools view -@ ${NTHREADS} -F 1804 -q $miniQ -b $SRX.bam \
            > $SRX.filter.bam && samtools index $SRX.filter.bam; \
            java -jar $PICARD/picard.jar MarkDuplicates \
            I=$SRX.filter.bam O=$SRX.dupmark.bam M=$SRX.dup.qc \
            VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false \
            ASSUME_SORTED=true; samtools view -@ ${NTHREADS} \
            -F 1804 -b $SRX.dupmark.bam > $SRX.final.bam; \
            samtools index $SRX.final.bam && rm $SRX.dupmark.bam "
    done | parallelRun

# alignment summary 
cd $myWorkDIR/align
echo -e "TF\tRep\tSample\t#Raw\t#Align\t%Align\t#Filtered\tFinal" \
    > $myWorkDIR/meta/align.table.txt 
cat $myWorkDIR/meta/chip-seq.txt | while read TF Rep SRX; do 
    [[ $Rep =~ '_control' ]] && continue; 
    if [[ -f $SRX.filter.bam ]]; then 
        TR=$(grep "reads; of these" $SRX.bt2.log | awk '{print $1}')
        UA=$(grep "aligned exactly 1 time" $SRX.bt2.log | awk '{print $1}')
        PA=$(grep "overall alignment rate" $SRX.bt2.log | awk '{print $1}')
        FA=$(samtools idxstats $SRX.filter.bam | awk '{S+=$3}END{print S}')
        LA=$(samtools idxstats $SRX.final.bam | awk '{S+=$3}END{print S}')
        echo -e "$TF\t$Rep\t$SRX\t$TR\t$UA\t$PA\t$FA\t$LA"
    fi 
    done >> $myWorkDIR/meta/align.table.txt 

# randomly split each alignment file (for ChIP sample) into 2 parts (self pseudo-replicates) with approximately equal number of reads 
cd $myWorkDIR/align
awk '(!/_control/){print $2, $3}' $myWorkDIR/meta/chip-seq.txt \
    | while read Rep SRX; do 
        [[ ! -f ${Rep}.pr1.bam || ! -f ${Rep}.pr2.bam ]] && \
            echo "NN=\$(samtools idxstats $SRX.final.bam \
                | awk '{S+=\$3}END{printf(\"%d\n\",S/2+0.9)}'); \
            samtools view -H $SRX.final.bam > $Rep.header; \
            samtools view $SRX.final.bam | shuf | split -d \
                -l \$NN - ${Rep}; \
            cat $Rep.header ${Rep}00 | samtools view -bS - \
                | samtools sort - ${Rep}.pr1 & \
            cat $Rep.header ${Rep}01 | samtools view -bS - \
                | samtools sort - ${Rep}.pr2; \
            samtools index ${Rep}.pr1.bam & \
            samtools index ${Rep}.pr2.bam; \
            rm -f ${Rep}0* $Rep.header "
    done | parallelRun 

# merge replicates (if any) into a pooled replicate 
# If only one replicate exists for a sample, make symbolic links of the pooled sample to the replicate. 
cd $myWorkDIR/align
awk -vFS="\t" -vOFS="\t" '{bam=$3".final.bam"; \
    if(($2~/_control/)){c=$1"_control"}else{c=$1"_rep0"} \
        if(!f[c][bam]){a[c]=a[c]?a[c]";"bam:bam} f[c][bam]++} \
        END{for(i in a){print i,a[i]}}' \
    $myWorkDIR/meta/chip-seq.txt \
    | while read TF bams; do 
        out=$TF".bam"
        if [[ $bams =~ ";" ]]; then 
            inbams=$(echo $bams | sed -e 's/;/ /g')
            [[ ! -f $out || ! -f "$out.bai" ]] && \
            echo "samtools merge $out $inbams; samtools index $out"
        else 
            ln -sf $bams $out && ln -sf $bams".bai" $out".bai";
        fi
    done | parallelRun 

# randomly split the pooled sample into 2 parts (pooled pseudo-replicates) as above 
cd $myWorkDIR/align
cut -f 1 $myWorkDIR/meta/macs2.meta | uniq | while read TF; do 
        BAM=$myWorkDIR/align/${TF}_rep0.bam 
        [[ ! -f $BAM ]] && continue; 
        OUT="${TF}_rep0"
        [[ ! -f ${OUT}.pr1.bam || ! -f ${OUT}.pr2.bam ]] && \
            echo "NN=\$(samtools idxstats $BAM \
                | awk '{S+=\$3}END{printf(\"%d\n\",S/2+0.9)}'); \
            samtools view -H $BAM > $OUT.header; \
            samtools view $BAM | shuf | split -d -l \$NN - ${OUT}; \
            cat $OUT.header ${OUT}00 | samtools view -bS - \
                | samtools sort - ${OUT}.pr1 & \
            cat $OUT.header ${OUT}01 | samtools view -bS - \
                | samtools sort - ${OUT}.pr2; \
            samtools index ${OUT}.pr1.bam & \
            samtools index ${OUT}.pr2.bam; \
            rm -f ${OUT}0* $OUT.header "
    done | parallelRun 

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

# extract quality metrics  
# QualityTag: Quality tag based on thresholded RSC 
#   (codes: -2:veryLow,-1:Low,0:Medium,1:High,2:veryHigh)
echo -e "TF\tRep\tSample\t#Reads\tFragLen\tNSC\tRSC\tQualityTag\tNRF\tPBC1\tPBC2" \
    > $myWorkDIR/meta/qc.table.txt 
cat $myWorkDIR/meta/chip-seq.txt | while read TF Rep SRX; do 
    [[ $Rep =~ '_control' ]] && continue; 
        if [[ -f $SRX.pbc.qc && -f $SRX.cc.qc ]]; then 
            CC=$(awk -vFS="\t" '{split($3,L,",");FL=L[1]>0?L[1]:(L[2]>0?L[2]:L[3]); \
            printf "%d\t%d\t%f\t%f\t%d",$2,FL,$9,$10,$11}' $SRX.cc.qc)
            PBC=$(awk '{printf "%f\t%f\t%f",$5,$6,$7}' $SRX.pbc.qc)
            echo -e "$TF\t$Rep\t$SRX\t$CC\t$PBC"
        fi 
    done >> $myWorkDIR/meta/qc.table.txt 

# **************** 11. MACS2 peak and signal **************** 
# prepare a metadata file for MACS2 
awk -vFS="\t" -vOFS="\t" '{a[$2]=$3;if(!($2~/_control/)){b[$2]=$1;}} \
    END{for(i in b){print b[i],i,a[i],a[i"_control"]}}' \
    $myWorkDIR/meta/chip-seq.txt | sort -k 1 -k 2 \
    > $myWorkDIR/meta/macs2.meta 

#=====================================
# 1. Generate narrow peaks (FDR<0.001)
# 2. Fold enrichment signal tracks
# 3. -log10(p-value) signal tracks
#=====================================
mkdir -p $myWorkDIR/macs2/signal && cd $myWorkDIR/macs2
# 1. peak calling for each replicate, use pooled control sample data as control 
cat $myWorkDIR/meta/macs2.meta | while read TF Rep Trmt Ctrl; do 
        ChIP=$myWorkDIR/align/$Trmt.final.bam 
        Ctrl=$myWorkDIR/align/${TF}_control.bam 
        OUT=${Rep} 
        [[ ! -f chrom.sizes ]] && samtools view -H $ChIP | \
        awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); \
            match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' \
            > chrom.sizes 
        GS=$(awk '{GS+=$2}END{print int(0.85*GS)}' chrom.sizes)
        [[ ! -f ${OUT}.peaks.bed ]] && \
        echo "macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} \
            -g $GS -p 1e-2 --mfold 2 20 --to-large; \
            maxS=\$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak \
                    | head -n 1 | cut -f 5); \
            minS=\$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak \
                    | tail -n 1 | cut -f 5); \
            awk -vOFS='\t' -vm=\$minS -vM=\$maxS '{ \
                \$5=int(((\$5-m)*(1000-10)/(M-m))+10); print}' \
                ${OUT}_peaks.narrowPeak > ${OUT}.peaks.bed; \
            bedToBigBed -type=bed6+4 ${OUT}.peaks.bed \
                chrom.sizes signal/${OUT}.peaks.bb; \
            rm -f ${OUT}_model.r ${OUT}_peaks.narrowPeak \
                ${OUT}_summits.bed ${OUT}_peaks.xls " 
    done | parallelRun 
# 2. peak calling for pseudo-replicates, using analogous commands are for the individual original replicates (above).
cd $myWorkDIR/macs2
# 2.1 self pseudo-replicates
cat $myWorkDIR/meta/macs2.meta | while read TF Rep Trmt Ctrl; do 
        Ctrl=$myWorkDIR/align/${TF}_control.bam 
        for pRep in pr1 pr2; do 
            ChIP=$myWorkDIR/align/$Rep.$pRep.bam 
            OUT=${Rep}_self_${pRep}
            [[ ! -f chrom.sizes ]] && samtools view -H $ChIP | \
            awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); \
                match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' \
                > chrom.sizes 
            GS=$(awk '{GS+=$2}END{print int(0.85*GS)}' chrom.sizes)
            [[ ! -f ${OUT}.peaks.bed ]] && \
            echo "macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} \
                -g $GS -p 1e-2 --mfold 2 20 --to-large; \
                mv ${OUT}_peaks.narrowPeak ${OUT}.peaks.bed; \
                rm -f ${OUT}_model.r ${OUT}_peaks.narrowPeak \
                    ${OUT}_summits.bed ${OUT}_peaks.xls " 
        done 
    done | parallelRun 
# 2.2 pooled pseudo-replicates
cut -f 1 $myWorkDIR/meta/macs2.meta | uniq | while read TF; do 
        Ctrl=$myWorkDIR/align/${TF}_control.bam 
        for pRep in pr1 pr2; do 
            ChIP=$myWorkDIR/align/${TF}_rep0.${pRep}.bam 
            [[ ! -f $ChIP ]] && continue; 
            OUT=${TF}_rep0_${pRep}
            [[ ! -f chrom.sizes ]] && samtools view -H $ChIP | \
            awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); \
                match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' \
                > chrom.sizes 
            GS=$(awk '{GS+=$2}END{print int(0.85*GS)}' chrom.sizes)
            [[ ! -f ${OUT}.peaks.bed ]] && \
            echo "macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} \
                -g $GS -p 1e-2 --mfold 2 20 --to-large; \
                mv ${OUT}_peaks.narrowPeak ${OUT}.peaks.bed; \
                rm -f ${OUT}_model.r ${OUT}_peaks.narrowPeak \
                    ${OUT}_summits.bed ${OUT}_peaks.xls " 
        done 
    done | parallelRun 
# 3. peak calling based on pooled chip and pooled control data 
mkdir -p $myWorkDIR/macs2/signal && cd $myWorkDIR/macs2
cut -f 1 $myWorkDIR/meta/macs2.meta | uniq | while read TF; do 
        ChIP=$myWorkDIR/align/${TF}_rep0.bam 
        Ctrl=$myWorkDIR/align/${TF}_control.bam 
        [[ ! -f $ChIP || ! -f $Ctrl ]] && continue; 
        OUT=${TF}_rep0 
        [[ ! -f chrom.sizes ]] && samtools view -H $ChIP | \
        awk -vOFS="\t" '(/^@SQ/){match($0,/SN:(\w+)/,SN); \
            match($0,/LN:([0-9]+)/,LN);print SN[1],LN[1]}' \
            > chrom.sizes 
        cSF=$(samtools idxstats $ChIP | awk '!(/*/){TL+=$3} \
            END{print TL/10e+6}')
        tSF=$(samtools idxstats $Ctrl | awk '!(/*/){TL+=$3} \
            END{print TL/10e+6}')
        SF=$(echo "$cSF $tSF" | awk '($1>$2){print $2} \
            ($1<=$2){print $1}')
        GS=$(awk '{GS+=$2}END{print int(0.85*GS)}' chrom.sizes)
        [[ ! -f ${OUT}.peaks.bed ]] && \
        echo "macs2 callpeak -f BAM -t $ChIP -c $Ctrl -n ${OUT} \
            -g $GS -p 1e-2 --mfold 2 20 -B --SPMR --to-large; \
            maxS=\$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak \
                    | head -n 1 | cut -f 5); \
            minS=\$(sort -k 5gr,5gr ${OUT}_peaks.narrowPeak \
                    | tail -n 1 | cut -f 5); \
            awk -vOFS='\t' -vm=\$minS -vM=\$maxS '{ \
                \$5=int(((\$5-m)*(1000-10)/(M-m))+10); print}' \
                ${OUT}_peaks.narrowPeak > ${OUT}.peaks.bed; \
            bedToBigBed -type=bed6+4 ${OUT}.peaks.bed \
                chrom.sizes signal/${OUT}.peaks.bb; \
            macs2 bdgcmp -m FE -t ${OUT}_treat_pileup.bdg \
                -c ${OUT}_control_lambda.bdg -o ${OUT}_FE.bdg; \
            sortBed -i ${OUT}_FE.bdg | slopBed -i - -g chrom.sizes \
                -b 0 | bedClip stdin chrom.sizes ${OUT}.fc.signal.bdg; \
            bedGraphToBigWig ${OUT}.fc.signal.bdg \
                chrom.sizes signal/${OUT}.fc.signal.bw; \
            macs2 bdgcmp -m ppois -t ${OUT}_treat_pileup.bdg \
                -c ${OUT}_control_lambda.bdg -o ${OUT}_ppois.bdg -S $SF; \
            sortBed -i ${OUT}_ppois.bdg | slopBed -i - -g chrom.sizes \
                -b 0 | bedClip stdin chrom.sizes ${OUT}.pval.signal.bdg; \
            bedGraphToBigWig ${OUT}.pval.signal.bdg \
                chrom.sizes signal/${OUT}.pval.signal.bw; \
            rm -f ${OUT}*.bdg ${OUT}_model.r ${OUT}_peaks.xls \
                ${OUT}_summits.bed ${OUT}_peaks.narrowPeak " 
    done | parallelRun 

# **************** 12. IDR **************** 
# examples: 
## Rscript batch-consistency-analysis.r $myWorkDIR/macs2/SVP_rep1.peaks.bed $myWorkDIR/macs2/SVP_rep2.peaks.bed -1 $myWorkDIR/macs2/idr/SVP_rep1_VS_SVP_rep2 0 F p.value 
## Rscript batch-consistency-analysis.r $myWorkDIR/macs2/SVP_rep1.peaks.bed $myWorkDIR/macs2/SVP_rep3.peaks.bed -1 $myWorkDIR/macs2/idr/SVP_rep1_VS_SVP_rep3 0 F p.value 
## Rscript batch-consistency-analysis.r $myWorkDIR/macs2/SVP_rep2.peaks.bed $myWorkDIR/macs2/SVP_rep3.peaks.bed -1 $myWorkDIR/macs2/idr/SVP_rep2_VS_SVP_rep3 0 F p.value 
## Rscript batch-consistency-plot.r 1 $myWorkDIR/macs2/idr/SVP_repAll $myWorkDIR/macs2/idr/SVP_rep1_VS_SVP_rep3
# 
mkdir -p $myWorkDIR/macs2/idr && cd $myWorkDIR/bin/idrCode/ 
cp $myWorkDIR/macs2/chrom.sizes genome_table.txt 
pDIR=$myWorkDIR/macs2
# 1. Run IDR analysis on original replicates: 
awk -vFS="\t" -vOFS="\t" '{a[$1]=a[$1]?a[$1]"\t"$2:$2;} \
    END{for(i in a){print i,a[i]}}' $myWorkDIR/meta/macs2.meta \
    | sort -k 1 -k 2 | tee $myWorkDIR/meta/idr.meta \
    | awk -vFS="\t" -vDIR=$pDIR '(NF>2){for(i=2;i<=NF;i++) \
        {for(j=i+1;j<=NF;j++){print "Rscript \
        batch-consistency-analysis.r "DIR"/"$i".peaks.bed \
        "DIR"/"$j".peaks.bed -1 "DIR"/idr/"$1"_"$i"_VS_"$j" \
        0 F p.value"}}}' \
    | parallelRun 
# 2. Run IDR analysis on self pseudo-replicates
cut -f 2 $myWorkDIR/meta/macs2.meta | while read Rep; do 
        echo "Rscript batch-consistency-analysis.r \
            $pDIR/${Rep}_self_pr1.peaks.bed \
            $pDIR/${Rep}_self_pr2.peaks.bed \
            -1 $pDIR/idr/${Rep}_self_pr1_VS_pr2 0 F p.value"
    done | parallelRun 
# 3. Run IDR analysis on pooled pseudo-replicates
cut -f 1 $myWorkDIR/meta/macs2.meta | uniq | while read TF; do 
    echo "Rscript batch-consistency-analysis.r \
            $pDIR/${TF}_rep0_pr1.peaks.bed \
            $pDIR/${TF}_rep0_pr2.peaks.bed \
            -1 $pDIR/idr/${TF}_rep0_pr1_VS_pr2 0 F p.value"
    done | parallelRun 

# Get thresholds to truncate peak lists and final set of peak calls 
mkdir -p $myWorkDIR/macs2/final 
oDIR=$myWorkDIR/macs2/idr/
cat $myWorkDIR/meta/macs2.meta | while read TF Rep MORE; do 
        nfdr=$(awk '($9>5)' $myWorkDIR/macs2/${Rep}.peaks.bed | wc -l)
        # self-consistency threshold 
        nps_Self=$( awk '$11 <= 0.05 {print $0}' \
            $oDIR/${Rep}_self_pr1_VS_pr2-overlapped-peaks.txt \
            | wc -l )
        # original replicate threshold 
        max_nps_Rep=
        if [[ $(ls $oDIR/${TF}*${Rep}*-overlapped-peaks.txt \
                2> /dev/null | wc -l) -gt 0 ]]; then
            echo "*** $TF with replicates" 1>&2
            numPeaks=
            for cmp in $oDIR/${TF}*${Rep}*-overlapped-peaks.txt; do 
                numPeaks="$numPeaks "$( awk '$11 <= 0.02 \
                    {print $0}' $cmp | wc -l )
            done 
            max_nps_Rep=$(echo $numPeaks | awk \
                '{for(i=1;i<=NF;i++){if(!m || $i>m) \
                {m=$i}}}END{print m}')
        fi
        # pooled-consistency threshold 
        nps_Rep0=$( awk '$11 <= 0.01 {print $0}' \
            $oDIR/${TF}_rep0_pr1_VS_pr2-overlapped-peaks.txt \
            | wc -l )
        echo -e "$TF\t$Rep\t$nps_Self\t$max_nps_Rep\t$nps_Rep0\t$nfdr"
    done | tee $myWorkDIR/meta/idr.threshold \
    | awk -vFS="\t" -vOFS="\t" '{n[$1]++; \
        s[$1]=s[$1]?s[$1]"|"$3:$3; \
        if(!r[$1]||r[$1]<$4){r[$1]=$4} \
        if(!p[$1]||p[$1]<$5){p[$1]=$5} } \
        END{for(i in n){o=p[i];if(!o||o<r[i]){o=r[i]} \
            print i,o,n[i],s[i],r[i],p[i]}}' \
    | sort -k 1 | tee $myWorkDIR/meta/idr.optimal.threshold \
    | while read TF optThresh MORE; do 
        sort -k8nr,8nr $myWorkDIR/macs2/${TF}_rep0.peaks.bed \
        | head -n $optThresh | bedtools sort -i \
        > $myWorkDIR/macs2/final/$TF.peaks.bed 
    done 

# **************** 13. Peak annotation **************** 
# peak annotation using ChIPseeker: 
# http://bioconductor.org/packages/release/bioc/html/ChIPseeker.html  
cut -f 1 $myWorkDIR/meta/macs2.meta | uniq \
    | awk '{a=$1":"$1".peaks.bed";s=s?s"\t"a:a} \
        END{print "FlowerTF\t"s}' \
    | peakAnnotation -i $myWorkDIR/macs2/final \
    -o $myWorkDIR/annotation \
    -g $myWorkDIR/db/tair10Genes.gff3 

# compare peaks from the original publication
cd $myWorkDIR/macs2/final 
for peak in *.peaks.bed; do 
    TF=$(echo $peak | sed -e 's/.peaks.bed//')
    A=$(cat $peak | wc -l)
    B=$(cat $myWorkDIR/publication/$TF.bed | wc -l)
    C=$(bedtools intersect -u -a $peak -b $myWorkDIR/publication/$TF.bed | wc -l)
    D=$(bedtools intersect -u -b $peak -a $myWorkDIR/publication/$TF.bed | wc -l)
    echo -e "$TF\t$A\t$C\t$B\t$D"
done > $myWorkDIR/meta/peak_overlap_with_publication.txt 

# **************** 14. Target genes **************** 
# Genomic regions were associated with genes if located 3 kb upstream of the start of the gene up to 1 kb downstream of the end of the gene

mkdir -p $myWorkDIR/target && cd $myWorkDIR/target 
cp $myWorkDIR/meta/DE.gene.table.txt ./
awk -vFS="\t" -vOFS="\t" '{if($6=="+") \
    {$2=$2-3000;$3=$3+1000;if($2<0){$2=0;}} \
    else{$2=$2-1000;$3=$3+3000;if($2<0){$2=0;}} \
    print}' $myWorkDIR/db/tair10Genes.bed \
    > gene.3u1d.bed 

cut -f 1 $myWorkDIR/meta/idr.meta | while read TF; do 
        bedtools intersect -u -a gene.3u1d.bed \
            -b $myWorkDIR/macs2/final/$TF.peaks.bed \
        | cut -f4,5 | awk -vTF=$TF '{print TF"\t"$0}' \
        > $TF.target.list 
    done 

# **************** 15. DE genes **************** 
cat *.target.list | awk -vFS="\t" '{a[$2][$1]++;S[$1]++} \
    END{n=asorti(S);h=".";for(i=1;i<=n;i++){h=h"\t"S[i]} \
    print h; for(g in a){o=g;for(i=1;i<=n;i++) \
    {o=o"\t"(a[g][S[i]]+0)} print o;}}' | sort -k1 \
    > target.summary.txt

# R code 
R 
options(stringsAsFactors=FALSE)
targets <- read.table('target.summary.txt', head=T, row.names=1)
DEgenes <- read.table('../meta/DE.gene.table.txt', head=T, row.names=1,check.names=F)
DEgenes <- apply(DEgenes, c(1,2), function(x){r="NS";if(x>0)r="Up";if(x<0)r="Down";r})
pltdata <- c()
for(TF in colnames(targets)){
    day <- gsub(".*_day", "day", TF)
    target <- rownames(targets)[targets[, TF] > 0]
    number <- apply(DEgenes[rownames(DEgenes) %in% target, ], 2, table)
    dat <- number[, grep(paste0(day, "$"), colnames(number))]
    pltdata <- rbind(pltdata, dat[c("Down", "NS", "Up")])
}
rownames(pltdata) <- colnames(targets)
colnames(pltdata) <- c("Down", "NS", "Up")
pdf("Fig. 4.pdf", width=8.27, heigh=5, pointsize=10)
bx <- barplot(t(pltdata), col=c("#66bd63", "#878787", "#f46d43"), 
        border=NA, ylab="# genes", legend.text=colnames(pltdata), 
        args.legend=list(x="top", ncol=3, border=NA))
text(bx, rowSums(pltdata), rowSums(pltdata), pos=3, xpd=T)
dev.off()


# **************** 14. Target genes **************** 
cat *.target.list | tee target.list | awk -vFS="\t" 'BEGIN{ \
    while(getline<"../meta/TF.txt"){a[$1]++}}(a[$2])' > TFtarget.list

# R code 
R 
options(stringsAsFactors=FALSE)
library(igraph)
library(RColorBrewer)
## display.brewer.all()
TFtargets <- read.table("TFtarget.list", head=FALSE, sep="\t")
TFtargets[TFtargets[, 3]=="", 3] <- TFtargets[TFtargets[, 3]=="", 2]
net.data <- TFtargets[, c(1,3)]
regulators <- unique(net.data[,1])
vertices <- unique(c(net.data[,2], regulators))
net <- graph.data.frame(net.data, directed=T, vertices=vertices)

node.color <- rep(adjustcolor("grey80",alpha=.4), vcount(net))
node.size <- rep(2, vcount(net))
node.label <- rep(NA, vcount(net))
names(node.color) <- names(node.size) <- names(node.label) <- vertices

node.degree <- table(net.data[,2])
topTarget <- names(which(node.degree==max(node.degree)))
node.size[names(node.degree)] <- 2+node.degree
node.size[regulators] <- 12
node.color[names(node.degree)] <- colorRampPalette(brewer.pal(n=9,name="Greys")[1:7])(max(node.degree))[node.degree]
node.color[regulators] <- "#ed6d45" 
node.label[c(regulators, topTarget)] <- c(regulators, topTarget)

label.cex <- node.size
label.cex[] <- 0.4
label.cex[regulators] <- 0.8

V(net)$color <- node.color
V(net)$label <- node.label
V(net)$size <- node.size
V(net)$shape <- "circle"
V(net)$label.cex <- label.cex

E(net)$width <- 1 
E(net)$arrow.size <- .4
E(net)$arrow.width <- .4

mkgs <- list(grep("AP1_day", vertices), grep("SEP3_day", vertices))

pdf("Fig. 5.pdf", width=8.27, heigh=8.27, pointsize=10)
set.seed(1234)
plot(net, layout=layout.auto, mark.groups=mkgs, mark.col=c("#36a9e0","#36a9e0"), mark.border=NA, edge.arrow.size=.8, vertex.label.font=3)
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
dev.off()




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

Rscript $myWorkDIR/pipeline/corrplot.R
