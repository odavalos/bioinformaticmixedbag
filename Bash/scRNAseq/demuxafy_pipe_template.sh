#!/bin/bash

# Running Demuxafy

#  BAM = BAM file
#  VCF = VCF file
#  BARCODES = Barcode file 
#  N = Number of samples
#  FASTA = Reference genome
#  THREADS = Number of threads
#  FREEMUXLET_OUTDIR = Output directory for Freemuxlet
#  SOUPORCELL_OUTDIR = Output directory for Souporcell
#  VIREO_OUTDIR = Output directory for Vireo
#  SCSPLIT_OUTDIR = Output directory for SCSplit



#  Input Paths
#  ===================================================================
BAM=
VCF=
BARCODES=
N=
FASTA=
THREADS=16

# Output Directories
FREEMUXLET_OUTDIR=
SOUPORCELL_OUTDIR=
VIREO_OUTDIR=
SCSPLIT_OUTDIR=
#  ===================================================================

# Print input paths
echo "----------------------------------------"
echo $'\n'
echo $' BAM: $BAM\n'
echo $' VCF: $VCF\n'
echo $' BARCODES: $BARCODES\n'
echo $' N: $N\n'
echo $' FASTA: $FASTA\n'
echo $' THREADS: $THREADS\n'
echo $' FREEMUXLET_OUTDIR: $FREEMUXLET_OUTDIR\n'
echo $' SOUPORCELL_OUTDIR: $SOUPORCELL_OUTDIR\n'
echo $' VIREO_OUTDIR: $VIREO_OUTDIR\n'
echo $' SCSPLIT_OUTDIR: $SCSPLIT_OUTDIR\n'
echo $'\n'
echo "----------------------------------------"

# export SINGULARITY_BIND="/scratch:/demuxafy/data"

# calculte run time
start=$(date +%s)


#  ===================================================================
#  Freemuxlet
#  ===================================================================
echo $'\nRunning Freemuxlet . . . \n'

# calculte run time
start_freemuxlet=$(date +%s)

# Popscle Pileup

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif popscle dsc-pileup \
	--sam $BAM \
	--vcf $VCF \
	--group-list $BARCODES \
	--out $FREEMUXLET_OUTDIR/pileup


# Popscle Freemuxlet

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif popscle freemuxlet \
	--plp $FREEMUXLET_OUTDIR/pileup \
	--out $FREEMUXLET_OUTDIR/freemuxlet \
	--group-list $BARCODES \
	--nsample $N


# Freemuxlet Summary

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif bash Freemuxlet_summary.sh $FREEMUXLET_OUTDIR/freemuxlet.clust1.samples.gz > $FREEMUXLET_OUTDIR/freemuxlet_summary.tsv

echo $'\nFreemuxlet complete! \n'

end_freemuxlet=$(date +%s)
runtime_freemuxlet=$((end_freemuxlet-start_freemuxlet))

#  ===================================================================
#  Souporcell
#  ===================================================================
echo $'\nRunning Souporcell . . . \n'

# calculte run time
start_souporcell=$(date +%s)

# Run Souporcell

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif souporcell_pipeline.py \
	-i $BAM \
	-b $BARCODES \
	-f $FASTA \
	-t $THREADS \
	-o $SOUPORCELL_OUTDIR \
	-k $N \
	--common_variants $VCF

# Souporcell Summary

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif bash souporcell_summary.sh $SOUPORCELL_OUTDIR/clusters.tsv > $SOUPORCELL_OUTDIR/souporcell_summary.tsv

echo $'\nSouporcell complete! \n'

end_souporcell=$(date +%s)
runtime_souporcell=$((end_souporcell-start_souporcell))

#  ===================================================================
#  Vireo
#  ===================================================================
echo $'\nRunning Vireo . . . \n'

# calculte run time
start_vireo=$(date +%s)

# CellSNP Pileup
singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif cellsnp-lite \
	-s $BAM \
	-b $BARCODES \
	-O $VIREO_OUTDIR \
	-R $VCF \
	-p 20 \
	--minMAF 0.1 \
	--minCOUNT 20 \
	--gzip

# Demultiplex with Vireo

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif vireo \
	-c $VIREO_OUTDIR \
	-o $VIREO_OUTDIR \
	-N $N

echo $'\nVireo complete! \n'

end_vireo=$(date +%s)
runtime_vireo=$((end_vireo-start_vireo))

#  ===================================================================
#  scSplit
#  ===================================================================
echo $'\nRunning scSplit . . . \n'

# calculte run time
start_scsplit=$(date +%s)

# Prepare bam file
singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif samtools view \
	-b \
	-S \
	-q 10 \
	-F 3844 $BAM > $SCSPLIT_OUTDIR/filtered_bam.bam

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif samtools rmdup $SCSPLIT_OUTDIR/filtered_bam.bam $SCSPLIT_OUTDIR/filtered_bam_dedup.bam

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif samtools sort \
	-o $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam $SCSPLIT_OUTDIR/filtered_bam_dedup.bam

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif samtools index $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam

# Call Sample SNVs
singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif freebayes \
	-f $FASTA \
	-iXu \
	-C 2 \
	-q 1 $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam > $SCSPLIT_OUTDIR/freebayes_var.vcf

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif vcftools \
	--gzvcf $SCSPLIT_OUTDIR/freebayes_var.vcf \
	--minQ 30 \
	--recode \
	--recode-INFO-all \
	--out $SCSPLIT_OUTDIR/freebayes_var_qual30

# Demultiplex with ScSplit
singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif scSplit count \
	-c $VCF \
	-v $SCSPLIT_OUTDIR/freebayes_var_qual30.recode.vcf \
	-i $SCSPLIT_OUTDIR/filtered_bam_dedup_sorted.bam \
	-b $BARCODES \
	-r $SCSPLIT_OUTDIR/ref_filtered.csv \
	-a $SCSPLIT_OUTDIR/alt_filtered.csv \
	-o $SCSPLIT_OUTDIR

singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif scSplit run \
	-r $SCSPLIT_OUTDIR/ref_filtered.csv \
	-a $SCSPLIT_OUTDIR/alt_filtered.csv \
	-n $N \
	-o $SCSPLIT_OUTDIR

singularity exec Demuxafy.sif scSplit genotype \
	-r $SCSPLIT_OUTDIR/ref_filtered.csv \
	-a $SCSPLIT_OUTDIR/alt_filtered.csv \
	-p $SCSPLIT_OUTDIR/scSplit_P_s_c.csv \
	-o $SCSPLIT_OUTDIR

# scSplit Summary
singularity exec --cleanenv -H $PWD -B $PWD Demuxafy.sif bash scSplit_summary.sh $SCSPLIT_OUTDIR/scSplit_result.csv > $SCSPLIT_OUTDIR/scSplit_summary.tsv

echo $'\nscSplit complete! \n'

end_scsplit=$(date +%s)
runtime_scsplit=$((end_scsplit-start_scsplit))

# total runtime
end = $(date +%s)
runtime = $((end-start))

# ===================================================================

# Print runtimes
echo "----------------------------------------"
echo $'Pipeline runtimes: \n'
echo "----------------------------------------"
echo $'Demuxafy: $runtime_demuxafy seconds \n'
echo $'Freemuxlet: $runtime_freemuxlet seconds \n'
echo $'Souporcell: $runtime_souporcell seconds \n'
echo $'Vireo: $runtime_vireo seconds \n'
echo $'scSplit: $runtime_scsplit seconds \n'
echo "----------------------------------------"
echo $'Total runtime: $runtime seconds \n'
echo "----------------------------------------"

