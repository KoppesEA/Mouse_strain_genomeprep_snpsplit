#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -J SNPsplit_mGRC39
#SBATCH --output=CastEiJ_snpsplit_genomegen_out-%A_%a.txt
#SBATCH --cpus-per-task=16
#SBATCH --mem=24g 

#BASH SCRIPT TO HAVE snpsplit CastEiJ SNPs split and genome generate

WKDIR=/ix1/mmann/KoppesEA/REF_Sequences/Mus_musculus
REFDIR=$WKDIR/GRCm39_ref  #needs to contain unzipped ref .fa
CASTDIR=$WKDIR/Cast_EiJ/CastEiJ_GRCm39  ##put both .bash script and vcf.gz into CASTDIR
VCF_FILE=$CASTDIR/mgp_REL2021_snps.vcf.gz 

echo "REFDIR is defined as: " $REFDIR
echo "CASTDIR is defined as: " $CASTDIR
echo "VCF_FILE is defined as: " $VCF_FILE

module load snpsplit/0.6.0

SNPsplit_genome_preparation --full_sequence --vcf_file $VCF_FILE --reference_genome $REFDIR --strain CAST_EiJ
