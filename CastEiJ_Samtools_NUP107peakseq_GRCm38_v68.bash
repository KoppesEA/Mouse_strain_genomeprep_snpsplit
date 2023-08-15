#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -J Samtools_mGRC38
#SBATCH --output=CastEiJ_Samtools_out-%A_%a.txt
#SBATCH --cpus-per-task=8 
#SBATCH --mem=12g  ##

## BASH SCRIPT TO EXTRACT GRCm38_68 COORDINATES DEFINED IN BEDFILE FROM REF(B6/C57), CAST AND N-MASKED GENOMES

## Define Variables
WKDIR=/ix1/mmann/KoppesEA/REF_Sequences/Mus_musculus
CASTDIR=$WKDIR/Cast_EiJ/CastEiJ_GRCm38
REF_GRCm38v68=$WKDIR/GRCm38_ref/GRCm38_68/Mus_musculus.GRCm38.68.dna.toplevel.fa  #
CAST_GRCm38v68=$CASTDIR/CAST_EiJ_full_sequence/Mus_musculus.GRCm38.68.dna.CASTEiJ.fa  ## snp-split generated .fa Cast and N-masked .txt files concatenated
NMASK_GRCm38v68=$CASTDIR/CAST_EiJ_N-masked/Mus_musculus.GRCm38.68.dna.NMASK.fa
COORD_LIST=$CASTDIR/Nup107_mm10peak_GRCm38_v68_imprint.txt
OUT_FA=$CASTDIR/Mmus_GRCm38_68_REFCASTNMASK_Nup107imprint.fa

## Import Nup107 peaks formatted as manually, remove ## comment out lines
coords=($(cat $COORD_LIST))
echo -e "Reporting Ref, CAST and N-masked sequences for:\n" $REF_GRCm38v68 "\n" $CAST_GRCm38v68 "\n" $NMASK_GRCm38v68 "\nbased on:\n" $COORD_LIST 

## Load Samtools
module load gcc/8.2.0
module load samtools/1.14

#`_ comment out if done before
## Samtools index .fa from Ref, Cast and N-masked GRCm38_68 snp-split generated genomes
echo "Indexing REF, CAST and NMASK genomes"
samtools index $REF_GRCm38v68
samtools index $CAST_GRCm38v68
samtools index $NMASK_GRCm38v68
#`

rm $OUT_FA
touch $OUT_FA
echo -e "removing and making:\n" $OUT_FA

## Samtools faidx extract chromosomal coordinates
while IFS=$'\t' read -r -a rec; do
    NAME=${rec[0]}
    COORD=${rec[1]}
    echo "Processing Record NAME: " $NAME
    echo "Processing Record COORD: " $COORD
    echo "Processing Record NAME: " $NAME >> $OUT_FA
    echo ">REFB6" $NAME >> $OUT_FA
    samtools faidx $REF_GRCm38v68 $COORD >> $OUT_FA
    echo ">CASTEiJ_" $NAME >> $OUT_FA
	samtools faidx $CAST_GRCm38v68 $COORD >> $OUT_FA
	echo ">NMASK_" $NAME >> $OUT_FA
	samtools faidx $NMASK_GRCm38v68 $COORD >> $OUT_FA    
done < $COORD_LIST
