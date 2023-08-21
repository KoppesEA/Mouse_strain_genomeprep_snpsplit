#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -J Samtools_mGRC38
#SBATCH --output=CastEiJ_Samtools_out-%A_%a.txt
#SBATCH --cpus-per-task=8 
#SBATCH --mem=12g  ##

## BASH SCRIPT TO EXTRACT GRCm39_109 COORDINATES DEFINED IN BEDFILE FROM REF(B6/C57), CAST AND N-MASKED GENOMES

## Define Variables
WKDIR=/ix1/mmann/KoppesEA/REF_Sequences/Mus_musculus
CASTDIR=$WKDIR/Cast_EiJ/CastEiJ_GRCm39
REF_GRCm39v109=$WKDIR/GRCm39_ref/Mus_musculus.GRCm39.dna.toplevel.fa #
CAST_GRCm39v109=$CASTDIR/CAST_EiJ_full_sequence/Mus_musculus.GRCm39.109.dna.CASTEiJ.fa  ## snp-split generated .fa Cast and N-masked .txt files concatenated
NMASK_GRCm39v109=$CASTDIR/CAST_EiJ_N-masked/Mus_musculus.GRCm39.109.dna.NMASK.fa
COORD_LIST=$CASTDIR/Nup_mm11_GRCm39_genecoord.txt
OUT_FA=$CASTDIR/Nup_mm11_GRCm39_genecoord.txt.fa
REF_base=`basename $REF_GRCm39v109 .fa`
REF_sort=$WKDIR/GRCm39_ref/${REF_base}_sorted.fa
CAST_base=`basename $CAST_GRCm39v109 .fa`
CAST_sort=$CASTDIR/CAST_EiJ_full_sequence/${CAST_base}_sorted.fa
NMASK_base=`basename $NMASK_GRCm39v109 .fa`
NMASK_sort=$CASTDIR/CAST_EiJ_N-masked/${NMASK_base}_sorted.fa


## Import Nup107 peaks formatted as manually, remove ## comment out lines
coords=($(cat $COORD_LIST))
echo -e "Reporting Ref, CAST and N-masked sequences for:\n" $REF_GRCm39v109 "\n" $CAST_GRCm39v109 "\n" $NMASK_GRCm39v109 "\nbased on:\n" $COORD_LIST 

## Load Samtools
module load gcc/8.2.0
module load samtools/1.14

#`_ comment out if done before
## Samtools index .fa from Ref, Cast and N-masked GRCm38_68 snp-split generated genomes
echo "Sorting and Indexing REF, CAST and NMASK genomes"
samtools sort -o $REF_sort $REF_GRCm39v109
samtools sort -o $CAST_sort $CAST_GRCm39v109
samtools sort -o $NMASK_sort $NMASK_GRCm39v109
samtools index $REF_sort
samtools index $CAST_sort
samtools index $NMASK_sort
#`

rm $OUT_FA
touch $OUT_FA
echo -e "removing and making:\n" $OUT_FA

## Samtools faidx extract chromosomal coordinates
## 3 column tsv with name coord strand as input
## while with if else loop to get reverse comp if feature on (-) strand as specified by faidx -i option
while IFS=$'\t' read -r -a rec; do
    NAME=${rec[0]}
    COORD=${rec[1]}
    STRAND=${rec[2]}
    echo "Processing Record NAME: " $NAME
    echo "Processing Record COORD: " $COORD ":" $STRAND
    echo "Processing Record NAME: " $NAME >> $OUT_FA
    echo ">REFB6" $NAME >> $OUT_FA
    if [ "$STRAND" -eq 1 ]; then
    	samtools faidx $REF_sort $COORD >> $OUT_FA
    	echo ">CASTEiJ_" $NAME >> $OUT_FA
		samtools faidx $CAST_sort $COORD >> $OUT_FA
		echo ">NMASK_" $NAME >> $OUT_FA
		samtools faidx $NMASK_sort $COORD >> $OUT_FA
	elif [ "$STRAND" -eq -1 ]; then
	    samtools faidx -i --mark-strand sign $REF_sort $COORD >> $OUT_FA
    	echo ">CASTEiJ_" $NAME >> $OUT_FA
		samtools faidx -i --mark-strand sign $CAST_sort $COORD >> $OUT_FA
		echo ">NMASK_" $NAME >> $OUT_FA
		samtools faidx -i --mark-strand sign $NMASK_sort $COORD >> $OUT_FA
	else
		echo "Strand-Specifity not set, expecting 3 column tsv with name coord strand as input, Processing as + strand"
		samtools faidx $REF_sort $COORD >> $OUT_FA
    	echo ">CASTEiJ_" $NAME >> $OUT_FA
		samtools faidx $CAST_sort $COORD >> $OUT_FA
		echo ">NMASK_" $NAME >> $OUT_FA
		samtools faidx $NMASK_sort $COORD >> $OUT_FA
	fi
done < $COORD_LIST