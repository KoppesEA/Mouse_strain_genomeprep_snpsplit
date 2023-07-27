# Mouse_strain_genomeprep_snpsplit
GRCm38 and GRCm39 strain-specific reference genome generated with snp-split

## Summarry
The following commands were run to download appropriate GRCm38 and GRCm39 genomes and then to generate Castenous (Cast_EiJ) strain-specific genomes. While this largely followed the snpsplit genome gen guide (see references) there were some specific issues getting the backdated GRCm38 to work. This was run July 26th 2023 on the PITT HTC.

### Obtain GRCm38 and GRCm39 genomes and place them in separate directories
`wget -O Mus_musculus.GRCm39.dna.toplevel.fa.gz https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz &`
`wget -c -O Mus_musculus.GRCm38.68.dna.toplevel.fa.gz ftp://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz`
The GRCm38.68 ref was required for the v5 SNP build (other builds did not work), I could not access via the sanger.ac.uk link from legacy SNPsplit forum. Then unzip each within their folder.

### Obtain GRCm38 and m39 strain specific SNPs in VCF format in another directory
`wget -O mgp_REL2021_snps.vcf.gz ftp://ftp.sra.ebi.ac.uk/vol1/ERZ840/ERZ8406694/mgp_REL2021_snps.vcf.gz &`
`wget -O mgp.v5.merged.snps_all.dbSNP142.vcf.gz https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz &`

### Run SNPsplit_genome_preparation using SLURM
Run the Bash scripts `CastEiJ_snpsplit_genomegen_GRCm39.bash` and `CastEiJ_snpsplit_genomegen_GRCm38_v68` including output for strain-specific genome

### Extraction of Genomic Sequences of regions of interestr of Reference, N-masked and CastEiJ using Samtools
TBD--inprogress

## References and Links
https://felixkrueger.github.io/SNPsplit/genome_prep/legacy/
https://ftp.ebi.ac.uk/pub/databases/mousegenomes/
https://www.sanger.ac.uk/data/mouse-genomes-project/
https://www.mousegenomes.org/snps-indels/
