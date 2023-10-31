# Build Reference Genome
## Reference genome for primary samples
```{bash}
## GENOME:  MM10 GRCm38p6
## GTF:     GENCODE vM17

module load cellranger/3.0.2

## Extracting rows with feature type 'transcript'
grep -P '\ttranscript\t' gencode.vM17.annotation.gtf > gencode.vM17.annotation.premrna.gtf

## Replacing feature type 'transcript' with 'exon'
sed -i 's/\ttranscript\t/\texon\t/' gencode.vM17.annotation.premrna.gtf

genomefa="MM10_GRCm38p6/GRCm38.p6.genome.selected.fa"
annotgtf="gencode.vM17.annotation.gtf"
newannotgtf="gencode.vM17.annotation.premrna.gtf"

## Generating genome index
cellranger mkref --genome=MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER \
                 --fasta=${genomefa} \
                 --genes=${newannotgtf}
```

<p>&nbsp;</p>

## Reference genome for amplified (WPRE, MCHERRY) samples
```{bash}
## GENOME:  MM10 GRCm38p6
## GTF:     GENCODE vM17

cat mm10_genome.fa wpre_mcherry_sequenced.fa > mm10_genome_wpre_mcherry.fa
cat mm10_genes.gtf wpre_mcherry_sequenced.gtf > mm10_genome_wpre_mcherry.gtf

genomefa="mm10_genome_wpre_mcherry.fa"
annotgtf2="mm10_genome_wpre_mcherry.gtf"

echo "Generating Genome Index for CellRanger"
cellranger mkref --genome=MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER_AMP \
                 --fasta=${genomefa} \
                 --genes=${annotgtf2}
```

<p>&nbsp;</p>

