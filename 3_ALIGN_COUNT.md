
# Alignment and Count 

## Dataset | Jan 2022 | Primary Samples
```{bash}
module add samtools/1.6
module load cellranger/3.0.2

GENOMEINDEX="MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER"
FQDIR="DEMUX_JAN2022"
CNTDIR="COUNT_JAN2022"

cd ${CNTDIR}

for lib in NK3BF1 NK3BF5 NK3BM3 NK3BM4 NK3BM8 NK3CF1 NK3CF5 NK3CF8
	do 
		echo "PROCESSING ", ${lib}

		cellranger count --id=${lib} \
						 --fastqs=${FQDIR} \
						 --sample=${lib} \
						 --transcriptome=${GENOMEINDEX} \
						 --expect-cells=10000 \
						 --chemistry=SC3Pv3
	
	done
```

<p>&nbsp;</p>


## Dataset | Jan 2022 | Amplified Samples
```{bash}
module add samtools/1.6
module load cellranger/3.0.2

GENOMEINDEX="MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER_AMP"
FQDIR="DEMUX_JAN2022"
CNTDIR="COUNT_JAN2022_AMP"

cd ${CNTDIR}

for lib in NK3BF1c NK3BF5c NK3BM3c NK3BM4c NK3BM8c NK3CF1c NK3CF5c NK3CF8c
	do 
		echo "PROCESSING ", ${lib}

		cellranger count --id=${lib} \
						 --fastqs=${FQDIR} \
						 --sample=${lib} \
						 --transcriptome=${GENOMEINDEX} \
						 --expect-cells=10000 \
						 --chemistry=SC3Pv3
	
	done
```

<p>&nbsp;</p>


## Dataset | Sep 2022 | Primary Samples
```{bash}
module add samtools/1.6
module load cellranger/3.0.2

GENOMEINDEX="MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER"
FQDIR="DEMUX_SEP2022"
CNTDIR="COUNT_SEP2022"

cd ${CNTDIR}

for lib in WT_1 WT_2 WT_3
	do 
	echo "PROCESSING ", ${lib}

	cellranger count --id=${lib} \
					 --fastqs=${FQDIR} \
					 --sample=${lib} \
					 --transcriptome=${GENOMEINDEX} \
					 --expect-cells=10000

	done
```

<p>&nbsp;</p>


## Dataset | Sep 2022 | Amplified Samples
```{bash}
module add samtools/1.6
module load cellranger/3.0.2

GENOMEINDEX="MM10_GRCm38p6_GENCODEvM17_PREMRNA_CELLRANGER_AMP"
FQDIR="DEMUX_SEP2022"
CNTDIR="COUNT_SEP2022_AMP"

cd ${CNTDIR}

for lib in WT_1_Amp WT_2_Amp WT_3_Amp
	do 
	echo "PROCESSING ", ${lib}

	cellranger count --id=${lib} \
					 --fastqs=${FQDIR} \
					 --sample=${lib} \
					 --transcriptome=${GENOMEINDEX} \
					 --expect-cells=10000

	done
```

<p>&nbsp;</p>

