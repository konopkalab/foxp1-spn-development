# Demultiplexing BCL Data

## DATASET | JAN 2022
### CELLRANGER
```{bash}
module load bcl2fastq/2.20.0
module load cellranger/3.0.2

cellranger mkfastq --id=DEMUX_JAN2022 \
                   --run=220113_A00672_0081_AH5VLKDSX3 \
                   --csv=samplesheet.csv \
                   --ignore-dual-index \
                   --use-bases-mask=Y151,I8n2,n10,Y151
```

### SAMPLESHEET
```{csv}
Lane,Sample,Index
2,NK3CF1,SI-GA-A1
2,NK3BM4,SI-GA-A2
2,NK3CF5,SI-GA-A3
2,NK3CF8,SI-GA-A4
2,NK3BF1,SI-GA-A5
2,NK3BM3,SI-GA-A6
2,NK3BF5,SI-GA-A7
2,NK3BM8,SI-GA-A8
2,NK3CF1_Amp,SI-GA-A9
2,NK3BM4_Amp,SI-GA-A10
2,NK3CF5_Amp,SI-GA-A11
2,NK3CF8_Amp,SI-GA-A12
2,NK3BF1_Amp,SI-GA-B1
2,NK3BM3_Amp,SI-GA-B2
2,NK3BF5_Amp,SI-GA-B3
2,NK3BM8_Amp,SI-GA-B4
```

<p>&nbsp;</p>

## DATASET | SEP 2022 
### BCL2FASTQ
```{bash}
module load bcl2fastq/2.20.0

FLOWCELL_DIR="220920_A00672_0105_BHY7VFDSX3"
OUTPUT_DIR="DEMUX_SEP2022"
INTEROP_DIR="00_B2FTEMP"
SAMPLE_SHEET_PATH="samplesheet.csv"
BASES_MASK="Y*,I*,I*,Y*"

bcl2fastq --use-bases-mask=${BASES_MASK} \
          --create-fastq-for-index-reads \
          --minimum-trimmed-read-length=8 \
          --mask-short-adapter-reads=8 \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 6 -w 6 \
          -R ${FLOWCELL_DIR} \
          --output-dir=${OUTPUT_DIR} \
          --interop-dir=${INTEROP_DIR} \
          --sample-sheet=${SAMPLE_SHEET_PATH}

```

### SAMPLESHEET
```{csv}
[Header]
EMFileVersion,4

[Reads]
28
151

[Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,Original_Sample_ID
2,WT_1,WT_1,TTATTCGAGG,AGCAGGACAG,HY7VFDSX3,WT_1
2,WT_2,WT_2,ATGGAGGGAG,AATGGGTTAT,HY7VFDSX3,WT_2
2,WT_3,WT_3,ACCAGACAAC,CCTAGTTCCT,HY7VFDSX3,WT_3
2,WT_1_Amp,WT_1_Amp,CGCGCACTTA,AGAATACAGG,HY7VFDSX3,WT_1_Amp
2,WT_2_Amp,WT_2_Amp,GCTACAAAGC,AGGGCACGTG,HY7VFDSX3,WT_2_Amp
2,WT_3_Amp,WT_3_Amp,TATCAGCCTA,AGGACGAAAC,HY7VFDSX3,WT_3_Amp
```

<p>&nbsp;</p>

# FastQC

### FastQC | Dataset | Jan 2022
```{bash}
module load fastqc/0.11.8

FQDIR="DEMUX_JAN2022"
cd ${FQDIR}

ls *.fastq.gz | sed "s/_001.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; fastqc %_001.fastq.gz'
```

<p>&nbsp;</p>

### FastQC | Dataset | Sep 2022
```{bash}
module load fastqc/0.11.8

FQDIR="DEMUX_SEP2022"
cd ${FQDIR}

ls *.fastq.gz | sed "s/_001.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; fastqc %_001.fastq.gz'
```

<p>&nbsp;</p>

