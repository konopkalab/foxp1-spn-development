# Removing Ambient RNA using CellBender

## Setup conda environment and install CellBender
```{bash}
# conda create -n CellBender python=3.7
# conda activate CellBender
# conda install -c anaconda pytables
# conda install pytorch torchvision -c pytorch
# git clone https://github.com/broadinstitute/CellBender.git
# pip install -e CellBender
# conda deactivate
```

## Run CellBender | Jan 2022 Dataset | Primary Samples
```{bash}
conda activate CellBender

ls | grep -P "^NK3BF" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000 --total-droplets-included 70000 --fpr 0.1'
ls | grep -P "^NK3BM" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000 --total-droplets-included 70000 --fpr 0.1'
ls | grep -P "^NK3CF" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000 --total-droplets-included 70000 --fpr 0.1'

conda deactivate
```

<p>&nbsp;</p>

## Run CellBender | Jan 2022 Dataset | Primary Samples
```{bash}
conda activate CellBender

cellbender remove-background --input WT_1_raw_feature_bc_matrix.h5 --output WT_1_CellBender_Out.h5 --expected-cells 10000 --total-droplets-included 50000 
cellbender remove-background --input WT_2_raw_feature_bc_matrix.h5 --output WT_2_CellBender_Out.h5 --expected-cells 10000 --total-droplets-included 50000
cellbender remove-background --input WT_3_raw_feature_bc_matrix.h5 --output WT_3_CellBender_Out.h5 --expected-cells 10000 --total-droplets-included 50000

conda deactivate

```

<p>&nbsp;</p>

