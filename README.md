
`Ayu` (not an acronym) is a subcellular location predictor for **marine** prokaryotic proteins, based on amino amino acid composition features. It is adapted for its use in large metagenomic datasets.


- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#documentation)
- [Usage](#usage)
- [Preprocessing](#preprocessing)
- [Running external programs: SignalP6, TMBed, IPC2](#running-external-programs-signalp6-tmbed-ipc2)
- [Prediction](#prediction)
- [Prediction output](#prediction-output)
- [Training and validation datasets] (#training-and-validation-datasets)
- [Citation](#citation)
- [License](#license)


# Overview
``Ayu`` aims to provide fast, reliable predictions of subecullar locations for proteins coded by marine prokaryotes. It is bases in XGboost, using as features the presence - absence of a signal peptide, amino acid composition information and transmembrane domain regions. The preprint can be found [here](https://www.researchsquare.com/article/rs-3585715/v1).

# Requirements
Ayu requires both [SignalP6.0](https://github.com/fteufel/signalp-6.0) and [TMBed](https://github.com/BernhoferM/TMbed) to be installed in order to use the `easy-workflow` command. If SignalP and TMBed are going to be run separately (See [Running SignalP6 and TMBed separately](#running-signalp6-and-tmbed-separately)), this is not required.

The python libraries required are:
```
biopython
pandas
numpy
scipy
scikit-bio
xgboost==1.7.6
```


# Installation

First, clone the repository and move to that folder:
```
git clone https://github.com/asierzaragoza/ayu && cd ayu
```

Then install with:
```
pip3 install .
```

# Usage
The program requires for prediction both internal (calculated by Ayu) and external (Calculated by other software) protein features.

## Preprocessing

Internal protein features have to be processed first with the command:
```
ayu preprocessing <fasta-file> <ayu_project_dir>
```
This command will also create an ayu project folder, which will be used in all following steps.


## Running external programs: SignalP6, TMBed, IPC2

Ayu provides a series of commands to run external programs. In order to use these commands, the program has to be in the $PATH or the path has to be provided by the ``--path`` flag:

```
ayu run_ipc2 <ayu_project_dir>
ayu run_tmbed <ayu_project_dir>
ayu run_signalp6 --path <path_to_signalp6_executable> <ayu_project_dir>
```
However, prediction of signal peptide and transmembrane regions is the most time- and resource-consuming step of the process. Therefore, some users might prefer to run TMBed and SignalP6 separately (or in a computer cluster) in order to speed up the prediction.

### SignalP6.0
The SignalP 6.0 results for ayu can be obtained with the command:
```
signalp6 --fastafile <protein_fasta_file> --output_dir <output_signalp_dir> --format none --mode fast
```
This command will result in a folder named <output_signalp_dir> with a `prediction_results.txt` file. The results can then be loaded into ayu with the command:
```
ayu load_signalp6 <output_dir>/prediction_results.txt <ayu_project_dir>
```
If SignalP6.0 has been ran multiple times (for example, if the fasta file has been divided into multiple smaller files for parallel processing), just run the command once for every output file:
```
ayu load_signalp6 <output_signalp_dir_1>/prediction_results.txt <ayu_project_dir>
ayu load_signalp6 <output_signalp_dir_2>/prediction_results.txt <ayu_project_dir>
```
### TMBed
The TMBed results for ayu can be obtained with the command:
```
tmbed predict -f <protein_fasta_file> -p <output_tmbed_file>  --out-format=3 --use-gpu
```
The resulting output file <output_tmbed_file> can be loaded into ayu with the command:
```
ayu load_tmbed <output_tmbed_file> <ayu_project_dir>
```
### IPC2
IPC pI predictions for ayu can be obtained with the command:
```
ipc2_protein_svr_predictor.py /<IPC2_models>/IPC2_protein_75_SVR_19.pickle <fasta_file> <ipc2_output>
```
The resulting file <output_ipc2_file> can be loaded into ayu with the command:
```
ayu load_ipc2 <output_ipc2_file> <ayu_project_dir>
```

## Prediction
Once all protein features have been processed, protein sublocations can be predicted with the command:
```
```
By default, Ayu will only run predictions for proteins with all the features already processed. To check which proteins do not have all the features, run the command: 
```
ayu predict --check_features <ayu_project_dir> <prediction_tsv_file>
```

This will output into STDOUT which proteins do not have all the features present in the project.

# Prediction output
The TSV file contains 4 columns:
[*] **prot_ID** protein header as found in the input fasta file.
[*] **cyto**: probability that the protein is located in the cytoplasm (0-1)
[*] **peri**: probability that the protein is located in the periplasm (0-1)
[*] **extr**: probability that the protein is secreted to the extracellular milieu (0-1)

Ayu does not provide a cut-off for predictions, but we recommend using a cut-off of at least 0.6.


# Training and validation datasets
Uniprot IDs of the proteins used for the validation and prediction datasets are contained in the `eval_train_sets` folder.

# Citation
If you use this tool for your research, please cite:


# License
This work is covered under the **GNU General Public License**.
