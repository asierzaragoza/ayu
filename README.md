
`Ayu` (not an acronym) is a subcellular location predictor for **marine** prokaryotic proteins, based on amino amino acid composition features. It is adapted for its use in large metagenomic datasets.


- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#documentation)
- [Usage](#usage)
- [Preprocessing](#preprocessing)
- [Running external programs: SignalP6, TMBed, IPC2](#running-external-programs-signalp6-tmbed-ipc2)
- [Prediction](#prediction)
- [Prediction output](#prediction-output)
- [Citation](#citation)
- [License](#license)


# Overview
``Ayu`` aims to provide fast, reliable

# Requirements
Ayu requires both [SignalP6.0](https://github.com/fteufel/signalp-6.0) and [TMBed](https://github.com/BernhoferM/TMbed) to be installed in order to use the `easy-workflow` command. If SignalP and TMBed are going to be run separately (See [Running SignalP6 and TMBed separately](#running-signalp6-and-tmbed-separately)), this is not required.

# Installation

# Usage

## Preprocessing

## Running external programs: SignalP6, TMBed, IPC2
Prediction of signal peptide and transmembrane regions is the most time- and resource-consuming step of the process. Therefore, some users might prefer to run TMBed and SignalP6 separately (or in a computer cluster) in order to speed up the prediction.

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

# Prediction output

# Citation
If you use this tool for your research, please cite:

# License
This work is covered under the **GNU General Public License**.
