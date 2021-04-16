# CADEGA
CADEGA is a framework for identifying subnetworks of genes whose mutations are most associated with an increase in drug sensitivity. This offers key insights towards work evaluating cancer dependency genes, via a similar framework as NETPHIX ( ). _This project is in pre-production._
## Setup
The set up process for CADEGA requires the following steps.
## Download
Download the uploaded scripts in the CADEGA folder.
## Installation of Sub-Components
The main project that this framework depends on is NetMix, available from GitHub ( ).
Note that NetMix has the following subcomponents itself: Linux or Windows Subsystem for Linux, NumPy, SciPy, h5py, and heinz.
## Use
To run a CADEGA analysis, currently one will have to run the subsets of data separately. Note that his is a consequence of the system being, as mentioned before, in pre-production.

(from Windows or Linux)
[drugset] = PADLL.Integrator.scanDrugSet()
PADLL.Main([drugset])

(from Linux)
./netmix.sh ../network.tsv ../inputs/CDInput##.tsv test/outputfile.txt

## On Inputs
The network.tsv is a gene-to-gene interaction network from the STRING database, all with significance ratings at or above 0.9.
The PADLL input drug sets are the standard .csv output from Cancer RX Gene, with sensitivies here being recorded as AUC values across patient cell-lines and tissues.
The NetMix input is a gene to score table (with no header), with normalized, sensitivity correlated scores attached for each type of gene alteration.
