# graph_freq_MS
graph frequency analysis of brain fingerprints for multiple sclerosis

#### Brief intro
This repository is associated with the paper "Coupling of brain activity and structural network in clinical isolated syndrome and multiple sclerosis: a graph frequency analysis study".
We adopted graph frequency analysis to quantify how much brain regional functional activity with different graph frequency are organized atop the underlying wiring diagram in clinical isolated syndrome and multiple sclerosis.

#### Brain segmentation and tractograph
`mrtrix3_connectome.py`: based on packages inside folder `mrtrix3`.


#### Surrogate functional signal creation
`surrogate_BOLD_create.py`: python script modified from [Preti's study](https://doi.org/10.1038/s41467-019-12765-7).

#### Parcel and ICN lookup table
`Yeo-hcpmmp-label`

#### Statistics 
`GSP_study_Neuroimmune.Rmd`: Rmarkdown script for group analysis

#### Brain genetic-imaging colocalization analysis
`AHBA_HCPMMP_MS.R`: R script for analysis and plotting.
`targets_associated_with_multiple_sclerosis.csv`: MS-associated genes from [the Open Targets Platform](https://www.targetvalidation.org/). 
