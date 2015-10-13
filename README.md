# SCITE



## Description


**SCITE** is a software to compute mutation histories of somatic cells.
Given noisy mutation profiles of single cells, **SCITE** performs a stochastic
search to find the Maximum Likelihood tree and to sample from the posterior
probability distribution. Tree reconstruction can be combined with an estimation
of the error rates in the mutation profiles.

**SCITE** is particularly designed for reconstructing mutation histories of
tumours based on mutation profiles obtained from single-cell exome sequencing experiments, but is in 
principle applicable to any type of (noisy) mutation profiles for which the
infinite sites assumption can be made.

## Availability

**SCITE** is freely available under a GPL3 license at https://gitlab.com/jahnka/SCITE

##    How to run **SCITE**



### Mac OS X

To compile the C/C++ program, open a terminal and go to the folder containing the source files, and type

	clang++ findBestTrees.cpp matrices.cpp mcmc.cpp output.cpp rand.cpp scoreTree.cpp treelist.cpp trees.cpp -o scite

This writes a file named `scite`. Assuming the sample data file dataKimSimon.csv is located in the same folder, `scite` can then be executed as follows

	./scite -i dataKimSimon.csv -n 18 -m 58 -r 1 -l 900000 -g 1.25 -fd 6.04e-5 -ad 0.21545 0.21545 -cc 1.299164e-05



### Linux/Unix

To compile the C/C++ program, open a terminal and go to the folder containing the source files, and type

	g++ findBestTrees.cpp matrices.cpp mcmc.cpp output.cpp rand.cpp scoreTree.cpp treelist.cpp trees.cpp -o scite

This writes a file named `scite`. Assuming the sample data file dataKimSimon.csv is located in the same folder, `scite` can then be executed as follows

	./scite -i dataKimSimon.csv -n 18 -m 58 -r 1 -l 900000 -g 1.25 -fd 6.04e-5 -ad 0.21545 0.21545 -cc 1.299164e-05


##  Input Files


### 1. Mutation Matrix


Each row specifies the mutation profile of a single cell, and each column
represents one mutation.

#### (a) Only absence/presence of mutation is distinguished
The entry at position [i,j] should be

* 0 if mutation j is not observed in cell i,
* 1 if mutation j is observed in cell i, or
* 3 if the data point is missing

	
#### (b) Heterozygous and homozygous mutations distinguished

* 0 if mutation j is not observed in cell i,
* 1 if heterozygous mutation j is observed in cell i
* 2 if homozygous mutation j is observed in cell i
* 3 if the data point is missing

### 2. Mutation names (optional)


A list specifying the names the mutations, e.g. the name of the gene in which
the mutation occurs. If not specified, the mutations are numbered from 1 to n.

### 3. The true tree (optional)


If the true mutation tree is known (as in simulation experiments), the file with
the true tree (in GraphViz format) can be specified to compare the
predicted trees internally with the true tree.

##  Output Files
ML trees are each written to files in GraphViz and Newick format. Files are numbered
consecutively (e.g. dataKimSimon_ml1.gv, dataKimSimon_ml1.newick, dataKimSimon_ml2.gv, dataKimSimon_ml2.newick, ...). The base name of the output file is derived from the name of the input file (unless a different name is specified via `-o <filename>`).


## Parameters

`-i <filename>`     Replace \<filename\> with the file containing the mutation matrix

`-n <INT>`  Replace \<INT\> with the number of mutations (columns) in the dataset.

`-m <INT>`  Replace \<INT\> with the  number of cells (rows) in the dataset.

`-r <INT>`  Set \<INT\> to the desired number of repetitions of the MCMC.

`-l <INT>`  Set \<INT\> to the desired chain length of each MCMC repetition

`-g <DOUBLE>` Set \<DOUBLE\> to the desired value of gamma for ML computation (gamma > 1: fast convergence, possibly local optimum; gamma < 1: slower convergence); gamma needs to equal 1 to be able to sample from the posterior distribution



##### In case only absence/presence of a mutation is distinguished

`-fd <DOUBLE>` Set \<DOUBLE\> to the estimated false positive rate (false discoveries) of the sequencing experiment.

`-ad <DOUBLE>` Set \<DOUBLE\> to the estimated false negative rate (allelic dropout) of the sequencing experiment.

##### In case heterozygous and homozygous mutations are distinguished

`-fd <DOUBLE>` Set \<DOUBLE\> to the estimated false positive rate (false calling of heterozygous mutation) of the sequencing. experiment

`-ad <DOUBLE> <DOUBLE>` Set the first \<DOUBLE\> to the estimated rate of missed heterozygous mutations in the sequencing experiment and set the second \<DOUBLE\> to the estimated rate of homozygous mutations called as heterozygous mutations.

`-cc <DOUBLE>` Set \<DOUBLE\> to the estimated rate of homozygous mutations missed in the sequencing experiment.



##### Optional parameters

`-s` Setting this option causes the sample attachment points (i. e. where the cells would attach to the tree) to be marginalized out.

`-a` When setting this option, **SCITE** adds the individual cells as additional nodes (leafs) to the reported trees. Cells are attached where they fit best in the given tree (with respect to the error rates). By default, only the mutation tree is reported.

`-t <filename>`  Replace \<filename\> with a file containing the true tree in GraphViz format.

`-names <filename>` Replace \<filename\> with a file listing the mutation names.

`-o <filename>`   Replace \<filename\> with the desired base of the output file to overwrite the default output file names. 



