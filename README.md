# SCITE
========

## Update
---------
Please be aware that the software (https://github.com/cbg-ethz/infSCITE) accompanying our publication

Kuipers J et al. Single-cell sequencing data reveals widespread recurrence and loss of mutational hits in the life histories of tumorus. Genome Research 2017; 27:1885-1894.

is an extended version of SCITE that has a couple of new features:

* Explicit modelling of doublets
* Option to specify one recurrent mutation/one mutation loss.
* Learning  of the false positive rate for panel sequencing data


## Description
--------------


**SCITE** is a software package to compute mutation histories of somatic cells.
Given noisy mutation profiles of single cells, **SCITE** performs a stochastic
search to find the Maximum Likelihood (ML) or Maximum aposterori (MAP) tree and/or to sample from the posterior
probability distribution. Tree reconstruction can be combined with an estimation
of the error rates in the mutation profiles.

**SCITE** is particularly designed for reconstructing mutation histories of
tumours based on mutation profiles obtained from single-cell exome sequencing experiments, but is in 
principle applicable to any type of (noisy) mutation profiles for which the
infinite sites assumption can be made.

## Availability
---------------

**SCITE** is freely available under a GPL3 license at https://gitlab.com/jahnka/SCITE

##    How to run **SCITE**
--------------------------



### Mac OS X

To compile the C/C++ program, open a terminal and go to the folder containing the source files, and type

	clang++ *.cpp -o scite

This writes a file named `scite`. With older compiler versions you may need to use the option `-std=c++11`.

Assuming the sample data file dataHou18.csv is located in the same folder, `scite` can then be executed as follows

	./scite -i dataHou18.csv -n 18 -m 58 -r 1 -l 900000 -fd 6.04e-5 -ad 0.21545 0.21545 -cc 1.299164e-05

This call computes the ML tree(s) for the given dataset and parameter settings. See below for other program options.

### Linux/Unix

To compile the C/C++ program, open a terminal and go to the folder containing the source files, and type

	g++ *.cpp -o scite
	
This writes a file named `scite`. With older compiler versions you may need to use the option `-std=c++11`.

Assuming the sample data file dataHou18.csv is located in the same folder, `scite` can then be executed as follows

	./scite -i dataHou18.csv -n 18 -m 58 -r 1 -l 900000 -fd 6.04e-5 -ad 0.21545 0.21545 -cc 1.299164e-05

This call computes the ML tree(s) for the given dataset and parameter settings. See below for other program options.

##  Input Files
---------------


### 1. Mutation Matrix


Each column specifies the mutation profile of a single cell, and each row
represents one mutation. Columns are separated by a white space character.

#### (a) Only absence/presence of mutation is distinguished
The entry at position [i,j] should be

* 0 if mutation i is not observed in cell j,
* 1 if mutation i is observed in cell j, or
* 3 if the data point is missing

The sample datasets dataNavin.csv and dataXu.csv have this format.	
#### (b) Heterozygous and homozygous mutations distinguished
The entry at position [i,j] should be

* 0 if mutation i is not observed in cell j,
* 1 if heterozygous mutation i is observed in cell j
* 2 if homozygous mutation i is observed in cell j
* 3 if the data point is missing

The sample datasets dataHou18.csv and dataHou78.csv have this format

### 2. Mutation names (optional)


A list specifying the names of the mutations, e.g. the name of the gene in which
the mutation occurs. For the sample datasets provided here, these files have the extension *.geneNames*. If no such file is specified, the mutations are numbered from 1 to n.

### 3. The true tree (optional)


If the true mutation tree is known (as in simulation experiments), the file with
the true tree (in GraphViz format) can be specified to compare the
predicted trees internally with the true tree.

##  Output Files
----------------

### 1. ML/MAP trees

ML/MAP trees are written to files in GraphViz and Newick format. Files are numbered
consecutively (e.g. dataHou18_ml1.gv, dataHou18_ml1.newick, dataHou18_ml2.gv, dataHou18_ml2.newick, ...). The base name of the output file is derived from the name of the input file (unless a different name is specified via `-o <filename>`).

### 2. Samples from the posterior distribution (optional)

When the `-p <INT>` option is set, **SCITE** samples from the posterior distribution, and writes the sampled trees (in parent vector format) together with their scores and learned error rates to a single file (one sample per line). The name of the output file is derived from the input file name using the ending *.sample*.



## Parameters
-------------

*	`-i <filename>`     Replace \<filename\> with the file containing the mutation matrix

*	`-n <INT>`  Replace \<INT\> with the number of mutations (rows) in the dataset.

* `-m <INT>`  Replace \<INT\> with the  number of cells (columns) in the dataset.

* 	`-r <INT>`  Set \<INT\> to the desired number of repetitions of the MCMC.

* 	`-l <INT>`  Set \<INT\> to the desired chain length of each MCMC repetition



##### In case only absence/presence of a mutation is distinguished

*	`-fd <DOUBLE>` Set \<DOUBLE\> to the estimated false positive rate (false discoveries) of the sequencing experiment.

*	`-ad <DOUBLE>` Set \<DOUBLE\> to the estimated false negative rate (allelic dropout) of the sequencing experiment.

##### In case heterozygous and homozygous mutations are distinguished

*	`-fd <DOUBLE>` Set \<DOUBLE\> to the estimated false positive rate (false calling of heterozygous mutation) of the sequencing. experiment

*	`-ad <DOUBLE> <DOUBLE>` Set the first \<DOUBLE\> to the estimated rate of missed heterozygous mutations in the sequencing experiment and set the second \<DOUBLE\> to the estimated rate of heterozygous mutations called as homozygous mutations (dropout of the wildtype allele).

*	`-cc <DOUBLE>` Set \<DOUBLE\> to the estimated rate of non-mutated sites called as homozygous mutations.   



## Optional parameters
----------------------

#### Program variants


*	`-s` Setting this option causes the sample attachment points (i. e. where the cells would attach to the tree) to be marginalized out.

*	`-p <INT>` When setting this option, **SCITE** samples from the posterior distribution, and writes the trees to a file using the parent vector format. The value of \<INT\> specifies how dense the sampling is. The name of the output file is derived from the input file name using the ending *.sample*.
To make sure that **SCITE** samples from the posterior distribution `-p <INT>` needs to be combined with `-s` (single cell attachment to the tree is marginalized out) and `-g 1` (gamma is set to 1).


*	`-transpose`	This changes the tree representation from mutation tree to rooted binary leaf-labelled tree, where the samples are the leaf labels and mutations are placed on the edges. Using this option can decrease the search space size when there are more mutations than samples. This only works for ML trees.



#### Error learning

*	`-e <double>`   Invokes the learning of error rate beta. Set \<double\> to a value between zero and one to specify the probability to chose the move for changing the error rate in the MCMC.

*	`-x <double>`   Scaling of the known error rate for the MH jump (default is 10).

*	`-sd <double>`  Prior standard deviation for AD error rate (default is 0.1). For the mean, the error rate specified after `-ad` is used.






#### Output

*	`-o <filename>`   Replace \<filename\> with the desired base of the output file to overwrite the default output file names.

*	`-names <filename>` Replace \<filename\> with a file listing the mutation names. By default the mutations are numbered from 1 to n by order of appearance in the inputfile.

*	`-a` When setting this option, **SCITE** adds the individual cells as additional nodes (leafs) to the reported trees. Cells are attached where they fit best in the given tree (with respect to the error rates). By default, only the mutation tree is reported.

*	`-max_treelist_size <INT>`	 This limits the number of co-optimal trees written to output files to the specified <INT> value.

#### Other parameters

*	`-g <DOUBLE>` For ML/MAP computation only: Set \<DOUBLE\> to the desired value of gamma (gamma > 1: more local exploration, possibly local optimum; gamma < 1: easier to explore the space, but less deeply). The default value of gamma is 1 which is necessary for the MCMC chain to converge to the posterior distribution.

*	`-seed <INT>`   Replace \<INT\> with a positive integer to be used as a fixed seed for the random number generator.

*	`-no_tree_list`		This turns off the collection of optimal trees. It can be used to speed up the program when only sampling from the posterior distribution.

*	`-t <filename>`  Replace \<filename\> with a file containing the true tree in GraphViz format.

*	`-move_probs <double> <double> <double>`   Changes the default probabilities for the three MCMC moves to the spedified values. The first move is *prune and re-attach*, the second is *swap node labels*, the third is *swap subtrees*. The default values are (0.55, 0.4, 0.05).

When combined with `-transpose` there are only two move types, *prune and re-attach* and *swap leaf labels* with default probabilities (0.4, 0.6).  Therefore the parameter format for changing these move probabilities is `-move_probs <double> <double>`.

