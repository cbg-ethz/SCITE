/*
 * findBestTrees_noR.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include "matrices.h"
#include "treelist.h"
#include "trees.h"
#include "output.h"
#include "mcmc.h"
#include "rand.h"
#include "scoreTree.h"

using namespace std;

int** getDataMatrix(int n, int m, string fileName);
double* getErrorRatesArray(double fd, double ad1, double ad2, double cc);
int readParameters(int argc, char* argv[]);
string getOutputFilePrefix(string fileName, string outFile);
string getFileName(string prefix, string ending);
string getFileName2(int i, string prefix, string ending);
vector<string> getGeneNames(string fileName, int nOrig);
vector<double> setMoveProbs();
int* getParentVectorFromGVfile(string fileName, int n);
int getMinDist(int* trueVector, std::vector<bool**> optimalTrees, int n);
void printGeneFrequencies(int** dataMatrix, int n, int m, vector<string> geneNames);

static double move1_prob = 0.1;
static double move2_prob = 0.65;
static double move3_prob = 0.25;
double errorRateMove = 0.0;
double chi = 10;
double priorSd = 0.1;
string fileName;      // data file
string outFile;       // the name of the outputfile, only the prefix before the dot
int n;                // number of genes
int m;                // number of samples
bool sumScore = false;
int rep;            // number of repetitions of the MCMC
int loops;          // number of loops within a MCMC
double gamma;
double fd;          // rate of false discoveries (false positives 0->1)
double ad1;          // rate of allelic dropout (false negatives 1->0)
double ad2 = 0.0;         // rate of allelic dropout (2->1)
double cc = 0.0;          // rate of falsely discovered homozygous mutations (0->2)
bool sample = false;
int sampleStep;
bool useGeneNames = false;        // use gene names in tree plotting
string geneNameFile;              // file where the gene names are listed.
bool trueTreeComp = false;      // set to true if true tree is given as parameter for comparison
string trueTreeFileName;      // optional true tree
bool attachSamples = false;       // attach samples to the tree
bool useFixedSeed = false;      // use a predefined seed for the random number generator
unsigned int fixedSeed = 1;

int main(int argc, char* argv[])
{

	/**  read parameters and data file  **/

	readParameters(argc, argv);
	int** dataMatrix = getDataMatrix(n, m, fileName);
	vector<double> moveProbs = setMoveProbs();
	double* errorRates = getErrorRatesArray(fd, ad1, ad2, cc);
	char scoreType = 'm';
	if(sumScore){ scoreType = 's';}


	/****************   begin timing  *********************/
		clock_t begin=clock();
	/****************************************************/

	std::vector<int*> optimalTrees;          // list of all optimal trees (as ancestor matrices)
	double bestScore = -DBL_MAX;                          // score of the optimal trees


	/* get the true parent vector from GraphViz file if available (for simulated data only)  */
	int* trueParentVec = NULL;
	if(trueTreeComp==true){
		trueParentVec = getParentVectorFromGVfile(trueTreeFileName, n);
	}

	std::vector<int*> sampleTrees;                    // list where tree samples are stored, if sampling based on posterior distribution is needed

	/* initialize the random number generator, either with a user defined seed, or a random number */
	if(useFixedSeed){
		srand(fixedSeed);
	}
	else{
		initRand();                                  // initialize random number generator
	}

	/***  Find best scoring trees by MCMC  ***/
	//cout << "running MCMC now...\n";
	if(errorRateMove == 0.0){
		bestScore = runMCMC(optimalTrees, errorRates, rep, loops, gamma, moveProbs, n, m, dataMatrix, scoreType, trueParentVec, sampleTrees, sampleStep, sample);
	}
	else{
		bestScore = runMCMCbeta(optimalTrees, errorRates, rep, loops, gamma, moveProbs, n, m, dataMatrix, scoreType, trueParentVec, sampleTrees, sampleStep, sample, chi, priorSd);
	}


	//cout << "best score: " << bestScore << "\n";

    /***  output results  ***/

	int trueN = n;

	string prefix = getOutputFilePrefix(fileName, outFile);
	for(int i=0; i<optimalTrees.size(); i++){
		string outputFile = getFileName2(i, prefix, ".newick");
		//cout << outputFile << "\n";
		int* parentVector = optimalTrees.at(i);
		vector<vector<int> > childLists = getChildListFromParentVector(parentVector, trueN);
		stringstream newick;
		newick << getNewickCode(childLists, trueN) << "\n";
		writeToFile(newick.str(), outputFile);
		outputFile = getFileName2(i, prefix, ".gv");
		double** logScores = getLogScores(fd, ad1, ad2, cc);
		bool** ancMatrix = parentVector2ancMatrix(parentVector, trueN);
		string output = getGraphVizFileContentNames(parentVector, trueN, getGeneNames(geneNameFile, n), attachSamples, ancMatrix, m, logScores, dataMatrix);
		free_boolMatrix(ancMatrix);
		delete [] logScores[0];
		delete [] logScores;
		writeToFile(output, outputFile);
	}

	//printSampleTrees(sampleTrees, n, getFileName(getOutputFilePrefix(fileName, outFile), ".sample"));


	emptyVectorFast(optimalTrees, n);
	delete [] errorRates;
	/****************   end timing  *********************/
  		clock_t end=clock();
  		double diffticks=end-begin;
  		double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
  		//cout << "Time elapsed: " << diffms << " ms"<< endl;
  	/****************************************************/

  	free_intMatrix(dataMatrix);
}




void printGeneFrequencies(int** dataMatrix, int n, int m, vector<string> geneNames){
	for(int i=0; i<n; i++){
		int freq = 0;
		for(int j=0; j<m; j++){
			if(dataMatrix[j][i]==1 || dataMatrix[j][i]==2){
				freq++;
			}
		}
		cout << freq << "\t" << geneNames.at(i) << "\n";
	}
}



int* getParentVectorFromGVfile(string fileName, int n){
	int* parentVector = new int[n];
	std::vector<std::string> lines;
	std::ifstream file(fileName.c_str());
	std::string line;
	while ( std::getline(file, line) ) {
	    if ( !line.empty() )
	        lines.push_back(line);
	}
	for(int i=0; i < lines.size(); i++){

		std::size_t found = lines[i].find(" -> ");
		if (found!=std::string::npos){
			int parent = atoi(lines[i].substr(0, found).c_str());
			int child = atoi(lines[i].substr(found+3).c_str());
			parentVector[child-1] = parent-1;
	   }
	}
	return parentVector;
}



int getMinDist(int* trueVector, std::vector<bool**> optimalTrees, int n){
	int minDist = n+1;
	for(int i=0; i<optimalTrees.size(); i++){
		int dist = getSimpleDistance(trueVector, ancMatrixToParVector(optimalTrees.at(i), n), n);
		minDist = min(minDist, dist);
	}
	return minDist;
}


string getOutputFilePrefix(string fileName, string outFile){
	if(outFile.empty()){
		int lastIndex = fileName.find_last_of(".");
		return fileName.substr(0, lastIndex);
	}
	return outFile;
}


string getFileName(string prefix, string ending){
	stringstream fileName;
	fileName << prefix << ending;
	return fileName.str();
}

string getFileName2(int i, string prefix, string ending){
	stringstream fileName;
	fileName << prefix << "_ml" << i << ending;
	return fileName.str();
}

int readParameters(int argc, char* argv[]){
	for (int i = 1; i < argc; ++i) {

		if (strcmp(argv[i], "-i") == 0) {
			if (i + 1 < argc) { fileName = argv[++i];}
		} else if (strcmp(argv[i], "-t") == 0) {
			if (i + 1 < argc) {
				trueTreeFileName = argv[++i];
				trueTreeComp = true;
			}
		} else if(strcmp(argv[i], "-o")==0) {
			if (i + 1 < argc) { outFile = argv[++i];}
		} else if(strcmp(argv[i], "-n")==0) {
			if (i + 1 < argc) { n = atoi(argv[++i]);}
		} else if(strcmp(argv[i], "-m")==0) {
			if (i + 1 < argc) { m = atoi(argv[++i]);}
		} else if(strcmp(argv[i], "-r") == 0) {
			if (i + 1 < argc) { rep = atoi(argv[++i]);}
		} else if(strcmp(argv[i], "-l")==0) {
			if (i + 1 < argc) { loops = atoi(argv[++i]);}
		} else if(strcmp(argv[i], "-g")==0) {
			if (i + 1 < argc) { gamma = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-fd")==0) {
			if (i + 1 < argc) { fd = atof(argv[++i]);}
		} else if(strcmp(argv[i],"-ad")==0) {
			if (i + 1 < argc) { ad1 = atof(argv[++i]);}
			if (i + 1 < argc){
				string next = argv[i+1];
				if(next.compare(0, 1, "-") != 0){
					ad2 = atof(argv[++i]);
				}
			}
		} else if(strcmp(argv[i],"-cc")==0) {
			if (i + 1 < argc) { cc = atof(argv[++i]);}
		} else if(strcmp(argv[i],"-e")==0) {
					if (i + 1 < argc) { errorRateMove = atof(argv[++i]);}
		} else if(strcmp(argv[i],"-x")==0) {
							if (i + 1 < argc) { chi = atof(argv[++i]);}
		} else if(strcmp(argv[i],"-sd")==0) {
									if (i + 1 < argc) { priorSd = atof(argv[++i]);}
		} else if (strcmp(argv[i], "-a")==0) {
			attachSamples = true;
		} else if(strcmp(argv[i], "-p")==0) {
			if (i + 1 < argc) {
				sampleStep = atoi(argv[++i]);
				sample = true;
			}
		}else if (strcmp(argv[i], "-names")==0) {
			useGeneNames = true;
			if (i + 1 < argc) { geneNameFile = argv[++i];}
		}else if (strcmp(argv[i], "-seed")==0) {
			useFixedSeed = true;
			if (i + 1 < argc) { fixedSeed = atoi(argv[++i]);}
		} else if (strcmp(argv[i],"-s")==0) {
			sumScore = true;
		} else {
			std::cerr << "unknown parameter " << argv[i] << std::endl;
			return 1;
		}
	}
	return 0;
}


vector<double> setMoveProbs(){
	vector<double> moveProbs;
	moveProbs.push_back(move1_prob*(1-errorRateMove));
	moveProbs.push_back(move2_prob*(1-errorRateMove));
	moveProbs.push_back(move3_prob*(1-errorRateMove));
	moveProbs.push_back(errorRateMove);
	return moveProbs;
}


int** getDataMatrix(int n, int m, string fileName){

    int** dataMatrix = init_intMatrix(n, m, -1);

    ifstream in(fileName.c_str());

    if (!in) {
    	cout << "2 Cannot open file " << fileName << "\n";
      cout << fileName;
      cout << "\n";
      return NULL;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            in >> dataMatrix[i][j];
        }
    }

    in.close();
    int** transposedMatrix = transposeMatrix(dataMatrix, n, m);
    free_intMatrix(dataMatrix);
    return transposedMatrix;
}


vector<string> getGeneNames(string fileName, int nOrig){

	vector<string> v;
	ifstream in(fileName.c_str());


	n = nOrig;

	if (!in) {
		//cout << "Cannot open gene names file " << fileName << ", ";
	    //cout << "using ids instead.\n";
	    vector<string> empty;
	    for(int i=0; i<=n; i++){
	    	stringstream id;
	    	id << i+1;
	    	empty.push_back(id.str());
	    }
	    return empty;
	}

	for (int i = 0; i < nOrig; i++) {
		string temp;
	    in >> temp;
	    v.push_back(temp);
	}
	v.push_back("Root"); // the root
	return v;
}


double* getErrorRatesArray(double fd, double ad1, double ad2, double cc){
	double* array = new double[4];
	array[0] = fd;
	array[1] = ad1;
	array[2] = ad2;
	array[3] = cc;
	return array;
}

