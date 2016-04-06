/*
 * mcmc.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include "matrices.h"
#include "treelist.h"
#include "trees.h"
#include "mcmc.h"
#include "scoreTree.h"
#include "scoreBinTree.h"
#include "rand.h"
#include "limits.h"
#include "output.h"
#include "mcmcBinTreeMove.h"
#include "mcmcTreeMove.h"
using namespace std;
//default_random_engine generator;

/* Run the MCMC to find the best scoring trees, returns best score, inserts all best trees to bestTrees */
/* noOfReps = how often to repeat the MCMC */
/* noOfLoops = how long to run the chain */
/* gamma = scaling factor for the scores to speed up convergence */

unsigned int optCount = 0;    // number of steps spent in optimal tree after burn-in phase
double burnInPhase = 0.25;    // first quarter of steps are burn in phase



/* This runs the MCMC for learning the tree and beta, or only the tree with a fixed beta, it samples from the posterior and/or records the optimal trees/beta */
std::string runMCMCbeta(vector<struct treeBeta>& bestTrees, double* errorRates, int noOfReps, int noOfLoops, double gamma, vector<double> moveProbs, int n, int m, int** dataMatrix, char scoreType, int* trueParentVec, int step, bool sample, double chi, double priorSd, bool useTreeList, char treeType){


	unsigned int optStatesAfterBurnIn = 0;
	int burnIn = noOfLoops*burnInPhase;
	int parentVectorSize = n;
	if(treeType=='t'){parentVectorSize = (2*m)-2;}                     // transposed case: binary tree, m leafs and m-1 inner nodes, root has no parent
	double betaPriorMean = errorRates[1] + errorRates[2];             // AD1 + AD2 the prior mean for AD error rate
	double betaPriorSd   = priorSd;                                     //  prior sd for AD error rate
	double bpriora = ((1-betaPriorMean)*betaPriorMean*betaPriorMean/(betaPriorSd*betaPriorSd)) - betaPriorMean;     // <-10.13585344 turn the mean and sd into parameters of the beta distribution
	double bpriorb = bpriora*((1/betaPriorMean)-1);           //<-13.38666556
	double jumpSd = betaPriorSd/chi;                          // chi: scaling of the known error rate for the MH jump; resulting jump sd
	//cout << "betaPriorMean: " << betaPriorMean << "\n";
	//cout << "betaPriorSd:   " << betaPriorSd << "\n";
	//cout << "bpriora:       " << bpriora << "\n";
	//cout << "bpriorb:       " << bpriorb << "\n";
	//printLogScores(logScores);

	int minDistToTrueTree = INT_MAX;             // smallest distance between an optimal tree and the true (if given)
	double bestTreeLogScore = -DBL_MAX;          // log score of T in best (T,beta)
	double bestScore = -DBL_MAX;                 // log score of best combination (T, beta)
	double bestBeta = betaPriorMean;
	stringstream sampleOutput;

	for(int r=0; r<noOfReps; r++){   // repeat the MCMC, start over with random tree each time, only best score and list of best trees is kept between repetitions

		//cout << "MCMC repetition " << r << "\n";
		int*   currTreeParentVec;
		if(treeType=='m'){currTreeParentVec = getRandParentVec(parentVectorSize);}                                     // start MCMC with random tree
		else{             currTreeParentVec = getRandomBinaryTree(m);}                                                 // transposed case: random binary tree

		bool** currTreeAncMatrix =  parentVector2ancMatrix(currTreeParentVec,parentVectorSize);
		double** currLogScores = getLogScores(errorRates[0], errorRates[1], errorRates[2], errorRates[3]);           // compute logScores of conditional probabilities
		double currBeta = betaPriorMean;                                                                                  // the current AD rate
		double currTreeLogScore;
		if(treeType=='m'){ currTreeLogScore = scoreTreeAccurate( n, m, currLogScores, dataMatrix, scoreType, currTreeParentVec);}
		else{              currTreeLogScore = getBinTreeScore(dataMatrix, n, m, currLogScores, currTreeParentVec);}
		double currBetaLogScore = (moveProbs[0]==0) ? 0.0 : logBetaPDF(currBeta, bpriora, bpriorb);                     // zero if beta is fixed
		double currScore = currTreeLogScore+currBetaLogScore;                                                         // combined score of current tree and current beta

		for(int it=0; it<noOfLoops; it++){                                     // run the iterations of the MCMC
        	if(it % 100000 == 0){ cout << "At mcmc repetition " << r+1 << "/" << noOfReps << ", step " << it << ": best tree score " << bestTreeLogScore << " and best beta " << bestBeta << " and best overall score " << bestScore << "\n";}

        	bool moveAccepted = false;                                           // Is the MCMC move accepted?
        	bool moveChangesBeta = changeBeta(moveProbs[0]);                     // true if this move changes beta, not the tree

        	if(moveChangesBeta){                                                                // new beta is proposed, log scores change tree is copy of current tree
        		double propBeta = proposeNewBeta(currBeta, jumpSd);
        		double** propLogScores = deepCopy_doubleMatrix(currLogScores, 4, 2);
        		updateLogScores(propLogScores, propBeta);
        		double propBetaLogScore = logBetaPDF(propBeta, bpriora, bpriorb);
        		double propTreeLogScore;
        		if(treeType=='m'){ propTreeLogScore = scoreTree( n, m, propLogScores, dataMatrix, scoreType, currTreeParentVec, bestTreeLogScore);}   // compute the new tree score for new beta
        		else{              propTreeLogScore = getBinTreeScore(dataMatrix, n, m, propLogScores, currTreeParentVec);}

        		if (sample_0_1() < exp((propTreeLogScore+propBetaLogScore-currTreeLogScore-currBetaLogScore)*gamma)){               // the proposed move is accepted
        			moveAccepted = true;
        			free_doubleMatrix(currLogScores);
        		    currTreeLogScore  = propTreeLogScore;                                       // update score of current tree
        		    currBeta = propBeta;                                                        // the current AD rate
        		    currBetaLogScore = propBetaLogScore;
        		    currScore = currTreeLogScore+currBetaLogScore;                          // combined score of current tree and current beta
        		    currLogScores = propLogScores;
        		}
        		else{
        			delete [] propLogScores[0];
        			delete [] propLogScores;
        		}
        	}
        	else{                                   // move changed tree
        		double nbhcorrection = 1.0;
        		int* propTreeParVec;
        		double propTreeLogScore;
        		if(treeType=='m'){ propTreeParVec = proposeNewTree(moveProbs, n, currTreeAncMatrix, currTreeParentVec, nbhcorrection);              // propose new tree and
        		                   propTreeLogScore = scoreTree( n, m, currLogScores, dataMatrix, scoreType, propTreeParVec, bestTreeLogScore);}    //  get the new tree score
        		else{              propTreeParVec = proposeNextBinTree(moveProbs, m, currTreeParentVec, currTreeAncMatrix);
        		                   propTreeLogScore = getBinTreeScore(dataMatrix, n, m, currLogScores, propTreeParVec);}

        		if (sample_0_1() < nbhcorrection*exp((propTreeLogScore-currTreeLogScore)*gamma)){                    // the proposed tree is accepted
        			moveAccepted = true;
        			free_boolMatrix(currTreeAncMatrix);                                            // discard outdated tree
        			delete[] currTreeParentVec;
        			currTreeAncMatrix = parentVector2ancMatrix(propTreeParVec,parentVectorSize); // update matrix of current tree
        			currTreeParentVec = propTreeParVec;                                         // update parent vector of current tree
        			currTreeLogScore  = propTreeLogScore;                                       // update score of current tree
        			currScore = currTreeLogScore+currBetaLogScore;
        		}
        		else{
        			delete [] propTreeParVec;            // discard proposed tree
        		}
        	}

        	/* If the true tree is given update the smallest distance between a currently best tree and the true tree */
        	if(trueParentVec){
        		minDistToTrueTree = updateMinDistToTrueTree(trueParentVec, currTreeParentVec, parentVectorSize, minDistToTrueTree, currScore, bestScore);
        	}

        	/* If the list of optimal trees is used, update it */
        	if(useTreeList){
        		updateTreeList(bestTrees, currTreeParentVec, parentVectorSize, currScore, bestScore, currBeta);
        	}

        	/* Sample from the posterior if required and past the burn-in phase */
        	if(sample && it>=burnIn && it % step == 0){
        		sampleOutput << sampleFromPosterior(currTreeLogScore, parentVectorSize, currTreeParentVec, moveProbs[0], currBeta, currScore);
        	}

        	/* Update best tree in case we have found a new best one */
        	if(currScore > bestScore){
        		optStatesAfterBurnIn = 0;                    // new opt state found, discard old count
        		bestTreeLogScore = currTreeLogScore;
        		bestScore = currScore;                 // log score of best combination (T, beta)
        		bestBeta = currBeta;
        	}

        	/* Update the number of MCMC steps we spent in an optimal state */
        	if(currScore == bestScore && it>=burnIn){
        		optStatesAfterBurnIn++;
        	}
        }
        delete [] currTreeParentVec;
        free_doubleMatrix(currLogScores);
        free_boolMatrix(currTreeAncMatrix);
	}                                              // last repetition of MCMC done

	unsigned int noStepsAfterBurnin = noOfReps*(noOfLoops-burnIn);
	cout.precision(17);
	cout << "best log score for tree:\t" << bestTreeLogScore <<  "\n";
	cout << "#optimal steps after burn-in:\t" << optStatesAfterBurnIn << "\n";
	cout << "total #steps after burn-in:\t" << noStepsAfterBurnin << "\n";
	cout << "%optimal steps after burn-in:\t" << (1.0*optStatesAfterBurnIn)/noStepsAfterBurnin << "\n";
	if(moveProbs[0]!=0.0){
		cout << "best value for beta:\t" << bestBeta << "\n";
		cout << "best log score for (T, beta):\t" << bestScore << "\n";
	}

	return sampleOutput.str();
}


double logBetaPDF(double x, double bpriora, double bpriorb){
	double logScore = log(tgamma(bpriora+bpriorb))+(bpriora-1)*log(x)+(bpriorb-1)*log(1-x)-log(tgamma(bpriora))-log(tgamma(bpriorb));    // f(x,a,b) = gamma(a+b)/(gamma(a)gamma(b)) * x^(a-1) * (1-x)^(b-1)
	return logScore;
}

/* a new value for the error probability beta is sampled from a normal distribution around the current beta */
double proposeNewBeta(double currBeta, double jumpSd){
	double sampledValue = sampleNormal(0, jumpSd);
	double propBeta = currBeta+sampledValue ;                   //rnorm(1,0,jumpsd)
	if(propBeta < 0){
		propBeta = abs(propBeta);
	}
	if(propBeta > 1){
		propBeta = propBeta - 2*(propBeta-1);
	}
    return propBeta;
}


/* samples a new value for beta from a normal distribution around the current value */
double sampleNormal(double mean, double sd) {
    double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
    double r = u * u + v * v;
    if (r == 0 || r > 1){
    	return sampleNormal(mean, sd);
    }
    double c = sqrt(-2 * log(r) / r);
    double value =  u * c;                       // value times sd and add the mean
    return (value * sd + mean);
}


/* prints out the current tree and beta to sample from the posterior distribution */
string sampleFromPosterior(double currTreeLogScore, int n, int* currTreeParentVec, double betaProb, double currBeta, double currScore){

	std::stringstream content;
	content << currTreeLogScore  << "\t";                 // logscore of current tree
	content << countBranches(currTreeParentVec, n);       // number of branches in current tree
	if(betaProb>0.0){
		content << "\t" << currBeta;                      // current beta
		content << "\t" << currScore;                     // current combined logscore for tree and beta
	}
	content << "\t";
	for(int i=0; i<n; i++){
		content << currTreeParentVec[i] << " ";
	}
	content << "\n";
	return content.str();
}


/* updates the minimum distance between any of the optimal trees and the true tree (if available) */
int updateMinDistToTrueTree(int* trueParentVec, int* currTreeParentVec, int length, int minDistToTrueTree, int currScore, int bestScore){

	int currDistToTrueTree = getSimpleDistance(trueParentVec, currTreeParentVec, length);

	if(currScore >= bestScore){
		return currDistToTrueTree;
	}

	if(currScore == bestScore && currDistToTrueTree < minDistToTrueTree){         // the current tree is closest to the true tree among the current optimal trees
		return currDistToTrueTree;
	}

	return minDistToTrueTree;
}



/* Returns the distance between two trees, where the distance is the number of nodes having different parents in the two trees  */
int getSimpleDistance(int* trueVector, int* predVector, int length){
	int dist = 0;
	for(int i=0; i<length; i++){
		if(trueVector[i]!=predVector[i]){
			dist++;
		}
	}
	return dist;
}
