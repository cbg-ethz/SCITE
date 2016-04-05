/*
 * mcmc.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#ifndef MCMC_H
#define MCMC_H

string runMCMCbeta(vector<struct treeBeta>& bestTrees, double* errorRates, int noOfReps, int noOfLoops, double gamma, vector<double> moveProbs, int n, int m, int** dataMatrix, char scoreType, int* trueParentVec, int step, bool sample, double chi, double priorSd, bool useTreeList, char treeType);
double logBetaPDF(double x, double bpriora, double bpriorb);
double proposeNewBeta(double currBeta, double jumpSd);
double sampleNormal(double mean, double sd);
string sampleFromPosterior(double currTreeLogScore, int n, int* currTreeParentVec, double betaProb, double currBeta, double currScore);
int updateMinDistToTrueTree(int* trueParentVec, int* currTreeParentVec, int length, int minDistToTrueTree, int currScore, int bestScore);
int getSimpleDistance(int* trueVector, int* predVector, int length);


#endif
