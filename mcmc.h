/*
 * mcmc.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#ifndef MCMC_H
#define MCMC_H

double runMCMC(std::vector<int*>& bestTrees, double* errorRates, int noOfReps, int noOfLoops, double beta, std::vector<double> moveProbs, int n, int m, int** dataMatrix, char scoreType, int* trueParentVec, std::vector<int*>& sampleTrees, int step, bool sample);

int* getNewParentVecFast(int* currTreeParentVec, int nodeToMove, int newParent, int n);
int* getNewParentVec_SwapFast(int* currTreeParentVec, int first, int second, int n);
int* proposeNextTreeFast(std::vector<double> moveProbs, int n, bool** currTreeAncMatrix, int* currTreeParentVec, int& nbhcorrection);
int proposeNextTree(std::vector<double> moveProbs, int n, bool** currTreeAncMatrix, int* currTreeParentVec, int*& propTreeParVec, bool**& propTreeAncMatrix);
std::vector<int> getDescendants(bool** ancMatrix, int node, int n);
std::vector<int> getNonDescendants(bool**& ancMatrix, int node, int n);
int* reorderToStartWithDescendant(int* nodestoswap, bool** currTreeAncMatrix);
int* getNewParentVec_Swap(int* currTreeParentVec, int first, int second, int n, int *propTreeParVec);
bool** getNewAncMatrix_Swap(bool** currTreeAncMatrix, int first, int second, int n, bool** propTreeAncMatrix);
int* getNewParentVec(int* currTreeParentVec, int nodeToMove, int newParent, int n, int *propTreeParVec);
bool** getNewAncMatrix(bool** currTreeAncMatrix, int newParent, std::vector<int> descendants, std::vector<int> possibleParents, int n, bool** propTreeAncMatrix);
int choseParent(std::vector<int> &possibleParents, int root);
int getSimpleDistance(int* trueVector, int* predVector, int n);


#endif
