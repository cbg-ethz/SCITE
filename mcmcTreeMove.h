/*
 * mcmcTreeMove.h
 *
 *  Created on: Mar 15, 2016
 *      Author: jahnka
 */

#ifndef MCMCTREEMOVE_H_
#define MCMCTREEMOVE_H_

int* proposeNewTree(std::vector<double> moveProbs, int n, bool** currTreeAncMatrix, int* currTreeParentVec, double& nbhcorrection);
int choseParent(std::vector<int> &possibleParents, int root);
int* getNewParentVecFast(int* currTreeParentVec, int nodeToMove, int newParent, int n);
int* getNewParentVec_SwapFast(int* currTreeParentVec, int first, int second, int n);
int* reorderToStartWithDescendant(int* nodestoswap, bool** currTreeAncMatrix);
int* getNewParentVec_Swap(int* currTreeParentVec, int first, int second, int n, int* propTreeParVec);
bool** getNewAncMatrix_Swap(bool** currTreeAncMatrix, int first, int second, int n, bool** propTreeAncMatrix);
int* getNewParentVec(int* currTreeParentVec, int nodeToMove, int newParent, int n, int *propTreeParVec);
bool** getNewAncMatrix(bool** currTreeAncMatrix, int newParent, std::vector<int> descendants, std::vector<int> possibleParents, int n, bool** propTreeAncMatrix);



#endif /* MCMCTREEMOVE_H_ */
