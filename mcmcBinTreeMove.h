/*
 * mcmcBinTreeMove.h
 *
 *  Created on: Mar 15, 2016
 *      Author: jahnka
 */

#ifndef MCMCBINTREEMOVE_H_
#define MCMCBINTREEMOVE_H_

int* proposeNextBinTree(std::vector<double> moveProbs, int m, int* currTreeParVec, bool** currTreeAncMatrix);
int pickNodeToMove(int* currTreeParentVec, int parentVectorLength);
int getSibling(int v, int* currTreeParVec, std::vector<std::vector<int> > &childLists);



#endif /* MCMCBINTREEMOVE_H_ */
