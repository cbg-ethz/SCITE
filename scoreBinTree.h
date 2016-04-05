/*
 * scoreBinTree.h
 *
 *  Created on: Mar 15, 2016
 *      Author: jahnka
 */

#ifndef SCOREBINTREE_H_
#define SCOREBINTREE_H_

double* getBinSubtreeScore(bool state, int* bft, std::vector<std::vector<int> > &childLists, int mut, int nodeCount, int m, int** obsMutProfiles, double ** logScores);
double getBinTreeMutScore(int* bft, std::vector<std::vector<int> > &childLists, int mut, int nodeCount, int m, int** obsMutProfiles, double ** logScores);
double getBinTreeScore(int** obsMutProfiles, int n, int m, double ** logScores, int* parent);



#endif /* SCOREBINTREE_H_ */
