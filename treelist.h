/*
 * treelist.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#include <string>
//#include <iostream>
//#include <sstream>

#ifndef TREELIST_H
#define TREELIST_H

bool isDuplicateTreeFast(std::vector<int*> &optimalTrees, int* newTree, int n);
bool isDuplicateTree(std::vector<bool**> &optimalTrees, bool** newTree, int n);
double updateListOfBestTrees(double currScore, double bestScore, bool**& propAncMatrix, std::vector<bool**> &optimalTrees, int n);
void emptyVector(std::vector<bool**> &optimalTrees, int n);
void emptyVectorFast(std::vector<int*> & optimalTrees, int n);
void foundBranchingTree(std::vector<bool**> treeList, int n);
void printParentVectors(std::vector<bool**> optimalTrees, int n, int m, double** logScores, int** dataMatrix);

#endif
