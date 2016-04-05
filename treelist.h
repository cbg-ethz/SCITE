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

struct treeBeta
{
	int* tree;
    double beta;
};

void updateTreeList(std::vector<struct treeBeta>& bestTrees, int* currTreeParentVec, int n, double currScore, double bestScore, double beta);
void resetTreeList(std::vector<struct treeBeta>& bestTrees, int* newBestTree, int n, double beta);
void emptyVectorFast(std::vector<struct treeBeta>& optimalTrees, int n);
void emptyTreeList(std::vector<int*>& optimalTrees, int n);
struct treeBeta createNewTreeListElement(int* tree, int n, double beta);
bool isDuplicateTreeFast(std::vector<struct treeBeta> &optimalTrees, int* newTree, int n);

#endif
