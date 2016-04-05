/*
 * mcmcBinTreeMove.cpp
 *
 *  Created on: Mar 14, 2016
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
#include "matrices.h"
#include "trees.h"
#include "rand.h"
#include "mcmcBinTreeMove.h"
#include "output.h"
using namespace std;


/* proposes a new binary tree by a single move from the current binary tree based on the move probabilities */
/* the old tree is kept as currTree, the new one is stored as propTreeParVec */
int* proposeNextBinTree(std::vector<double> moveProbs, int m, int* currTreeParVec, bool** currTreeAncMatrix){

	int movetype = sampleRandomMove(moveProbs);      // pick the move type according to move probabilities
	int parVecLength = (2*m)-2;               // 2m-1 nodes, but the root has no parent

	//cout << "move prob 0: " << moveProbs[0] << "\n";
	//cout << "move prob 1: " << moveProbs[1] << "\n";
	//cout << "move prob 2: " << moveProbs[2] << "\n";
	vector<vector<int> >childLists = getChildListFromParentVector(currTreeParVec, parVecLength);
	int* propTreeParVec  = deepCopy_intArray(currTreeParVec, parVecLength);

	if(movetype==1){                                                       /* type 1: prune and re-attach */
		//cout << "move type is prune and re-attach in binary tree\n";
		int v = pickNodeToMove(currTreeParVec, parVecLength);
		int p = currTreeParVec[v];
		int sib = getSibling(v, currTreeParVec, childLists);             // get the sibling of node v and attach it to the
		propTreeParVec[sib] = currTreeParVec[p];                         // grandparent of v, as the parent of v is moved along with v

		std::vector<int> possibleSiblings = getNonDescendants(currTreeAncMatrix, p, parVecLength);    // get list of possible new siblings of v

		if(possibleSiblings.size()==0){
			cerr << "Error: No new sibling found for node " << v << " for move type 1 in binary tree.\n"; // Should never occur. Something wrong with the tree.
			printGraphVizFile(currTreeParVec, parVecLength);
		}

		int newSibling = possibleSiblings[pickRandomNumber(possibleSiblings.size())]; // pick a new sibling from remaining tree (root can not be a sibling)
		propTreeParVec[newSibling] = p;                                               // make the new sibling a child of v's parent
		propTreeParVec[p] = currTreeParVec[newSibling];                            // make the parent of v the child of the new sibling's former parent
	}
    else{                                                                 /* type 2: swap two node labels  */
    	//cout << "move type is swap node labels in binary tree\n";
    	int v =  rand() % m;                                            // get random leaf to swap (only the first m nodes are leafs)
    	int w =  rand() % m;                                            // get second random leaf to swap
    	propTreeParVec[v] = currTreeParVec[w];                         // and just swap parents
    	propTreeParVec[w] = currTreeParVec[v];
    }
    return propTreeParVec;
}


/* returns a node where the prune and re-attach step starts */
int pickNodeToMove(int* currTreeParentVec, int parentVectorLength){
	bool validNode = false;
	int rootId = parentVectorLength;
	while(!validNode){
		int v = pickRandomNumber(parentVectorLength);   // pick a node for the prune and re-attach step;
		if(currTreeParentVec[v]!=rootId){               // it has to be a node whose parent is not the root, as node and parent are moved together
			return v;
		}                                      // for a binary tree with more than two leafs this can not be an infinite loop
	}
}


/* returns the (unique) sibling of node v. The sibling has to exist because tree is binary */
int getSibling(int v, int* currTreeParVec, vector<vector<int> > &childLists){

	if(childLists.at(currTreeParVec[v]).at(0) != v){
		return childLists.at(currTreeParVec[v]).at(0);
	}
	else{
		return childLists.at(currTreeParVec[v]).at(1);
	}
}

