/*
 * treelist.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */


#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <math.h>
#include <queue>
#include "matrices.h"
#include "treelist.h"
#include "rand.h"

using namespace std;


void updateTreeList(vector<struct treeBeta>& bestTrees, int* currTreeParentVec, int n, double currScore, double bestScore, double beta){

	if(currScore > bestScore){
		//cout << "tree list of size " << bestTrees.size() << " emptied\n";
		resetTreeList(bestTrees, currTreeParentVec, n, beta);                              // empty the list of best trees and insert current tree

	}
	else if (currScore == bestScore){
		if(!isDuplicateTreeFast(bestTrees, currTreeParentVec, n)){               // if the same tree was not previously found
			treeBeta newElem = createNewTreeListElement(currTreeParentVec, n, beta);
			bestTrees.push_back(newElem);        // add it to list
		}
	}
}


/* removes all elements from the vector and inserts the new best tree */
void resetTreeList(vector<struct treeBeta>& bestTrees, int* newBestTree, int n, double beta){
	emptyVectorFast(bestTrees, n);                                         // empty the list of best trees
	treeBeta newElem = createNewTreeListElement(newBestTree, n, beta);
	bestTrees.push_back(newElem);                // current tree is now the only best tree
}


/* removes all elements from the vector */
void emptyVectorFast(std::vector<struct treeBeta>& optimalTrees, int n){
    for(int i=0; i<optimalTrees.size(); i++){
    	delete [] optimalTrees[i].tree;
	}
    optimalTrees.clear();
}

/* removes all elements from the vector */
void emptyTreeList(std::vector<int*>& optimalTrees, int n){
    for(int i=0; i<optimalTrees.size(); i++){
    	delete [] optimalTrees[i];
	}
    optimalTrees.clear();
}

/* creates a new tree/beta combination */
struct treeBeta createNewTreeListElement(int* tree, int n, double beta){
	treeBeta newElem;
	newElem.tree = deepCopy_intArray(tree, n);
	newElem.beta = beta;
	return newElem;
}

/* returns true if the same tree was found before */
bool isDuplicateTreeFast(std::vector<struct treeBeta> &optimalTrees, int* newTree, int n){
    for(int k=0; k<optimalTrees.size(); k++){
      bool same = true;
      for(int i=0; i<n; i++){
    	  if(newTree[i] != optimalTrees[k].tree[i]){
              same = false;
              break;
          }
      }
      if(same == true){
        return true;
      }
    }
    return false;
}







///* returns true if the tree has a branching point */
//bool hasBranching(int* parents, int n){
//	bool* isParent = new bool[n+1];
//	for(int i=0; i<n; i++){
//		isParent[i]=false;
//	}
//	for(int i=0; i<n; i++){
//		if(isParent[parents[i]]==true){
//			return true;
//		}
//		isParent[parents[i]] = true;
//	}
//	return false;
//}
//
//void foundBranchingTree(std::vector<bool**> treeList, int n){
//	bool isBranching = false;
//	vector<int> branchingTrees;
//	for(int i=0; i<treeList.size(); i++){
//		int* parVector = ancMatrixToParVector(treeList[i], n);
//		if(hasBranching(parVector, n)){
//			branchingTrees.push_back(i);
//			isBranching = true;
//		}
//	}
//	if(isBranching==true){
//		cout << branchingTrees.size() << " out of " << treeList.size() << "trees are branching: ";
//		for (int i=0; i<branchingTrees.size(); i++){
//		    cout << branchingTrees[i] << " ";
//		}
//		cout << "\n";
//	}
//	else{
//		printf("no branching trees in list\n");
//	}
//}

//bool isDuplicateTree(std::vector<bool**> &optimalTrees, bool** newTree, int n){
//    for(int k=0; k<optimalTrees.size(); k++){
//      bool same = true;
//      for(int i=0; i<n && same==true; i++){
//        for(int j=0; j<n && same==true; j++){
//            if(newTree[i][j] != optimalTrees[k][i][j]){
//              same = false;
//            }
//        }
//      }
//      if(same == true){
//        return true;
//      }
//    }
//    return false;
//}

//double updateListOfBestTrees(double currScore, double bestScore, bool**& propAncMatrix, std::vector<bool**> &optimalTrees, int n){
//
//    if(currScore > bestScore){                                    // a tree with better score is found
//        bestScore = currScore;
//        emptyVector(optimalTrees, n);     // empty outdated list of best trees
//        optimalTrees.push_back(propAncMatrix);	               // add new optimal tree to list
//  	}
//  	else if(currScore == bestScore){                            // another instance with current best score is found
//  		  optimalTrees.push_back(propAncMatrix);	                    // add the new optimal tree to the list
//  	}
//    else{
//        free_boolMatrix(propAncMatrix);
//    }
//	  return bestScore;
//    return 0;
//}

///* empty vector */
//void emptyVector(std::vector<bool**> &optimalTrees, int n){
//    for(int i=optimalTrees.size()-1; i>=0; i--){
//    	free_boolMatrix(optimalTrees.at(i));
//    	optimalTrees.at(i) = NULL;
//	  }
//	  optimalTrees.clear();
//}
