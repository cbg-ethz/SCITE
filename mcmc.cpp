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
#include <fstream>
#include <sstream>
#include "matrices.h"
#include "treelist.h"
#include "trees.h"
#include "mcmc.h"
#include "scoreTree.h"
#include "rand.h"

using namespace std;


/* Run the MCMC to find the best scoring trees, returns best score, inserts all best trees to bestTrees */
/* noOfReps = how often to repeat the MCMC */
/* noOfLoops = how long to run the chain */
/* gamma = scaling factor for the scores to speed up convergence */



double runMCMC(std::vector<int*>& bestTrees, double* errorRates, int noOfReps, int noOfLoops, double gamma, std::vector<double> moveProbs, int n, int m, int** dataMatrix, char scoreType, int* trueParentVec, std::vector<int*>& sampleTrees, int step, bool sample){

	double ** logScores = getLogScores(errorRates[0], errorRates[1], errorRates[2], errorRates[3]);    // compute logScores of conditional probabilities
	//printLogScores(logScores);
	initRand();                                  // initialize random number generator
	int minDist = INT_MAX;                        // initialize distance to true tree if given
	double bestTreeLogScore = -DBL_MAX;          // initialize best tree score

	for(int r=0; r<noOfReps; r++){   // repeat the MCMC, start over with random tree each time, only best score and list of best trees is kept between repetitions

		//cout << "MCMC repetition " << r << "\n";
		int*   currTreeParentVec = getRandParentVec(n);                                                                // start MCMC with random tree
		bool** currTreeAncMatrix =  parentVector2ancMatrix(currTreeParentVec,n);
		double currTreeLogScore = scoreTreeAccurate( n, m, logScores, dataMatrix, scoreType, currTreeParentVec);       // get score of random tree

        for(int it=0; it<noOfLoops; it++){        // run the iterations of the MCMC
        	//if(it % 100000 == 0){ printf("%d iterations %f score \n", it, currTreeLogScore);}

        	int nbhcorrection = 1;
        	int* propTreeParVec = proposeNextTreeFast(moveProbs, n, currTreeAncMatrix, currTreeParentVec, nbhcorrection); // propose a new tree from neighborhood of current tree
        	double propTreeLogScore = scoreTree( n, m, logScores, dataMatrix, scoreType, propTreeParVec, bestTreeLogScore);   //   compute the score of the new tree

        	if(sample==true && it % step == 0){
        		sampleTrees.push_back(deepCopy_intArray(propTreeParVec, n));
        	}

        	if (sample_0_1() < nbhcorrection*exp((propTreeLogScore-currTreeLogScore)*gamma)){       // the proposed tree is accepted
  			    free_boolMatrix(currTreeAncMatrix);                                            // discard outdated tree
  			    delete[] currTreeParentVec;
    			currTreeAncMatrix = parentVector2ancMatrix(propTreeParVec,n);               // update matrix of current tree
    			currTreeParentVec = propTreeParVec;                                         // update parent vector of current tree
    			currTreeLogScore  = propTreeLogScore;                                       // update score of current tree

    			int currDist = -1;                                                         // current distance to the true tree (if given)
    			if(trueParentVec!=NULL){
    				currDist = getSimpleDistance(trueParentVec, currTreeParentVec, n);
    			}
    			if(currTreeLogScore > bestTreeLogScore){                                 // the current tree has better score than the previously best tree
    			    bestTreeLogScore = currTreeLogScore;                                // update best score

    			    emptyVectorFast(bestTrees, n);                                         // empty the list of best trees (ancestor matrices)
    			    minDist = currDist;
    			    bestTrees.push_back(deepCopy_intArray(currTreeParentVec, n));  // current tree is now the only best tree
    			}
    			else if(currTreeLogScore == bestTreeLogScore){                                 // the current tree is equally good as the best tree so far
    			    if(!isDuplicateTreeFast(bestTrees, currTreeParentVec, n)){               // if the same tree was not previously found,
    			    	if(currDist==minDist || trueParentVec==NULL){
    			    		bestTrees.push_back(deepCopy_intArray(currTreeParentVec, n));     // add it to list
    			    	}
    			    	else if(currDist<minDist && trueParentVec != NULL){
    			    		emptyVectorFast(bestTrees, n);                                          // empty the list of best trees (ancestor matrices)
    			    		minDist = currDist;
    			    		bestTrees.push_back(deepCopy_intArray(currTreeParentVec, n));
    			    	}
    			    }
    			}
  			 }
  			 else{                                     // the proposed tree was not accepted
  			    delete[] propTreeParVec;            // discard proposed tree
  			 }
        }                                          // end of this round of MCMC
        delete [] currTreeParentVec;
        free_boolMatrix(currTreeAncMatrix);
	}                                              // last repetition of MCMC done

	delete [] logScores[0];
	delete [] logScores;
 	return bestTreeLogScore;
}





/* creates the new parent vector after pruning and reattaching subtree */
int* getNewParentVecFast(int* currTreeParentVec, int nodeToMove, int newParent, int n){
	int* propTreeParVec = deepCopy_intArray(currTreeParentVec, n);        // deep copy of parent matrix to keep old one
    propTreeParVec[nodeToMove] = newParent;                       // only the parent of the moved node changes
    return propTreeParVec;
}


/* creates the new parent vector after swapping two nodes */
int* getNewParentVec_SwapFast(int* currTreeParentVec, int first, int second, int n){

	int* propTreeParVec = deepCopy_intArray(currTreeParentVec, n);
    for(int i=0; i<n; i++){           // set vector of proposed parents
      if(propTreeParVec[i] == first && i!=second){
      	propTreeParVec[i] = second;              // update entries with swapped parents
      }
      else if(propTreeParVec[i] == second && i!=first){
      	propTreeParVec[i] = first;
      }
    }
    int temp = propTreeParVec[first];
    propTreeParVec[first] = propTreeParVec[second];  // update parents of swapped nodes
    propTreeParVec[second] = temp;
    if(propTreeParVec[first]==first){propTreeParVec[first]=second;}    // this is needed to ensure that tree is connected, the above fails
    if(propTreeParVec[second]==second){propTreeParVec[second]=first;}  // if first is parent of second, or vice versa
    return propTreeParVec;
}

int* proposeNextTreeFast(std::vector<double> moveProbs, int n, bool** currTreeAncMatrix, int* currTreeParentVec, int& nbhcorrection){

	int* propTreeParVec;
	int movetype = sampleRandomMove(moveProbs);      // pick the move type
	nbhcorrection = 1;                               // reset the neighbourhood correction

	if(movetype==1){       /* prune and re-attach */
		//cout << "move type is prune and reattach\n";
		int nodeToMove = pickRandomNumber(n);   // pick a node to move with its subtree
		std::vector<int> possibleparents = getNonDescendants(currTreeAncMatrix, nodeToMove, n);      // possible attachment points
		int newParent = choseParent(possibleparents, n);                                             // randomly pick a new parent among available nodes, root (n+1) is also possible parent
		propTreeParVec = getNewParentVecFast(currTreeParentVec, nodeToMove, newParent, n);           // create new parent vector
	}
    else if(movetype==2){   /* swap two node labels  */
    	//cout << "move type is swap node labels\n";
    	int* nodestoswap = sampleTwoElementsWithoutReplacement(n);
        propTreeParVec    = getNewParentVec_SwapFast(currTreeParentVec, nodestoswap[0], nodestoswap[1], n);
        delete [] nodestoswap;
    }
    else{       // if(movetype==3){    /*  swap two subtrees  */
	   //cout << "move type is swap subtrees\n";
        int* nodestoswap = sampleTwoElementsWithoutReplacement(n);                    // pick the node that will be swapped
        nodestoswap = reorderToStartWithDescendant(nodestoswap, currTreeAncMatrix);   // make sure we move the descendant first (in case nodes are in same lineage)
        int nodeToMove = nodestoswap[0];                                              // now we move the first node chosen and its descendants
        int nextnodeToMove = nodestoswap[1];                                          // next we need to move the second node chosen and its descendants
        delete [] nodestoswap;

        if(currTreeAncMatrix[nextnodeToMove][nodeToMove]==0){                      // the nodes are in different lineages -- simple case

        	propTreeParVec = deepCopy_intArray(currTreeParentVec, n);         // deep copy of parent matrix to keep old one
        	propTreeParVec[nodeToMove] = currTreeParentVec[nextnodeToMove];        // and exchange the parents of the nodes
        	propTreeParVec[nextnodeToMove] = currTreeParentVec[nodeToMove];
        }
        else{                                                                      // the nodes are in the same lineage -- need to avoid cycles in the tree
        	propTreeParVec = deepCopy_intArray(currTreeParentVec, n);         // deep copy of parent vector to keep old one
        	propTreeParVec[nodeToMove] = currTreeParentVec[nextnodeToMove];        // lower node is attached to the parent of the upper node
        	std::vector<int> descendants     = getDescendants(currTreeAncMatrix, nodeToMove, n);   // all nodes in the subtree of the lower node
        	bool** propTreeAncMatrix = parentVector2ancMatrix(propTreeParVec, n);
        	std::vector<int> nextdescendants = getDescendants(propTreeAncMatrix, nextnodeToMove, n);
        	free_boolMatrix(propTreeAncMatrix);
        	propTreeParVec[nextnodeToMove] = descendants[pickRandomNumber(descendants.size())];  // node to move is attached to a node chosen uniformly from the descendants of the first node
        	nbhcorrection = descendants.size()/nextdescendants.size(); // neighbourhood correction needed for MCMC convergence, but not important for simulated annealing
        }
    }
    return propTreeParVec;
}

int proposeNextTree(std::vector<double> moveProbs, int n, bool** currTreeAncMatrix, int* currTreeParentVec, int*& propTreeParVec, bool**& propTreeAncMatrix){

	int movetype = sampleRandomMove(moveProbs);   // pick the move type
	int nbhcorrection = 1;               // reset the neighbourhood correction
	if(movetype==1){       /* prune and re-attach */
		//cout << "move type is prune and reattach\n";
		int nodeToMove = pickRandomNumber(n);   // pick a node to move with its subtree

		//  descendants are those in row nodeToMove with a 1, possible parents are the others
		std::vector<int> descendants  = getDescendants(currTreeAncMatrix, nodeToMove, n);                   // includes the node itself
		std::vector<int> possibleparents = getNonDescendants(currTreeAncMatrix, nodeToMove, n);             // possible attachment points

		int newParent = choseParent(possibleparents, n);                                    // randomly pick a new parent among available nodes, root (n+1) is also possible parent
		propTreeAncMatrix = getNewAncMatrix(currTreeAncMatrix, newParent, descendants, possibleparents, n, propTreeAncMatrix);   // create new ancestor matrix, TODO: avoid deepCopy of old matrix
		propTreeParVec    = getNewParentVec(currTreeParentVec, nodeToMove, newParent, n, propTreeParVec);                     // create new parent vector
	}
    else if(movetype==2){   /* swap two node labels  */
    //	cout << "move type is swap node labels\n";
    	int* nodestoswap = sampleTwoElementsWithoutReplacement(n);
        propTreeAncMatrix = getNewAncMatrix_Swap(currTreeAncMatrix, nodestoswap[0], nodestoswap[1], n, propTreeAncMatrix);  // create ancestor matrix after swapping
        propTreeParVec    = getNewParentVec_Swap(currTreeParentVec, nodestoswap[0], nodestoswap[1], n, propTreeParVec);
        delete [] nodestoswap;
    }
   else{       // if(movetype==3){    /*  swap two subtrees  */
	  // cout << "move type is swap subtrees\n";
        int* nodestoswap = sampleTwoElementsWithoutReplacement(n);                   // pick the nodes whose incoming edges will be cut and swapped around
        nodestoswap = reorderToStartWithDescendant(nodestoswap, currTreeAncMatrix);      // make sure we move the descendant first (in case nodes are in same lineage)
        int nodeToMove = nodestoswap[0];                                              // now we move the first node chosen and its descendants
        //cout << "picked nodes to prune tree: " << nodestoswap[0]+1 << " and " << nodestoswap[1]+1 << "\n";
        std::vector<int> descendants  = getDescendants(currTreeAncMatrix, nodeToMove, n);   // includes the node itself
        std::vector<int> possibleparents = getNonDescendants(currTreeAncMatrix, nodeToMove, n);  // possible attachment points
        int newParent = currTreeParentVec[nodestoswap[1]];                                           // the new parent is the parent of the other node
        //cout << "new parent: " << newParent+1 << " of " <<  nodestoswap[0]+1 << "\n";
        bool** propTreeAncMatrixTemp = getNewAncMatrix(currTreeAncMatrix, newParent, descendants, possibleparents, n, propTreeAncMatrix);
        propTreeParVec = getNewParentVec(currTreeParentVec, nodeToMove, newParent, n, propTreeParVec);
       // print_intArray(propTreeParVec, n);
        int nextnodeToMove = nodestoswap[1];  // next we need to move the second node chosen and its descendants
        delete [] nodestoswap;
        std::vector<int> nextdescendants     = getDescendants(propTreeAncMatrixTemp, nextnodeToMove, n);  // includes the node itself
        std::vector<int> nextpossibleparents = getNonDescendants(propTreeAncMatrixTemp, nextnodeToMove, n);

        if(currTreeAncMatrix[nextnodeToMove][nodeToMove]==1){             // if the two nodes used to be in the same lineage (first descendant of second) , the second
        	newParent = descendants[pickRandomNumber(descendants.size())];  // node to move is attached to a node chosen uniformly from the descendants of the first node
        	nbhcorrection = descendants.size()/nextdescendants.size(); // neighbourhood correction needed for MCMC convergence, but not important for simulated annealing
        }
        else  // if the nodes used to be in different lineages, just attach to the second node to initial parent of the first node
        {
        	newParent = currTreeParentVec[nodeToMove];
        }

        propTreeAncMatrix = getNewAncMatrix(propTreeAncMatrixTemp, newParent, nextdescendants, nextpossibleparents, n, deepCopy_boolMatrix(propTreeAncMatrixTemp, n, n));
        free_boolMatrix(propTreeAncMatrixTemp);
        propTreeParVec[nextnodeToMove] = newParent; // update the parent vector of the moved node
   }
	return nbhcorrection;
}

/* returns all nodes that are descendants of the given node */
/* note: ancMatrix is 1 at [i,j] if i is an ancestor of j in the tree */
std::vector<int> getDescendants(bool** ancMatrix, int node, int n){
  std::vector<int> descendants;
  for(int i=0; i<n; i++){
  	if(ancMatrix[node][i]==true){
			descendants.push_back(i);
		}
	}
	return descendants;
}

/* returns all nodes that are not descendants of the given node */
/* i.e. ancestors and nodes in a different branch of the tree   */
/* note: ancMatrix is 0 at [i,j] if i is not an ancestor of j in the tree */
std::vector<int> getNonDescendants(bool**& ancMatrix, int node, int n){
	std::vector<int> ancestors;
	for(int i=0; i<n; i++){
		if(ancMatrix[node][i]==false){
			ancestors.push_back(i);
		}
	}
	return ancestors;
}


/* re-orders the nodes so that the descendant is first in case the nodes are in the same lineage */
int* reorderToStartWithDescendant(int* nodestoswap, bool** currTreeAncMatrix){
    if(currTreeAncMatrix[nodestoswap[0]][nodestoswap[1]]==true){  // make sure we move the descendent first
        int temp = nodestoswap[0];
        nodestoswap[0] = nodestoswap[1];
        nodestoswap[1] = temp;
    }
    return nodestoswap;
}


/* creates the new parent vector after swapping two nodes */
int* getNewParentVec_Swap(int* currTreeParentVec, int first, int second, int n, int* propTreeParVec){
    for(int i=0; i<n; i++){           // set vector of proposed parents
      if(propTreeParVec[i] == first && i!=second){
      	propTreeParVec[i] = second;              // update entries with swapped parents
      }
      else if(propTreeParVec[i] == second && i!=first){
      	propTreeParVec[i] = first;
      }
    }

    int temp = propTreeParVec[first];
    propTreeParVec[first] = propTreeParVec[second];  // update parents of swapped nodes
    propTreeParVec[second] = temp;
    if(propTreeParVec[first]==first){propTreeParVec[first]=second;}    // this is needed to ensure that tree is connected, the above fails
    if(propTreeParVec[second]==second){propTreeParVec[second]=first;}  // if first is parent of second, or vice versa
    return propTreeParVec;
}


/* creates the new ancestor matrix after swapping two nodes, old matrix is kept */
bool** getNewAncMatrix_Swap(bool** currTreeAncMatrix, int first, int second, int n, bool** propTreeAncMatrix){

    for(int i=0; i<n; i++){                                       // swap columns
    	bool temp = propTreeAncMatrix[i][first];
      	propTreeAncMatrix[i][first] = propTreeAncMatrix[i][second];
      	propTreeAncMatrix[i][second] = temp;
      }
      for(int i=0; i<n; i++){                                // swap rows
      	bool temp = propTreeAncMatrix[first][i];
      	propTreeAncMatrix[first][i] = propTreeAncMatrix[second][i];
      	propTreeAncMatrix[second][i] = temp;
      }
     return propTreeAncMatrix;
 }


/* creates the new parent vector after pruning and reattaching subtree */
int* getNewParentVec(int* currTreeParentVec, int nodeToMove, int newParent, int n, int *propTreeParVec){
    propTreeParVec[nodeToMove] = newParent;                       // only the parent of the moved node changes
    return propTreeParVec;
}


/* creates the new ancestor matrix after pruning and reattaching subtree, old matrix is kept */
bool** getNewAncMatrix(bool** currTreeAncMatrix, int newParent, std::vector<int> descendants, std::vector<int> possibleParents, int n, bool** propTreeAncMatrix){

    if(newParent<n){    // replace the non-descendants of the node and its descendants by the ancestors of the new parent
		for(int i=0; i<possibleParents.size(); i++){
  			for(int j=0; j<descendants.size(); j++){
  				propTreeAncMatrix[possibleParents[i]][descendants[j]] = currTreeAncMatrix[possibleParents[i]][newParent];
  			}
  		}
    }
    else
    {     // if we attach to the root, then they have no further ancestors
        for(int i=0; i<possibleParents.size(); i++){
  			for(int j=0; j<descendants.size(); j++){
        		propTreeAncMatrix[possibleParents[i]][descendants[j]] = 0;
        	}
        }
    }
    return propTreeAncMatrix;
}


/* picks a parent randomly from the set of possible parents, this set includes the root (n+1) */
int choseParent(std::vector<int> &possibleParents, int root){
	possibleParents.push_back(root);                           // add root, as it is also possible attachement point
    int chosenParentPos = pickRandomNumber(possibleParents.size());  // choose where to append the subtree
    int newParent = possibleParents[chosenParentPos];
	possibleParents.pop_back();    // remove root from list of possible parents as it is treated as special case later on
	return newParent;
}


/* distance is number of wrong parents */
int getSimpleDistance(int* trueVector, int* predVector, int n){
	int dist = 0;
	for(int i=0; i<n; i++){
		if(trueVector[i]!=predVector[i]){
			dist++;
		}
	}
	return dist;
}




