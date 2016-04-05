/*
 * rand_C.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

//#include <string>
//#include <random>
#include "rand.h"
#include <iostream>
#include <random>
#include <stdlib.h>
#include <time.h>
#include "matrices.h"

using namespace std;


/*****    functions for sampling random numbers inside C++  *****/
void initRand(){
	time_t t;
	time(&t);
	srand((unsigned int)t);              // initialize random number generator
	//srand(1);
}


/* This function gets a number of nodes n, and creates a random pruefer code for a rooted tree with n+1 nodes (root is always node n+1) */
int* getRandTreeCode(int n){                // as usual n is the number of mutations

	int nodes = n+1;                        // #nodes = n mutations plus root (wildtype)
	int codeLength = nodes-2;
	int* code = new int[codeLength];
	for(int i=0; i<codeLength; i++){
		code[i] = rand() % nodes;
	}
	return code;
}

bool changeBeta(double prob){
	 double percent = (rand() % 100)+1;    // between 1 and 100
	 if(percent <= prob*100){
		 return true;
	 }
	 return false;
}

int sampleRandomMove(std::vector<double> prob){ // picks randomly one of the tree moves based on the move probabilities

    double percent = (rand() % 100)+1;    // between 1 and 100
    double probSum = prob[1];
    for(int i=1; i<prob.size()-1; i++){    // start at index 1; the probability at prob[0] is for changing the error rate (which is treated separately)
        if(percent <= probSum*100){
          return i;
        }
        probSum += prob[i+1];
    }
    return prob.size()-1;
}


bool samplingByProb(double prob){
	double percent = rand() % 100;
	if(percent <= prob*100){
		return true;
	}
	return false;
}


int* sampleTwoElementsWithoutReplacement(int n){

    int* result = new int[2];
	  result[0] = rand() % n;
	  result[1] = result[0];
    while(result[0]==result[1]){
      result[1] = rand() % n;
    }
	  return result;
}

int pickRandomNumber(int n){

    return (rand() % n);
}

double sample_0_1(){

  //return (((double) rand()+0.5) / ((RAND_MAX+1)));
  return ((double) rand() / RAND_MAX);
}

int getElemFromQueue(int index, std::vector<int> queue){
	int elem = queue.at(index);
	if (index != queue.size() - 1)
	{
		queue[index] = std::move(queue.back());
	}

	//cout << queue.size() << " elements in queue in subroutine\n";
	return elem;
}

// This creates the parent vector of a random binary tree. Entries 0...m-1 are for the leafs.
// Entries m....2m-3 are for the inner nodes except the root, the root has index 2m-2 which has no parent
// and therefore has no entry in the parent vector
int* getRandomBinaryTree(int m){
	int parentCount = (2*m)-2;     // the m leafs have a parent and so have m-2 of the m-1 inner nodes
	int* leafsAndInnerNodesParents = init_intArray(parentCount, -1);

	std::vector<int> queue;
	for(int i=0; i<m; i++){queue.push_back(i);}   // add the m leafs to the queue
	//cout << queue.size() << " elements in queue\n";
	int innerNodeCount = m;
	while(queue.size()>1){
		int pos = pickRandomNumber(queue.size());
		int child1 = queue.at(pos);
		if (pos != queue.size() - 1){queue[pos] = std::move(queue.back());}
		queue.pop_back();
		//cout << queue.size() << " elements in queue\n";

		pos = pickRandomNumber(queue.size());
		int child2 = queue.at(pos);
		if (pos != queue.size() - 1){queue[pos] = std::move(queue.back());}
		queue.pop_back();
		//cout << queue.size() << " elements in queue\n";

		leafsAndInnerNodesParents[child1] = innerNodeCount;
		leafsAndInnerNodesParents[child2] = innerNodeCount;
		queue.push_back(innerNodeCount);
		innerNodeCount++;
	}
	return leafsAndInnerNodesParents;
}

