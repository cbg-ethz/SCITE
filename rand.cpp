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


int sampleRandomMove(std::vector<double> prob){ // picks randomly one of the moves based on the move probabilities

    double percent = rand() % 100;
    double sum = prob[0];
    for(int i=0; i<prob.size()-1; i++){
        if(percent <= sum*100){
          return i+1;
        }
        sum += prob[i+1];
    }
    return prob.size();
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


