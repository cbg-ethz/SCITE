/*
 * scoreTree.cpp
 *
 *  Created on: Aug 15, 2015
 *      Author: jahnka
 */
#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <float.h>
#include <math.h>
#include <cmath>
#include <queue>
#include "matrices.h"
#include "treelist.h"
#include "trees.h"
#include "scoreTree.h"

using namespace std;

double epsilon = 0.000000000001;  // how much worse than the current best a score can be to still use the accurate score computation


/****     Tree scoring     ****/


/* Computes the score of a new candidate tree. First a fast approximate score is computed, then if the new score is better,   */
/* or slightly worse than the best score so far, a more accurate but more costly score computation is done.                    */
double scoreTree(int n, int m, double** logScores, int** dataMatrix, char type, int* parentVector, double bestTreeLogScore){

	double approx = scoreTreeFast(n, m, logScores, dataMatrix, type, parentVector);   // approximate score

	if(approx > bestTreeLogScore-epsilon){                                                  // approximate score is close to or better
		return scoreTreeAccurate(n, m, logScores, dataMatrix, type, parentVector);   // than the current best score, use accurate
	}                                                                                // score computation

	return approx;                                                              // otherwise the approximate score is sufficient
}



/****     Fast (approximate) tree scoring     ****/

/* computes an approximate score for a tree. This is fast, but rounding errors can occur  */
double scoreTreeFast(int n, int m, double** logScores, int** dataMatrix, char type, int* parentVector){

	double result = -DBL_MAX;
	int* bft = getBreadthFirstTraversal(parentVector, n);   // get breadth first traversal for simple scoring
	                                                        // by updating the parent score
	if(type=='m'){
		result = maxScoreTreeFast(n, m, logScores, dataMatrix, parentVector, bft);  // score by best attachment point per sample
	}
	else if(type=='s'){
		result = sumScoreTreeFast(n, m, logScores, dataMatrix, parentVector, bft);  // score by summing over all attachment points
	}

	delete [] bft;
	return result;
}

/* computes an approximate scoring for a tree using the max attachment score per sample */
double maxScoreTreeFast(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft){

    double treeScore = 0.0;

  	for(int sample=0; sample<m; sample++){                                                      // for all samples get
  		double* scores = getAttachmentScoresFast(parent,n, logScores, dataMatrix[sample], bft);  // all attachment scores
  		treeScore +=  getMaxEntry(scores, n+1);
  		delete [] scores;
  	}

  	return treeScore;    // sum over the best attachment scores of all samples is tree score
}


/* computes an approximate scoring for a tree summing the score over all attachment points per sample */
double sumScoreTreeFast(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft){

	double sumTreeScore = 0.0;

	for(int sample=0; sample<m; sample++){
		double* scores = getAttachmentScoresFast(parent, n, logScores, dataMatrix[sample], bft); // attachments scores of sample to each node
		double bestMaxTreeScore = getMaxEntry(scores, n+1);                                     // best of the scores (used to compute with score differences rather than scores)

		double sumScore = 0.0;
		for(int i=0; i<=n; i++){                                                 // sum over all attachment scores, exp is necessary as scores are actually log scores
			sumScore += exp(scores[bft[i]]-bestMaxTreeScore);                   // subtraction of best score to calculate with score differences (smaller values)
		}
		delete [] scores;
		sumTreeScore += log(sumScore)+bestMaxTreeScore;                     // transform back to log scores and change from score differences to actual scores
	}
	return sumTreeScore;
}


/* computes the attachment scores of a sample to all nodes in the tree (except root) */
double* getAttachmentScoresFast(int*parent, int n, double** logScores, int* dataVector, int*bft){

	double* attachmentScore = init_doubleArray(n+1, -DBL_MAX);
	attachmentScore[n] = rootAttachementScore(n, logScores, dataVector);
	for(int i=1; i<=n; i++){                                                              // try all attachment points (nodes in the mutation tree)
		int node = bft[i];
		attachmentScore[node] = attachmentScore[parent[node]];
		attachmentScore[node] -= logScores[dataVector[node]][0];
		attachmentScore[node] += logScores[dataVector[node]][1];
	}
	return attachmentScore;
}

/* computes the log score for attaching a sample to the root node (this score is equal for all trees) */
double rootAttachementScore(int n, double** logScores, int* dataVector){
	double score = 0.0;
	for(int gene=0; gene<n; gene++){                          // sum over log scores for all other nodes in tree
		score += logScores[dataVector[gene]][0];      // none of them is ancestor of the sample as it is attached to root
	}
	return score;
}





/****  accurate score computation (minimizes rounding errors)  ****/

double scoreTreeAccurate(int n, int m, double** logScores, int** dataMatrix, char type, int* parentVector){

	double result = -DBL_MAX;
	int* bft = getBreadthFirstTraversal(parentVector, n);
	if(type=='m'){
		result = maxScoreTreeAccurate(n, m, logScores, dataMatrix, parentVector, bft);
	}
	else if(type=='s'){
		result = sumScoreTreeAccurate(n, m, logScores, dataMatrix, parentVector, bft);
	}

	delete [] bft;
	return result;
}

double maxScoreTreeAccurate(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft){

    int** treeScoreMatrix = init_intMatrix(4, 2, 0);
  	for(int sample=0; sample<m; sample++){
  		int** bestAttachmentMatrix =  getBestAttachmentScoreAccurate(init_intMatrix(4, 2, 0), parent, n, logScores, dataMatrix[sample], bft);
  		treeScoreMatrix = sumMatrices(treeScoreMatrix, bestAttachmentMatrix, 4, 2);
  		free_intMatrix(bestAttachmentMatrix);
  	}
  	double treeScore = getTrueScore(treeScoreMatrix, logScores);
  	free_intMatrix(treeScoreMatrix);
  	return treeScore;
}

/* computes the log score for the complete tree using the sumScore scheme, where likelihoods of all attachment points of a sample are added */
double sumScoreTreeAccurate(int n, int m, double** logScores, int** dataMatrix, int* parent, int* bft){

	double sumTreeScore = 0.0;

	for(int sample=0; sample<m; sample++){
		sumTreeScore += getSumAttachmentScoreAccurate(parent, n, logScores, dataMatrix[sample], bft);
	}
	return sumTreeScore;
}

/* computes the best attachment score for a sample to a tree */
int** getBestAttachmentScoreAccurate(int** scoreMatrix, int* parent, int n, double** logScores, int* dataVector, int* bft){

	int*** attachmentScoreMatrix = getAttachmentMatrices(parent, n, dataVector, bft);   // matrix to keep attachment scores for each sample (not summing up
	                                                                                            // components to avoid rounding errors)
	double bestScore =  -DBL_MAX;
	int** bestScoreMatrix = NULL;

	for(int i=0; i<n+1; i++){                                                                   // now get true attachment scores and find best score among all attachment points
		double newScore = getTrueScore(attachmentScoreMatrix[i], logScores);
		if(bestScore <= newScore){
			bestScoreMatrix = attachmentScoreMatrix[i];
			bestScore = newScore;
		}
	}
	scoreMatrix = sumMatrices(scoreMatrix, bestScoreMatrix, 4, 2);
	delete_3D_intMatrix(attachmentScoreMatrix, n+1);
	return scoreMatrix;
}

/* computes the sum score for attaching a sample to all nodes */
double getSumAttachmentScoreAccurate(int* parent, int n, double** logScores, int* dataVector, int* bft){

	int*** attachmentScoreMatrix = getAttachmentMatrices(parent, n, dataVector, bft);   // matrix to keep attachment scores for each sample (not summing up
		                                                                                            // components to avoid rounding errors)
	double* attachmentScore = getTrueScores(attachmentScoreMatrix, n, logScores);                   // get the true attachment scores from the attachment matrices
	double bestScore = getMaxEntry(attachmentScore, n+1);                                            // identify best attachment score
	double sumScore = 0.0;
	for(int parent = 0; parent<n+1; parent++){                                                        // get score for attaching to the other nodes in the tree
		sumScore += exp(attachmentScore[parent]-bestScore);
	}
	delete_3D_intMatrix(attachmentScoreMatrix, n+1);
	delete [] attachmentScore;
	return log(sumScore)+bestScore;
}

/* computes the attachment scores of a sample to all nodes in a tree, score is a matrix counting the number of different match/mismatch score types */
int*** getAttachmentMatrices(int* parent, int n, int* dataVector, int* bft){
	int*** attachmentScoreMatrix = new int**[n+1];            // matrix to keep attachment scores for each sample (not summing up components to avoid rounding errors)

	// start with attaching node to root (no genes mutated)
	attachmentScoreMatrix[n] = init_intMatrix(4, 2, 0);
	for(int gene=0; gene<n; gene++){
		attachmentScoreMatrix[n][dataVector[gene]][0]++;
	}

	// now get scores for the other nodes due to bft traversal in an order such that attachment matrix of parent is already filled
	for(int i=1; i<n+1; i++){
		int node = bft[i];
		attachmentScoreMatrix[node] = deepCopy_intMatrix(attachmentScoreMatrix[parent[node]], 4, 2);
		attachmentScoreMatrix[node][dataVector[node]][0]--;
		attachmentScoreMatrix[node][dataVector[node]][1]++;
	}
	return attachmentScoreMatrix;
}

double* getTrueScores(int*** matrix, int n, double** logScores){
	double* scores = new double[n+1];
	for(int node=0; node<=n; node++){
		scores[node] = getTrueScore(matrix[node], logScores);
	}
	return scores;
}


/* computes the attachment score of a sample to a tree from a matrix */
/* representation of the score (counting the # 0->1, 1->0, 0->0, ...) */
double getTrueScore(int** matrix, double** logScores){
	double score = 0.0;
	for(int j=0; j<4; j++){
		for(int k=0; k<2; k++){
			double product = matrix[j][k] * logScores[j][k];
			//cout << "[" << j << "][" << k << "] = " << matrix[j][k] << " * " << logScores[j][k] << "[" << j <<"][" << k << "]\n";
			score = score + product;
		}
	}
	return score;
}

/***********************         Scoring Tables            ***************************/

/* computes a table of the log scores of observing one genotype, given that another genotype */
/* is the true one; for three observed types (+missing observation) and two true states */
double** getLogScores(double FD, double AD1, double AD2, double CC){

  double** logScores = init_doubleMatrix(4, 2, 0.0);
  logScores[0][0] = log(1.0-CC-FD);  // observed 0, true 0
	logScores[1][0] = log(FD);         // observed 1, true 0
	if(CC!=0.0){
		logScores[2][0] = log(CC);         // observed 2, true 0
	}
	else{                                //  to capture case where CC=0, because log(0) = -infinity
		logScores[2][0] = 0.0;           // CC=0 should only occur when there are no 2's (double mutations) in the matrix
	}
	logScores[3][0] = log(1.0);          // value N/A,  true 0
	logScores[0][1] = log(AD1);      // observed 0, true 1
	logScores[1][1] = log(1.0-(AD1+AD2));     // observed 1, true 1
	if(AD2 != 0.0){
		logScores[2][1] = log(AD2);     // observed 2, true 1
	}
	else{
		logScores[2][1] = 0;
	}
	logScores[3][1] = log(1.0);          // value N/A,  true 1
	return logScores;
}

/* updates the log scores after a new AD rate was accepted in the MCMC */
void updateLogScores(double** logScores, double newAD){

	double newAD1 = newAD;      // the default case: there are no homozygous mutation observed
	double newAD2 = 0.0;

	if(logScores[2][1] != 0){          // other case: homozygous mutations were called
		newAD1 = newAD/2;          // for simplicity we set both dropout rates to 1/2 of the new value
		newAD2 = newAD/2;          // learning the rates separately could also be done
	}

	logScores[0][1] = log(newAD1);          // observed 0, true 1
	logScores[1][1] = log(1.0-(newAD));     // observed 1, true 1
	if(newAD2 != 0.0){
		logScores[2][1] = log(newAD2);     // observed 2, true 1
	}
	else{
		logScores[2][1] = 0;
	}
	logScores[3][1] = log(1.0);          // value N/A,  true 1
}


double** getScores(double FD, double AD1, double AD2, double CC){

  double** scores = init_doubleMatrix(4, 2, 0.0);
  scores[0][0] = 1.0-CC-FD;  // observed 0, true 0
  scores[1][0] = FD;         // observed 1, true 0
  scores[2][0] = CC;         // observed 2, true 0
  scores[3][0] = 1.0;          // value N/A,  true 0
  scores[0][1] = AD1;          // observed 0, true 1
  scores[1][1] = 1.0-(AD1+AD2);     // observed 1, true 1
  scores[2][1] = AD2;               // observed 2, true 1
  scores[3][1] = 1.0;               // value N/A,  true 1
  return scores;
}

void printLogScores(double** logScores){
	cout.precision(70);
	for(int i=0; i<4; i++){
		for(int j=0; j<2; j++){
			cout << logScores[i][j] << "\t";
		}
		cout << "\n";
	}
}
