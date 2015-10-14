/*
 * output.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 */

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <float.h>
#include "output.h"
#include "scoreTree.h"
#include "matrices.h"

using namespace std;



/* writes the given string to file */
void writeToFile(string content, string fileName){
	ofstream outfile;
	outfile.open (fileName.c_str());
	outfile << content;
	outfile.close();
}

/* creates the content for the GraphViz file from a parent vector, using numbers as node labels (from 1 to n+1) */
std::string getGraphVizFileContentNumbers(int* parents, int n){
	std::stringstream content;
	content << "digraph G {\n";
	content << "node [color=deeppink4, style=filled, fontcolor=white];\n";
	for(int i=0; i<n; i++){
		content << parents[i]+1  << " -> "  << i+1 << ";\n";      // plus 1 to start gene labeling at 1 (instead of 0)
	}
	content <<  "}\n";
	return content.str();
}


/* creates the content for the GraphViz file from a parent vector, using the gene names as node labels */
std::string getGraphVizFileContentNames(int* parents, int n, vector<string> geneNames, bool attachSamples, bool** ancMatrix, int m, double** logScores, int** dataMatrix){
	std::stringstream content;
	content << "digraph G {\n";
	content << "node [color=deeppink4, style=filled, fontcolor=white];\n";

	for(int i=0; i<n; i++){
		content << geneNames[parents[i]] << " -> "  << geneNames[i]  << ";\n";
	}

	if(attachSamples==true){

		content << "node [color=lightgrey, style=filled, fontcolor=black];\n";
		std::string attachment = getBestAttachmentString(ancMatrix, n, m, logScores, dataMatrix, geneNames);
		content << attachment;
	}
	content <<  "}\n";
	return content.str();
}

/* creates the attachment string for the samples, the optimal attachment points are recomputed from scratch based on error log Scores */
std::string getBestAttachmentString(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix, vector<string> geneNames){
	bool** matrix = attachmentPoints(ancMatrix, n, m, logScores, dataMatrix);
	std::stringstream a;
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			if(matrix[i][j]==true){
				a << geneNames[i] << " -> s" << j << ";\n";
			}
		}
	}
	return a.str();
}

/* This is a re-computation of the best attachment points of the samples to a tree for printing the tree with attachment points */
/*   gets an ancestor matrix and returns a bit matrix indicating the best attachment points of each sample based on the error log scores */
bool** attachmentPoints(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix){

    double treeScore = 0.0;
    bool ** attachment = init_boolMatrix(n, m, false);
  	for(int sample=0; sample<m; sample++){       // foreach sample
  		double bestAttachmentScore = 0.0;     // currently best score for attaching sample
  		for(int gene=0; gene<n; gene++){   // start with attaching node to root (no genes mutated)
  			bestAttachmentScore += logScores[dataMatrix[sample][gene]][0];
  		}
  		for(int parent=0; parent<n; parent++){      // try all attachment points (genes)
  		    double attachmentScore=0.0;
  		    for(int gene=0; gene<n; gene++){     // sum up scores for each gene, score for zero if gene is not an ancestor of parent, score for one else wise
  		    	attachmentScore += logScores[dataMatrix[sample][gene]][ancMatrix[gene][parent]];
  		    }
  		    if(attachmentScore > bestAttachmentScore){
  		        bestAttachmentScore = attachmentScore;
  		    }
  		}
  		for(int parent=0; parent<n; parent++){      // try all attachment points (genes)
  		 	double attachmentScore=0.0;
  		 	for(int gene=0; gene<n; gene++){     // sum up scores for each gene, score for zero if gene is not an ancestor of parent, score for one else wise
  		 		attachmentScore += logScores[dataMatrix[sample][gene]][ancMatrix[gene][parent]];
  		 	}
  		  	if(attachmentScore == bestAttachmentScore){
  		  		attachment[parent][sample] = true;
  		  	}
  		}
  		treeScore += bestAttachmentScore;
  	}
  	return attachment;
}


/* prints all trees in list of optimal trees to the console, first as parent vector, then as GraphViz file */
void printParentVectors(vector<bool**> optimalTrees, int n, int m, double** logScores, int** dataMatrix){
	for(int i=0; i<optimalTrees.size(); i++){
		int* parents = ancMatrixToParVector(optimalTrees[i], n);
		print_intArray(parents,n);
		//print_boolMatrix(attachmentPoints(optimalTrees[i], n, m, logScores, dataMatrix), n, m);
		printGraphVizFile(parents, n);
	}
}


/* prints the GraphViz file for a tree to the console */
void printGraphVizFile(int* parents, int n){
	cout << "digraph G {\n";
	cout << "node [color=deeppink4, style=filled, fontcolor=white];\n";
	for(int i=0; i<n; i++){
		cout << parents[i] << " -> " << i << "\n";
	}
	cout << "}\n";
}

void printSampleTrees(vector<int*> list, int n, string fileName){
	if(list.size()==0){ return;}
	std::stringstream a;
	for(int i=0; i<list.size(); i++){
		for(int j=0; j<n; j++){
			a << list[i][j];
			if(j<n-1){
				a  << " ";
			}
		}
		a << "\n";
	}
	writeToFile(a.str(), fileName);
	cout << "Trees written to: " << fileName;
}

/* prints the score of the tree predicted by the Kim&Simon approach for the given error log scores */
void printScoreKimSimonTree(int n, int m, double** logScores, int** dataMatrix, char scoreType){
	int parent[] = {2, 4, 17, 2, 9, 9, 2, 2, 4, 18, 2, 1, 2, 2, 9, 2, 2, 11};
	double KimSimonScore = scoreTree(n, m, logScores, dataMatrix, scoreType, parent, -DBL_MAX);
	cout.precision(20);
	cout << "KimSimonScore: " << KimSimonScore << "\n";
}


