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

/* Score contribution by a specific mutation when placed at the root, that means all samples should have it */
/* This is the same for all trees and can be precomputed */
double binTreeRootScore(int** obsMutProfiles, int mut, int m, double ** logScores){
	double score = 0.0;
	for(int sample=0; sample<m; sample++){
		score += logScores[obsMutProfiles[sample][mut]][1];
	}
	return score;
}

/* computes the best placement of a mutation, the highest one if multiple co-opt. placements exist*/
int getHighestOptPlacement(int** obsMutProfiles, int mut, int m, double ** logScores, bool** ancMatrix){

	int nodeCount = (2*m)-1;
	int bestPlacement = (2*m)-2;   // root
	double bestPlacementScore = binTreeRootScore(obsMutProfiles, mut, m, logScores);
	//cout << bestPlacementScore << " (root)\n";
	//print_boolMatrix(bool** array, int n, int m);
	for(int p=0; p<nodeCount-1; p++){                           // try all possible placements (nodes in the mutation tree)

		double score = 0.0;                   // score for placing mutation at a specific node
		for(int sample=0; sample<m; sample++){
			//cout << p << " " << sample << "\n";
			if(ancMatrix[p][sample] == 1){
				score += logScores[obsMutProfiles[sample][mut]][1]; // sample should have the mutation
			}
			else{
				score += logScores[obsMutProfiles[sample][mut]][0]; // sample should not have the mutation
			}
		}
		if(score > bestPlacementScore){
			bestPlacement = p;
			bestPlacementScore = score;
			//cout << bestPlacementScore << " (non-root)\n";
		}
		else if (score == bestPlacementScore && ancMatrix[p][bestPlacement] == true){
			bestPlacement = p;
		}
	}

	//if(bestPlacement == (2*m)-2){
	//	cout<< "best placed at root\n";
	//	getchar();
	//}
	return bestPlacement;
}

/* computes the best placement of a mutation, the highest one if multiple co-opt. placements exist*/
int* getHighestOptPlacementVector(int** obsMutProfiles, int n, int m, double ** logScores, bool** ancMatrix){
	int* bestPlacements = init_intArray(n, -1);
	for(int mut=0; mut<n; mut++){                                                               // for all mutations get
		bestPlacements[mut] = getHighestOptPlacement(obsMutProfiles, mut, m, logScores, ancMatrix);         // bestPlacementScore
	 }
	//print_intArray(bestPlacements, n);
	return bestPlacements;
}

vector<string> getBinTreeNodeLabels(int nodeCount, int* optPlacements, int n, vector<string> geneNames){
	vector<string> v;
	int count = 0;
	for(int i = 0; i < nodeCount; i++){
		v.push_back("");
	}

	for(int mut=0; mut<n; mut++){
		string toAppend;
		if(v.at(optPlacements[mut]) == ""){
			toAppend = geneNames.at(mut);
			count++;
		}
		else{
			toAppend = ", " + geneNames.at(mut);
			count++;
		}
		//cout << "        " << j << "\n";
		//cout << "                     "<< optPlacements[j] << "\n";
		v.at(optPlacements[mut]) += toAppend;
	}
	if(v.at(nodeCount-1) == ""){
		v.at(nodeCount-1) = "root";
	}
	for(int i = 0; i < nodeCount; i++){
		if(v.at(i).find(" ") != string::npos){
			v.at(i) = "\"" + v.at(i) + "\"";
		}
	}
	//cout << "added mutations " << count << "\n";
	return v;
}

/* returns the lca of a node that has a non-empty label, the root is assumed to always have a label */
int getLcaWithLabel(int node, int* parent, vector<string> label, int nodeCount){
	int root = nodeCount -1;
	int p = parent[node];;
	while(p != root && label[p]==""){
		p = parent[p];
	}
	return p;
}

std::string getGraphVizBinTree(int* parents, int nodeCount, int m, vector<string> label){
	std::stringstream content;
	content << "digraph G {\n";
	content << "node [color=deeppink4, style=filled, fontcolor=white, fontsize=20, fontname=Verdana];\n";
	for(int i=m; i<nodeCount-1; i++){
		if(label[i] != ""){
			int labelledLCA = getLcaWithLabel(i, parents, label, nodeCount);
			content << label[labelledLCA] << " -> " << label[i] << ";\n";
//		if(label[parents[i]] == ""){
//			content  << parents[i] << " -> ";
//		}
//		else{
//			content << label[parents[i]] << " -> ";
//		}
//		if(label[i] == ""){
//		  content  << i << ";\n";
//		}
//		else{
//			content << label[i] << ";\n";

		}
	}
	content << "node [color=lightgrey, style=filled, fontcolor=black];\n";
	for(int i=0; i<m; i++){
		int labelledLCA = getLcaWithLabel(i, parents, label, nodeCount);
		content << label[labelledLCA] << " -> " << "s" << i << ";\n";



//		if(label[parents[i]] == ""){
//			content << parents[i] << " -> ";
//		}
//		else{
//			content << label[parents[i]] << " -> ";
//		}
//
//		content << "s" << i << ";\n";


	}
	content <<  "}\n";
	return content.str();
}



string getMutTreeGraphViz(vector<string> label, int nodeCount, int m, int* parent){
	stringstream nodes;
	stringstream leafedges;
	stringstream edges;
	for(int i=0; i<m; i++){
		if(label.at(i) != ""){
			nodes << "s" << i << "[label=\"s" << i << "\"];\n";                 // si [label="si"];
			nodes        << i << "[label=\"" << label.at(i) << "\"];\n";                 //   i [label="i"];
			leafedges << "s" << i << " -> " << i << ";\n";
			edges <<        i << " -> " << getLcaWithLabel(i, parent, label, nodeCount) << ";\n";
		}
		else{
			nodes << i << "[label=\"s" << i << "\"];\n";
			leafedges << i << " -> " << getLcaWithLabel(i, parent, label, nodeCount) << ";\n";
		}
	}

	stringstream str;

	str << "digraph g{\n";
	str << nodes.str();
	str << "node [color=deeppink4, style=filled, fontcolor=white];	\n";
	str << edges.str();
	str << "node [color=lightgrey, style=filled, fontcolor=black];  \n";
	str << leafedges.str();
	str << "}\n";
	return str.str();
}

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
	for(int i=0; i<=n; i++){
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
    bool ** attachment = init_boolMatrix(n+1, m, false);
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
  		bool rootAttachment = true;
  		for(int parent=0; parent<n; parent++){
  			if(attachment[parent][sample] == true){
  				rootAttachment = false;
  				break;
  			}
  		}
  		if(rootAttachment == true){
  			attachment[n][sample] = true;
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


