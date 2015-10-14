/*
 * output.h
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

void writeToFile(std::string content, std::string fileName);
std::string getGraphVizFileContentNumbers(int* parents, int n);
std::string getGraphVizFileContentNames(int* parents, int n, std::vector<std::string> geneNames, bool attachSamples, bool** ancMatrix, int m, double** logScores, int** dataMatrix);
std::string getBestAttachmentString(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix, std::vector<std::string> geneNames);
bool** attachmentPoints(bool ** ancMatrix, int n, int m, double** logScores, int** dataMatrix);
void printSampleTrees(std::vector<int*> list, int n, std::string fileName);
void printGraphVizFile(int* parents, int n);
void printScoreKimSimonTree(int n, int m, double** logScores, int** dataMatrix, char scoreType);

#endif /* OUTPUT_H_ */
