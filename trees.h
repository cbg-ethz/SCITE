/*
 * trees.h
 *
 *  Created on: Oct 12, 2015
 *      Author: jahnka
 */

#ifndef TREES_H_
#define TREES_H_
using namespace std;

vector<vector<int> > getChildListFromParentVector(int* parents, int n);
string getNewickCode(vector<vector<int> > list, int root);
int* prueferCode2parentVector(int* code, int codeLength);
int* getBreadthFirstTraversal(int* parent, int n);
bool** parentVector2ancMatrix(int* parent, int n);
int* getRandParentVec(int n);
bool* getInitialQueue(int* code, int codeLength);
int* getLastOcc(int* code, int codeLength);
int getNextInQueue(bool* queue, int pos, int length);
void updateQueue(int node, bool* queue, int next);
int updateQueueCutter(int node, bool* queue, int next);
int* starTreeVec(int n);
bool** starTreeMatrix(int n);

#endif /* TREES_H_ */
