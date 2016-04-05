/*
 * rand.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#ifndef RAND_H
#define RAND_H

#include <vector>
//#include <string>

void initRand();
bool changeBeta(double prob);
int sampleRandomMove(std::vector<double> prob);
int* sampleTwoElementsWithoutReplacement(int n);
int pickRandomNumber(int n);
double sample_0_1();
int* getRandTreeCode(int n);
bool samplingByProb(double prob);
int* getRandomBinaryTree(int m);
#endif
