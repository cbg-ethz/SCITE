/*
 * matrices.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include "matrices.h"

using namespace std;

/*****    basic functions on 1D and 2D arrays   *****/

double getMaxEntry(double* array, int n){
	double maxEntry = -DBL_MAX;
	for(int i=0; i<n; i++){
		maxEntry = max(maxEntry, array[i]);
	}
	return maxEntry;
}

int** sumMatrices(int** first, int** second, int n, int m){
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			first[i][j] += second[i][j];
		}
	}
	return first;
}

int ** transposeMatrix(int** matrix, int n, int m){
	int ** transposed = allocate_intMatrix(m, n);
	for(int i=0; i<m; i++){
		for(int j=0; j<n;j++){
			transposed[i][j] = matrix[j][i];
		}
	}
	return transposed;
}

void addToMatrix(int** first, int** second, int n, int m){
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			first[i][j] += second[i][j];
		}
	}
}

int* ancMatrixToParVector(bool** anc, int n){
	int* parVec = new int[n];
	for(int i=0; i<n; i++){
		parVec[i] = n;
	}
	for(int i=0; i<n; i++){
		for(int k=0; k<n; k++){
			if(k!=i && anc[k][i]==true){  // k is true ancestor of i
				bool cand = true;
				for(int l=0; l<n; l++){
					if(l!=i && l!=k && anc[l][i]==true && anc[k][l]==true){   // k is ancestor of l, and l is ancestor of i
						cand = false;                                        // k is no longer candidate for being parent of i
						break;
					}
				}
				if(cand==true){           // no ancestor of i is descendant of k -> k is parent of i
						parVec[i] = k;
				}
			}

		}
	}
	return parVec;
}

/*   allocation  */
double** allocate_doubleMatrix(int n, int m){

    double** matrix = new double*[n];
    matrix[0] = new double[n*m];
	  for (int i=1; i<n; ++i)
    {
        matrix[i] = matrix[i-1] + m;
    }
    return matrix;
}

int** allocate_intMatrix(int n, int m){

    int** matrix = new int*[n];
    matrix[0] = new int[n*m];
    for (int i=1; i<n; ++i)
    {
        matrix[i] = matrix[i-1] + m;
    }
    return matrix;
}

bool** allocate_boolMatrix(int n, int m){

    bool** matrix = new bool*[n];
    matrix[0] = new bool[n*m];
    for (int i=1; i<n; ++i)
    {
        matrix[i] = matrix[i-1] + m;
    }
    return matrix;
}


/*  initialization  */

int* init_intArray(int n, int value){
	int* array = new int[n];
	for(int i=0; i<n; i++){
		array[i] = value;
	}
	return array;
}

double* init_doubleArray(int n, double value){
	double* array = new double[n];
	for(int i=0; i<n; i++){
		array[i] = value;
	}
	return array;
}

bool* init_boolArray(int n, bool value){
	bool* array = new bool[n];
	for(int i=0; i<n; i++){
		array[i] = value;
	}
	return array;
}

double** init_doubleMatrix(int n, int m, double value){

	  double** matrix = allocate_doubleMatrix(n, m);     // allocate

    for (int i=0; i<n; ++i)             // initialize
    {
         for (int j=0; j<m; ++j)
      {
        	matrix[i][j] = value;
    	}
    }
    return matrix;
}

int** init_intMatrix(int n, int m, int value){

    int** matrix = allocate_intMatrix(n, m);  // allocate

    for (int i=0; i<n; ++i)            // initialize
    {
         for (int j=0; j<m; ++j)
      {
        	matrix[i][j] = value;
    	}
    }
    return matrix;
}

void reset_intMatrix(int** matrix, int n, int m, int value){

    for (int i=0; i<n; ++i)            // initialize
    {
         for (int j=0; j<m; ++j)
      {
        	matrix[i][j] = value;
    	}
    }
}


bool** init_boolMatrix(int n, int m, bool value){

    bool** matrix = allocate_boolMatrix(n, m);     // allocate

    for (int i=0; i<n; ++i)             // initialize
    {
         for (int j=0; j<m; ++j)
      {
        	matrix[i][j] = value;
    	}
    }
    return matrix;
}


/*  deallocation  */

void delete_3D_intMatrix(int*** matrix, int n){
	for(int i=0; i<n; i++){
			delete [] matrix[i][0];
			delete [] matrix[i];
		}
		delete [] matrix;
}


void free_boolMatrix(bool** matrix){
    delete [] matrix[0];
    delete [] matrix;
}

void free_intMatrix(int** matrix){
    delete [] matrix[0];
    delete [] matrix;
}

void free_doubleMatrix(double** matrix){
    delete [] matrix[0];
    delete [] matrix;
}


/*  deep copying  */

bool** deepCopy_boolMatrix(bool** matrix, int n, int m){
    bool** deepCopy = init_boolMatrix(n,m, false);
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
    	      deepCopy[i][j] = matrix[i][j];
	      }
    }
    return deepCopy;
}

int** deepCopy_intMatrix(int** matrix, int n, int m){
    int** deepCopy = init_intMatrix(n,m, -1);
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
    	      deepCopy[i][j] = matrix[i][j];
	      }
    }
    return deepCopy;
}

double** deepCopy_doubleMatrix(double** matrix, int n, int m){
    double** deepCopy = init_doubleMatrix(n,m, -1);
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
        {
    	      deepCopy[i][j] = matrix[i][j];
	      }
    }
    return deepCopy;
}

int* deepCopy_intArray(int* array, int n){
	  int* deepCopy = new int[n];
	  for (int i=0; i<n; ++i)
    {
        deepCopy[i] = array[i];
    }
    return deepCopy;
}

double* deepCopy_doubleArray(double* array, int n){
	double* deepCopy = new double[n];
	for (int i=0; i<n; ++i)
    {
        deepCopy[i] = array[i];
    }
    return deepCopy;
}

bool identical_boolMatrices(bool** first, bool** second, int n, int m){
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			//cout << "[" << i << "," << j << "] ";
			if(first[i][j] != second[i][j]){
				cout << "matrices differ!!!!!!!!!!!!!!!!\n";
				getchar();
				return false;
			}
		}
		//cout << "\n";
	}
	return true;
}

/*  printing  */

void print_boolMatrix(bool** array, int n, int m){
	  for(int i=0; i<n; i++){
		    for(int j=0; j<m; j++){
			      cout << array[i][j] << " ";
		    }
		    cout << "\n";
	  }
}

void print_doubleMatrix(double** matrix, int n, int m){
	  for(int i=0; i<n; i++){
  		  for(int j=0; j<m; j++){
  			    cout << matrix[i][j] << "\t";
  		  }
  		  cout << "\n";
  	}
}


void print_intArray(int* array, int n){
	for (int i=0; i<n; ++i)
    {
       cout << array[i] << " ";
    }
    cout << "\n";
}

void print_intMatrix(int** matrix, int n, int m, char del){
	  for(int i=0; i<n; i++){
  		  for(int j=0; j<m; j++){
  			    cout << matrix[i][j] <<  del;
  		  }
  		  cout << "\n";
  	}
}



