/* This file contains some useful functions for 3D Geometry manipulation */ 
#ifndef GEOMETRY3D_H
#define GEOMETRY3D_H

#include "Header.h"

/* This function computes direction cosines of a vector given its 
 * direction ratios. */ 
void computeDirectionCosines(vector<double> &, vector<double> &); 

void computeDirectionRatios(vector<double> &, vector<double> &, vector<double> &);

/* This functions computes the unit vector along the normal to a 
 * plane formed by two lines A and B */
void computeNormal(vector<double> &, vector<double> &, vector<double> &);


/* This function computes the cross product of two vectors A and B */
void computeCrossProduct(vector<double> &, vector<double> &, vector<double> &);


/* This function computes the dot product of two vectors A and B */
void computeDotProduct(vector<double> &, vector<double> &, double &); 

/*compute angle between 2 vectors A and B*/
double computeAngle(vector<double> &, vector<double> &);

double normAminusB(vector<double> &, vector<double> &);
double vectorNorm(vector<double> &);

/* This function computes the rotation matrix for rotation through an 
 * angle theta about a line whose direction cosines are given by dcosines 
*/
void computeRotationMatrix(vector<double> &, double, vector<vector<double>> &);

/* This function rotates initial_vector to final_vector using the 
 * rotation matrix */
void rotateVector(vector<vector<double>> &, vector<double> &, vector<double> &); 

void multiplyVectors(vector<vector<double>> &, vector<vector<double>> &,
                     vector<vector<double>> &); 

double computeEuclideanDistance(vector<double> &, vector<double> &);

#endif

