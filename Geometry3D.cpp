/* This file contains some useful functions for 3D Geometry manipulation */ 
#include "Geometry3D.h"
#define SQ(x) (x*x)

extern vector<double> ORIGIN;

/* This function computes direction cosines of a vector given its 
 * direction ratios. */ 
void computeDirectionCosines(vector<double> &dratios, vector<double> &dcosines)
{
	double sqsum = 
		1/(sqrt( SQ(dratios[0]) + SQ(dratios[1]) + SQ(dratios[2])) ); 
	for( int i = 0 ; i < 3 ; i++ ) dcosines[i] = dratios[i] * sqsum ; 
}

/* This function computes the direction ratios of two vectors. */
void computeDirectionRatios(vector<double> &A, vector<double> &B,
                            vector<double> &dratios)
{
  for (int i=0; i<3; i++) {
    dratios[i] = A[i] - B[i];
  }
}


/* This function computes the cross product of two vectors A and B */
void computeCrossProduct(vector<double> &dratiosA, vector<double> &dratiosB, 
                         vector<double> &dratnormal)
{ 
	dratnormal[0] = (dratiosA[1]*dratiosB[2])-(dratiosA[2]*dratiosB[1]);
	dratnormal[1] = (dratiosA[2]*dratiosB[0])-(dratiosA[0]*dratiosB[2]);
	dratnormal[2] = (dratiosA[0]*dratiosB[1])-(dratiosA[1]*dratiosB[0]);
}


/* This functions computes the unit vector along the normal to a plane 
 * formed by two lines A and B */
void computeNormal(vector<double> &dratiosA, vector<double> &dratiosB, 
                   vector<double> &dcosnormal)
{
	vector<double> dratnormal(3,0) ; 
	computeCrossProduct(dratiosA,dratiosB,dratnormal) ; 
	computeDirectionCosines( dratnormal, dcosnormal ) ; 
} 


/* This function computes the dot product of two vectors A and B */
void computeDotProduct(vector<double> &dratiosA, vector<double> &dratiosB, 
                       double& dotproduct)
{
	dotproduct = 0.0 ; 
	for( int i = 0 ; i < 3 ; i++ ) dotproduct += dratiosA[i]*dratiosB[i];
}

double normAminusB(vector<double> &A, vector<double> &B)
{
	double dratios[3] ;
	dratios[0] = A[0] - B[0] ;
	dratios[1] = A[1] - B[1] ;
	dratios[2] = A[2] - B[2] ;

	return sqrt( SQ(dratios[0]) + SQ(dratios[1]) + SQ(dratios[2]) )  ;
}

double vectorNorm( vector<double> &A ) 
{
	return normAminusB( A, ORIGIN ) ;
}

double computeAngle(vector<double> &A, vector<double> &B)
{

	double dp  = 0.0 ;
       	computeDotProduct( A, B, dp) ;

	double normA = normAminusB( A, ORIGIN ) ;
	double normB = normAminusB( B, ORIGIN ) ;
	//cout << dp << " " << normA << " " << normB << endl ;
	double cost = (dp/ (normA*normB)) ;
	  if( cost > 1 ) 
	  {
		  cost = 1 ;
	  }
	  else if( cost < -1 ) 
	  {
		  cost = -1 ;
	  }
	  double theta = acos( cost ) ; 

	  //find sign 
	  /*
    double normalAB[3] ;
	  double sign = 0.0 ;
	  computeNormal( A, B, normalAB ) ; 

	  computeBoxProduct( normalAB, A, B, sign ) ; 
	  theta = ( sign > 0 ) ? theta : -theta; 
    */
	  //theta *= 180/M_PI ; 

	return ( theta )	;
}

/* This function computes the rotation matrix for rotation through an angle 
 * theta about a line whose direction cosines are given by dcosines */
void computeRotationMatrix(vector<double> &dcosines, double theta, 
                           vector<vector<double>> &rotation_matrix)
{
	double cost = cos(theta) ; 
	double sint = sin(theta) ; 
	double a1 = dcosines[0] * sint ; 
	double a2 = dcosines[1] * sint ; 
	double a3 = dcosines[2] * sint ; 
	double b1 = dcosines[1] * dcosines[2] * ( 1 - cost ) ; 
	double b2 = dcosines[2] * dcosines[0] * ( 1 - cost ) ; 
	double b3 = dcosines[0] * dcosines[1] * ( 1 - cost ) ; 
	for( int i = 0 ; i < 3 ; i++ ) {
		rotation_matrix[i][i] = cost 
			+ ( dcosines[i] * dcosines[i] *(1-cost) ) ; 
	}
	rotation_matrix[0][1] = b3 - a3 ; 
	rotation_matrix[1][0] = b3 + a3 ; 
	rotation_matrix[2][0] = b2 - a2 ; 
	rotation_matrix[0][2] = b2 + a2 ; 
	rotation_matrix[1][2] = b1 - a1 ; 
	rotation_matrix[2][1] = b1 + a1 ; 
}

/* This function rotates initial_vector to final_vector using the rotation 
 * matrix */
void rotateVector(vector<vector<double>> &rotation_matrix, 
                  vector<double> &initial_vector, vector<double> &final_vector) 
{
	for( int i = 0 ; i < 3 ; i++ ) { 
		final_vector[i] = 0.0 ; 
		for( int j = 0 ; j < 3 ; j++ ) {
			final_vector[i] += (initial_vector[j] 
					* rotation_matrix[i][j] ) ;
		}
	}	
}

/* This function multiplies two matrices. */
void multiplyVectors(vector<vector<double>> &A, vector<vector<double>> &B,
                     vector<vector<double>> &product) 
{
  for (int i=0; i<A.size(); i++) {
    rotateVector(B,A[i],product[i]);
  }
}

double computeEuclideanDistance(vector<double> &v1, vector<double> &v2)
{
  d = 0;
  for (int i=0; i<3; i++) {
    d += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  }
  return sqrt(d);
}

