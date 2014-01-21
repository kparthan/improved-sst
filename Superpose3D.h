#ifndef SUPERPOSE3D_H
#define SUPERPOSE3D_H

#include "Header.h"
#include "suffstat.h"

/* A superposition class to find the optimal least-square 
 * transformation of two corresponding vectors sets. */
class Superpose3DClass {
 private:
   vector<vector<double> > vecSetA ; // movingset 
   vector<vector<double> > vecSetB ; // fixedset

   suffStatClass stats; //computes sufficientStatistics
   double rotacenterA[3] ; // movingset center-of-mass
   double rotacenterB[3] ; // fixedset  center-of-mass

   /* These are used when updating superposition using sufficient statistics */
   suffStatClass stats1, stats2;
   double cm1A[3], cm1B[3]; // First part centers of mass (two sets)
   double cm2A[3], cm2B[3]; // second parts centers of mass (two sets)
   double mu1_xm, mu1_ym, mu1_zm;
   double mu1_xp, mu1_yp, mu1_zp;
   double mu2_xm, mu2_ym, mu2_zm;
   double mu2_xp, mu2_yp, mu2_zp;
   suffStatClass updated_stats;

   double rmsd ;
   double eulerRotMat[3][3] ;

   double quat[4][4] ; // A 4x4 quaternion matrices
   double eigenValues[4] ;
   double eigenVectors[4][4] ;

   double INFINITYVAL ;
   size_t nVecs ;
   bool success_flag ;

   void computeRotationalCenter();
   void assignRotationalCenter(vector<double>, vector<double>);
   
   void updateCenterOfMasses();

   void computeQuaternionMatrix() ;
   void computeSufficientStatistics();
   void computeQuaternionMatrixUsingSufficientStats();
   void updateQuaternionMatrixUsingSufficientStats();
   void updateSufficientStatistics();

   void eigsrt() ;
   bool diagonalize() ;
   void computeRotationMatrix() ;
   bool computeLSqFit() ;
   bool computeLSqFit(vector<double>&rcA,vector<double>&rcB) ;
   bool updateLSqFit();

   //debug
   void printRotationalCenters() ;
   void printEigens() ;
   void printQuaternionMatrix() ;
   void printRotationMatrix() ;

 public:
   Superpose3DClass(vector<vector<double> >&, vector<vector<double> >&) ;
   Superpose3DClass(vector<vector<double> >&, vector<vector<double> >&, size_t) ;
   Superpose3DClass(vector<vector<double> >&, vector<vector<double> >&, vector<double>&, vector<double>&) ;
   Superpose3DClass(vector<vector<double> >&, vector<vector<double> >&, size_t, vector<double>&, vector<double>&) ;
   Superpose3DClass(suffStatClass&, suffStatClass&);
   Superpose3DClass(vector<double>&, vector<double>&);
   Superpose3DClass(suffStatClass&, vector<double>&, vector<double>&);
   double getRMSD() ;	
   void transformVector( vector<double> & ) ;
   void transformVectors( vector<vector<double> > & ) ;
   void copyRotationMatrixInto(double [][3]);
   void copyRotationalCentersInto(double [3], double[3]);
   suffStatClass getSufficientStatistics();
};
#endif
