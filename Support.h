#ifndef SUPPORT_H
#define SUPPORT_H

#include "Header.h"
#include "Protein.h"

struct Parameters
{
  string file;              // pdb file
  int force;                // flag to force build the angular profile
  string profiles_dir;      // path to the directory containing the profiles
  int read_profiles;        // flag to read through the existing profiles
  int update_bins;          // flag to update bins for MATLAB visualization
  double res;               // resolution of the bins
  int mixture_model;        // flag to run mixture modelling algorithm
  int infer_num_components; // flag to infer the number of components
  int num_components;       // the number of mixture components
  int update_weights_new;   // flag to update weights using modified rule
};

// general functions
struct Parameters parseCommandLineInput (int, char **); 
void Usage (const char *, options_description &);
bool checkFile(string &);
void writeToFile(vector<Point<double>> &, const char *);
string extractName(string &);
array<double,3> convertToSpherical(Point<double> &);
array<double,3> convertToCartesian(double, double, double);
void print(ostream &, array<double,3> &);
void vonMisesDistribution_2DPlot(array<double,3> &);
double ratioBesselFunction(double);
double ratioBesselFunction_firstDerivative(double);
double ratioBesselFunction_secondDerivative(double);
double ratioBesselFunction_thirdDerivative(double);

// Protein functions
string getPDBFilePath(string &);
string getSCOPFilePath(string &);
void buildAngularProfile(struct Parameters &);
bool checkIfSphericalProfileExists(string &);
ProteinStructure *parsePDBFile(string &);

void computeEstimators(struct Parameters &);
pair<array<double,3>,double> readProfiles(struct Parameters &);
vector<array<double,2>> gatherData(struct Parameters &);
void updateLogFile(string &, double, int);
void updateMeanDirection(array<double,3> &, double *, Protein &);
void updateBins(vector<vector<int>> &, double, Protein &);
void outputBins(vector<vector<int>> &, double);

#endif

