#ifndef SUPPORT_H
#define SUPPORT_H

#include "Header.h"
#include "Protein.h"
//#include "Component.h"
//#include "IdealModel.h"

struct Parameters
{
  string file;                  // pdb file
  int force;                    // flag to force build the angular profile
  string profiles_dir;          // path to the directory containing the profiles
  int constrain_kappa;          // flag to constrain kappa
  int read_profiles;            // flag to read through the existing profiles
  int heat_map;                 // flag to update bins for MATLAB visualization
  double res;                   // resolution of the bins
  int mixture_model;            // flag to run mixture modelling algorithm
  int infer_num_components;     // flag to infer the number of components
  int fit_num_components;       // the number of mixture components
  int update_weights_new;       // flag to update weights using modified rule
  // parameters to visualize mixture components
  int load_mixture;             // flag to load mixture details
  string mixture_file;          // path to the mixture file
  int sample_generation;        // using component weights/random sample size
  int num_samples;              // sample size to be generated
  // parameters to simulate the mixture modelling
  int simulation;               // flag to simulate
  int simulate_num_components;  // # of components to simulate
  // parameters to run sst
  int sst;                      // flag to run sst
  string structure;             // protein structure file
  int orientation;              // orientation to be used in the adaptive
                                // encoding scheme
};

struct Estimates
{
  vector<double> unit_mean;
  double kappa;
}

// general functions
void getHomeAndCurrentDirectory();
struct Parameters parseCommandLineInput (int, char **); 
void Usage (const char *, options_description &);
bool checkFile(string &);
void writeToFile(vector<vector<double>> &, const char *);
string extractName(string &);
void initializeMatrix(vector<vector<double>> &, int, int);
void cartesian2spherical(vector<double> &, vector<double> &);
void spherical2cartesian(vector<double> &, vector<double> &);
void point2vector(Point<double> &, vector<double> &);
void scaleToAOM(double *);
template <typename RealType> RealType minimum(RealType, RealType);
void print(ostream &, vector<double> &);
//void vonMisesDistribution_2DPlot(array<double,3> &);
double ratioBesselFunction(double);
double ratioBesselFunction_firstDerivative(double);
double ratioBesselFunction_secondDerivative(double);
double ratioBesselFunction_thirdDerivative(double);
double computeConstantTerm(int);
double getLatticeConstant(int);
double angleInRadians(double);

// Protein functions
string getPDBFilePath(string &);
string getSCOPFilePath(string &);
void buildAngularProfile(struct Parameters &);
bool checkIfSphericalProfileExists(string &);
ProteinStructure *parsePDBFile(string &);
void convertToCanonicalForm(vector<vector<double>> &, vector<vector<double>> &,
                            vector<vector<double>> &);
//Matrix<double> alignWithZAxis(vector<double> &, vector<double> &);
//array<double,3> applyIdealModelTransformation(Matrix<double> &, vector<double> &, vector<double> &);
//
//void computeEstimators(struct Parameters &);
//void modelOneComponent(struct Parameters &, pair<array<double,3>,double> &);
//void modelMixture(struct Parameters &, vector<array<double,3>> &);
//pair<array<double,3>,double> readProfiles(struct Parameters &);
//vector<array<double,3>> gatherData(struct Parameters &);
void updateLogFile(string &, double, int);
//void updateMeanDirection(array<double,3> &, double *, Protein &);
//void updateBins(vector<vector<int>> &, double, Protein &);
//void outputBins(vector<vector<int>> &, double);
//void visualizeMixtureComponents(struct Parameters &);
//void simulateMixtureModel(struct Parameters &);
//vector<double> generateRandomWeights(int, double);
//vector<Component> generateRandomComponents(int, int);
//void plotMessageLengthAgainstComponents(vector<int> &, vector<double> &, int);
//
//// sst functions
//void assignSecondaryStructure(string, string, int);
//vector<IdealModel> loadIdealModels();
//
#endif

