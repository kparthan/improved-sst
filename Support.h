#ifndef SUPPORT_H
#define SUPPORT_H

#include "Header.h"
#include "Protein.h"

struct Parameters
{
  string file;
  int force;
  string profiles_dir;
  int read_profiles;
  int update_bins;
  double res;
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
void updateLogFile(string &, double, int);
void updateMeanDirection(array<double,3> &, double *, Protein &);
void updateBins(vector<vector<int>> &, double, Protein &);
void outputBins(vector<vector<int>> &, double);

#endif

