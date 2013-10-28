#ifndef SUPPORT_H
#define SUPPORT_H

#include "Header.h"
#include "Protein.h"

struct Parameters
{
  string file;
  int force;
  string profiles_dir;
  double res;
  int read_profiles;
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

// Protein functions
string getPDBFilePath(string &);
string getSCOPFilePath(string &);
void buildAngularProfile(struct Parameters &);
bool checkIfSphericalProfileExists(string &);
ProteinStructure *parsePDBFile(string &);
array<double,3> computeVonMisesEstimates(array<double,3> &, double);

array<double,3> readProfiles(string &, double);
void updateLogFile(string &, double, int);
void updateEstimator(array<double,3> &, double *, Protein &);
void updateBins(vector<vector<int>> &, double, Protein &);
void outputBins(vector<vector<int>> &, double);
void vonMisesDistribution_2DPlot(array<double,3> &);

#endif

