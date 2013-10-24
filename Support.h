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
array<double,3> computeSphericalValues(Point<double> &);
array<double,3> computeVonMisesEstimates(array<double,3> &, double);
void print(array<double,3> &);
void generateHeatMap(array<double,3> &);
array<double,3> convertToCartesian(double, double, double);

// Protein functions
string getPDBFilePath(string &);
string getSCOPFilePath(string &);
void buildAngularProfile(struct Parameters &);
bool checkIfSphericalProfileExists(string &);
ProteinStructure *parsePDBFile(string &);
void updateLogFile(string &, double, int);

array<double,3> readProfiles(string &, double);
void updateEstimator(array<double,3> &, double *, Protein &);
void updateBins(vector<vector<int>> &, double, Protein &);

#endif

