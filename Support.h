#ifndef SUPPORT_H
#define SUPPORT_H

#include "Header.h"

struct Parameters
{
  string file;
};

// general functions
struct Parameters parseCommandLineInput (int, char **); 
void Usage (const char *, options_description &);
bool checkFile(string &);
void writeToFile(vector<Point<double>> &, const char *);
string extractName(string &);
string getPDBFilePath(string &);
string getSCOPFilePath(string &);

// Protein functions
void buildAngularProfile(struct Parameters &);
ProteinStructure *parsePDBFile(string &);

#endif

