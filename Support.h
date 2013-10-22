#ifndef SUPPORT_H
#define SUPPORT_H

#include "Header.h"
#include "Segmentation.h"
#include "Polygon.h"
#include "DistanceHistogram.h"
#include "Angles.h"
#include "Lengths.h"
#include "KnotInvariants.h"

struct Parameters
{
  string file;                      // path to the structure file
  int structure;                    // Protein or General or Test
  int portion_to_fit;               // Fit entire protein or a segment
  int print;                        // verbose or not 
  vector<string> end_points;        // end points of a segment
  vector<int> controls;             // # of intermediate control points
  string control_string;            // controls string
  int constrain_sigma;              // flag to constrain sigma
  double min_sigma;                 // lower limit of sigma
  double max_sigma;                 // upper limit of sigma
  int constrain_segment_length;     // flag to constrain segment length
  int max_segment_length;           // maximum allowed segment length
  int encode_deviations;            // mechanism to encode deviations
  int segmentation_only;            // flag to only do the segmentation
  int force_segmentation;           // flag to redo segmentation
  int force_profile;                // flag to rebuild profile 
  int segmentation;                 // type of segmentation
  int profile;                      // distance histogram or knot invariants
  int comparison;                   // comparison flag set/unset 
  vector<string> comparison_files;  // path to the structure files used
                                    // in comparison
  int database_comparison;          // flag for database comparison
  string database;                  // database file
  int record;                       // flag to record all experimental results
  // parameters used for profiling using dihedral angles
  int align_type;                   // alignment type
  int scoring_function;             // scoring function
  double gap_open_penalty;          // gap open penalty
  double gap_extension_penalty;     // gap extension penalty
  double gap_penalty;               // gap penalty used in basic alignment 
                                    // method
  double max_angle_diff;            // maximum allowed angle separation
                                    // when aligning two dihedral angles
  // parameters used for profiling using distance histogram 
  double scale;                     // scale the length/coordinates to obtain
                                    // # of samples
  int num_samples_on_curve;         // # of samples for histogram comparison
  int sampling_method;              // uniform/random generation of samples
  double increment_r;               // increment in r used in histogram method
  // parameters used for profiling using knot invariants
  string method;                    // general/specific
  int standardize;                  // flag for standardizing the invariants
  int standardize_status;           // whether the standardization exists or not
  string standardize_names;         // structure file to standardize premeasures
  string standardize_parameters;    // file to which the parameters are saved
  int construct_polygon;            // polygon construction heuristic 
  int num_sides;                    // # of sides in the approximating polygon
  int max_order;                    // maximum order of knot invariants
};

// segmentation functions
bool checkIfSegmentationExists(string &, string &);
Segmentation buildSegmentationProfile(struct Parameters &);
string getPDBFilePath(string &);
string getSCOPFilePath(string &);
ProteinStructure *parsePDBFile(string &);
vector<Point<double>> parseFile(string &);
Segmentation testFit(struct Parameters &);
vector<Point<double>> getProteinCoordinates(struct Parameters &);
Segmentation proteinFit(struct Parameters &);
Segmentation generalFit(struct Parameters &);
void updateRuntime(string, Segmentation &);

// angles functions
bool checkIfAnglesExist(string &, string &);
Polygon<double> getRepresentativePolygon(struct Parameters &, Segmentation &);
Angles buildAnglesProfile(struct Parameters &, Segmentation &);
double computeDihedralAngle(Line<double> &, Line<double> &);
void updateRuntime(string, Angles &, double);
void updateResults(struct Parameters &, vector<vector<double>> &);

// lengths functions
bool checkIfLengthsExist(string &, string &);
Lengths buildLengthsProfile(struct Parameters &, Segmentation &);
double computeMidPointsDistance(Line<double> &, Line<double> &);
void updateRuntime(string, Lengths &, double);

// histograms functions
bool checkIfHistogramExists(string &);
DistanceHistogram buildHistogramProfile(struct Parameters &, Segmentation &);
vector<double> getRValuesList(double, double);
void plotMultipleHistograms(vector<DistanceHistogram> &, vector<double> &, 
                            vector<string> &);
void printHistogramResults(vector<DistanceHistogram> &, vector<double> &, 
                           vector<string> &);

// knot invariants functions
bool checkIfKnotInvariantsExist(string &, string &);
KnotInvariants buildKnotInvariantsProfile(struct Parameters &, Segmentation &);
void standardizePremeasures(struct Parameters &);
pair<vector<double>,vector<double>> loadStandardizedParameters(string &);
void
updateExpectations(int &, vector<double> &, vector<double> &, vector<double> &);
template <typename RealType>
double exteriorAngle(Vector<double> &, Vector<double> &, Vector<double> &);
double sumExteriorAngles(Line<double> &, Line<double> &);
double computeEuclideanDistance(Vector<double> &, Vector<double> &);
void updateRuntime(string, int, double); 
void updateResults(vector<double> &, vector<double> &);

// sst functions
Angles buildSSTProfile(struct Parameters &);
vector<vector<string>> parse_segmentation(ProteinStructure *, string &);
Angles 
construct_angular_profiles(ProteinStructure *, vector<vector<string>> &, string &);
Lengths
construct_lengths_profiles(ProteinStructure *, vector<vector<string>> &, string &);
vector<Line<double>>
get_list_of_lines(ProteinStructure *, vector<vector<string>> &, string &);

// dssp functions
Angles buildDSSPProfile_Angles(struct Parameters &);
Lengths buildDSSPProfile_Lengths(struct Parameters &);
vector<pair<string,string>> split_segments(vector<vector<string>> &);
KnotInvariants build_DSSP_KnotInvariants(struct Parameters &);

// general functions
struct Parameters parseCommandLineInput (int, char **); 
void Usage (const char *, options_description &);
void build(struct Parameters &);
string getControlString(vector<int> &);
bool checkFile(string &);
string extractName(string &);
void writeToFile(vector<array<double,3>> &, const char*);
void compareStructuresList(struct Parameters &);
void compareDatabaseStructures(struct Parameters &);
vector<string> parseDatabase(string &);
void errorLog(vector<string> &);

template <typename RealType>
vector<RealType> sort(vector<RealType> &);
template <typename RealType>
vector<int> sortedListIndex(vector<RealType> &);
template <typename RealType>
void quicksort(vector<RealType> &, vector<int> &, int, int);
template <typename RealType>
int partition(vector<RealType> &, vector<int> &, int, int);

double estimateMean(vector<double> &);
double standardDeviation(vector<double> &, double);
double estimateVariance (vector<double> &, double);
double estimateVariance (vector<double> &);
int sign(double);
double cubeRoot(double);
template <typename RealType>
RealType absoluteMaximum(vector<RealType> &);
double minimum(vector<double> &);
double minimum(vector<vector<double>> &);
double bernstein(int, int, double);
double getMaximumDistance(vector<array<double,3>> &);

double getSelfAlignmentScore(Angles &, struct Parameters &);
double getSelfAlignmentScore(Angles &, Lengths &, struct Parameters &);
pair<double,double> computeNormalizedAlignmentScore(double, double, double);

#endif

