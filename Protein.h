#ifndef PROTEIN_H
#define PROTEIN_H

#include "Header.h"
#include "Mixture.h"
#include "OptimalFit.h"

class Protein
{
  private:
    //! Protein identifier
    string name;

    //! List of suitable chain identifiers
    vector<string> chains;

    //! Protein Structure
    ProteinStructure *structure;

    //! Translation point
    Point<double> initial_translation_vector;

    //! Stores the coordinates
    vector<vector<Point<double>>> coordinates;

    //! Distances between the successive residues
    vector<vector<double>> distances;

    //! List of spherical coordinates (r \neq 1,theta,phi)
    //! theta,phi measured in degrees
    vector<vector<array<double,3>>> spherical_coordinates;

    //! Gets the list of all spherical coordinates
    //! theta,phi measured in degrees
    vector<array<double,3>> all_spherical_coordinates;

    //! Times taken
    double cpu_time,wall_time;

    //! Optimal model matrix
    vector<vector<OptimalFit>> optimal_model;

    //! Optimal code length matrix
    vector<vector<double>> optimal_code_length;

  protected:
    //! Checks for a chain break
    bool checkChainBreak(string &, vector<Atom> &);

    //! Gets the number of chain breaks
    int getNumberOfChainBreaks(string &, vector<Atom> &);

    //! Translate the protein so that its first point is the origin
    void translateProteinToOrigin(vector<Point<double>> &);

    //! Computes the transformation at an index
    vector<Point<double>> computeTransformation(int, int);

    //! Reads the profile from a file
    void read_profile(string &);

    //! Initializes code length matrices
    void initializeCodeLengthMatrices(int);

    //! Prints the optimal code length matrix to a file
    void printCodeLengthMatrix(int);

  public: 
    //! Null constructor
    Protein();

    //! Constructor
    Protein(ProteinStructure *,string &);

    //! Computes the successive distances between atoms
    void computeSuccessiveDistances();

    //! Gets the list of spherical coordinates
    void computeSphericalTransformation();

    //! Loads the spherical system
    void load(string &);

    //! Loads the spherical system
    void load(path &);

    //! Saves the spherical system
    void save();

    //! Gets the list of all spherical coordinates
    vector<array<double,3>> getSphericalCoordinatesList();

    //! Get the CPU time
    double getCPUTime();

    //! Gets the number of usable chains
    int getNumberOfChains();

    //! Gets the size of spherical coordinates list 
    int getNumberOfSphericalCoordinates();

    //! Computes the mean direction
    array<double,3> computeMeanDirection();

    //! Computes the message length using sphere model
    double computeMessageLengthUsingSphereModel();

    //! Computes the message length using null model
    double computeMessageLengthUsingNullModel(Mixture &);

    //! Compression using ideal models
    void compressUsingIdealModels(Mixture &, int);

    //! Computes the optimal code length matrix
    void computeCodeLengthMatrix(vector<IdealModel> &, Mixture &, int, int);

    //! Computes the optimal segmentation using dynamic programming
    pair<double,vector<int>> computeOptimalSegmentation(int);

};

#endif

