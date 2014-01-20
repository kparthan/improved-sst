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
    //!i Checks for a chain break
    bool checkChainBreak(string &, vector<Atom> &);

    //! Gets the number of chain breaks
    int getNumberOfChainBreaks(string &, vector<Atom> &);

    //! Computes the transformation at an index
    vector<Point<double>> computeTransformation(int, int);

    //! Reads the profile from a file
    void read_profile(string &);

    //! Initializes code length matrices
    void initializeCodeLengthMatrices(int);

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

    //! Computes the optimal code length matrix
    void computeCodeLengthMatrix();
};

#endif

