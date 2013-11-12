#ifndef PROTEIN_H
#define PROTEIN_H

#include "Header.h"

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

    //! List of spherical coordinates (r \neq 1,theta,phi)
    vector<vector<array<double,3>>> spherical_coordinates;

    //! Gets the list of all spherical coordinates
    vector<array<double,3>> all_spherical_coordinates;

    //! Times taken
    double cpu_time,wall_time;

  protected:
    //!i Checks for a chain break
    bool checkChainBreak(string &, vector<Atom> &);

    //! Gets the number of chain breaks
    int getNumberOfChainBreaks(string &, vector<Atom> &);

    //! Computes the transformation at an index
    vector<Point<double>> computeTransformation(int, int);

    //! Reads the profile from a file
    void read_profile(string &);

  public: 
    //! Null constructor
    Protein();

    //! Constructor
    Protein(ProteinStructure *,string &);

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

};

#endif

