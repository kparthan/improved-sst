#ifndef PROTEIN_H
#define PROTEIN_H

#include "Support.h"

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

    //! List of spherical coordinates
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

    //! Computes the spherical coordinate values
    array<double,3> computeSphericalValues(vector<Point<double>> &);

  public: 
    //! Null constructor
    Protein();

    //! Constructor
    Protein(ProteinStructure *,string &);

    //! Gets the list of spherical coordinates
    void computeSphericalTransformation();

    //! Loads the spherical system
    void load(string &);

    //! Saves the spherical system
    void save();

    //! Gets the list of all spherical coordinates
    vector<array<double,3>> getSphericalCoordinatesList();

    //! Get the CPU time
    double getCPUTime();

    //! Gets the number of usable chains
    int getNumberOfChains();

};

#endif

