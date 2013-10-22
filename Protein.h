#ifndef PROTEIN_H
#define PROTEIN_H

#include "Support.h"

class Protein
{
  private:
    //! Protein Structure
    ProteinStructure *protein_structure;

    //! Stores the coordinates
    vector<Point<double>> coordinates;

    //! List of spherical coordinates
    vector<array<double,3>> spherical_coordinates;

  protected:
    //! Computes the transformation at an index
    vector<Point<double>> computeTransformation(int);

    //! Computes the spherical coordinate values
    array<double,3> computeSphericalValues(vector<Point<double>> &);

  public: 
    //! Constructor
    Protein(ProteinStructure *);

    //! Gets the list of spherical coordinates
    void getSphericalCoordinatesList();

};

#endif

