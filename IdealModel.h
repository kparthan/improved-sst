#ifndef IDEAL_MODEL_H
#define IDEAL_MODEL_H

#include "Header.h"

class IdealModel
{
  private:
    //! The ideal ProteinStructure
    ProteinStructure *model;

    //! Residues of the model
    vector<Point<double>> residues;

    //! Length of the model
    int length;

    //! Name of the model
    string name;

  public:
    //! Null constructor
    IdealModel();
  
    //! Constructor
    IdealModel(int, string);

    //! Constructor
    IdealModel(ProteinStructure *, string);

    //! Alters the length of the ideal model
    void setLength(int);

    //! Overloading = operator
    IdealModel operator=(const IdealModel &);

    //! Gets the condensed list of residues
    vector<Point<double>> getResidues(int);

    //! Returns the name of the ideal model
    string getName();

};

#endif

