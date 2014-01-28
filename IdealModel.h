#ifndef IDEAL_MODEL_H
#define IDEAL_MODEL_H

#include "Header.h"

class IdealModel
{
  private:
    //! The ideal ProteinStructure
    ProteinStructure *model;

    //! Residues of the model
    vector<vector<double>> residues;

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
    vector<vector<double>> getResidues(int);

    //! Returns the protein structure
    ProteinStructure *getStructure();

    //! Returns the name of the ideal model
    string getName();

    //! Prints the details of the optimal fit
    void printModelInfo();

};

#endif

