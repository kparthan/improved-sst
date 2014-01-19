#ifndef IDEAL_MODEL_H
#define IDEAL_MODEL_H

#include "Header.h"

class IdealModel
{
  private:
    //! Length of the model
    int length;

    //! Name of the model
    string name;

  public:
    //! Null constructor
    IdealModel();
  
    //! Constructor
    IdealModel(int, string);

    //! Overloading = operator
    IdealModel operator=(const IdealModel &);

};

#endif

