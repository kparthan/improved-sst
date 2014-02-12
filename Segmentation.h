#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include "Header.h"

class Segmentation
{
  private:
    //! the structure
    ProteinStructure *structure;

    //! the internal segmentation
    vector<int> segments;

    //! corresponding ideal models
    vector<string> model_names;

  public:
    //! Null constructor
    Segmentation();

    //! Constructor
    Segmentation(ProteinStructure *, vector<int> &, vector<string> &);
};

#endif

