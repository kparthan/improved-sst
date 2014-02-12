#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include "Header.h"

class Segmentation
{
  private:
    //! the structure
    ProteinStructure *structure;

    //! the chain
    string chain_id;

    //! the internal segmentation
    vector<array<int,2>> segments;

    //! corresponding ideal models
    vector<string> model_names;

  public:
    //! Null constructor
    Segmentation();

    //! Constructor
    Segmentation(ProteinStructure *, string, vector<int> &, vector<string> &);

    //! Gets the strand segments
    vector<array<int,2>> getStrandSegments();

    //! Postprocess the strand assignment
    void postprocess();
};

#endif

