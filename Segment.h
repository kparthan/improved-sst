#ifndef SEGMENT_H
#define SEGMENT_H

#include "Support.h"
#include "OptimalFit.h"
#include "Mixture.h"

class Segment
{
  private:
    //! End indices of the segment
    int start,end;

    //! Length of the segment
    int num_residues;

    //! Cartesian coordinates of the protein 
    vector<Point<double>> cartesian;

    //! Cartesian coordinates of the segment
    vector<vector<double>> observed_residues;

    //! Spherical coordinates of the protein 
    vector<array<double,3>> spherical;

    //! Initial two distances
    vector<double> radii;

  public:
    //! Constructor
    Segment(int, int, vector<Point<double>> &, vector<array<double,3>> &);

    //! Set initial two distances
    void setInitialDistances(double, double);

    //! Fit null model
    OptimalFit fitNullModel(Mixture &);

    //! Fit an ideal model
    OptimalFit fitIdealModel(IdealModel &, Mixture &);

    //!
    void transform();
};

#endif

