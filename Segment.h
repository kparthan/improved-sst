#ifndef SEGMENT_H
#define SEGMENT_H

#include "Header.h"
#include "OptimalFit.h"

class Segment
{
  private:
    //! End indices of the segment
    int start,end;

    //! Cartesian coordinates of the segment
    vector<Point<double>> cartesian;

    //! Spherical coordinates of the segment
    vector<array<double,3>> spherical;

    //! Initial two distances
    double radii[2];

  public:
    //! Constructor
    Segment(int, int, vector<Point<double>> &, vector<array<double,3>> &);

    //! Set initial two distances
    void setInitialDistances(double, double);

    //! Fit null model
    OptimalFit fitNullModel();

    //! Fit an ideal model
    OptimalFit fitIdealModel(IdealModel &);
};

#endif

