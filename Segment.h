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

    //! Cartesian coordinates of the protein chain
    vector<vector<double>> cartesian_coordinates;

    //! Cartesian coordinates of the current segment
    vector<vector<double>> observed_residues;

    //! Spherical coordinates of the protein chain
    vector<vector<double>> spherical_coordinates;

    //! Unit coordinates of the protein chain
    vector<vector<double>> unit_coordinates;

    //! Initial two distances
    vector<double> radii;

  public:
    //! Constructor
    Segment(int, int, vector<vector<double>> &, vector<vector<double>> &,
            vector<vector<double>> &);

    //! Set initial two distances
    void setInitialDistances(double, double);

    //! Fit null model
    OptimalFit fitNullModel(Mixture &);

    //! Fit null model (for non-adaptive)
    OptimalFit fitNullModel(Mixture &, double &);

    //! Fit an ideal model using the adaptive model for mixture
    OptimalFit fitIdealModel(IdealModel &, Mixture &, int);

    //! Fit an ideal model using the non adaptive model for mixture
    OptimalFit fitIdealModel(IdealModel &, Mixture &, Component &, double);

    //! Computes the current mean and direction in the adaptive superposition
    void getCurrentMeanAndDirection(vector<double> &, vector<vector<double>> &,
    vector<double> &, vector<double> &, vector<double> &, vector<double> &, 
    int, vector<double> &, vector<double> &);

};

#endif

