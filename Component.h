#ifndef COMPONENT_H
#define COMPONENT_H

#include "Header.h"
#include "VonMises3D.h"

class Component
{
  private:
    //! Resultant mean direction (Cartesian)
    vector<double> mean_direction;

    //! Sample size
    double N;

    //! rbar = |R|/N
    double rbar,R;

    //! Flag to constrain kappa
    int constrain_kappa;

    //! ML/MML estimate of (unit) mean direction
    vector<double> unit_mean;

    //! Directions (theta,phi) in radians of the unit mean
    array<double,2> mu;

    //! ML and MML estimates of kappa
    double kappa_ml,kappa_mml;

    //! Representative von mises distribution
    VonMises3D von_mises;

  public:
    //! Null constructor
    Component();

    //! Constructor
    Component(vector<double> &, double, int);

    //! Constructor
    Component(vector<double> &, double);

    // Update mu
    void updateMu(vector<double> &);

    //! Overloading = operator
    Component operator=(const Component &);

    //! Minimize message length
    void minimizeMessageLength();

    //! Minimize message length
    void minimizeMessageLength(vector<double> &);

    //! Computes the resultant length
    void computeMeanResultantLength();

    //! Computes the kappa estimates
    void estimateKappas();

    //! Estimate MML mean parameter
    void estimateVonMisesMean();

    //! Estimate ML kappa
    double estimateKappa_ML();
 
    //! Estimate MML kappa
    double estimateKappa_MML(double);

    //! Computes the message length
    double computeComponentLength(double);

    //! Computes the first part of the message
    double computeFirstPart(double);

    //! Computes the second part of the message
    double computeSecondPart(double);

    //! Computes the rate of change of message length
    double computeFirstDerivative(double);

    //! Computes the second derivative of message length
    double computeSecondDerivative(double);

    //! Computes the likelihood value
    double likelihood(vector<double> &);

    //! Computes the likelihood value
    double likelihood(vector<double> &, double);

    //! Computes the probability of the component parameters
    double computeParametersProbability();

    //! Computes the prior density on the parameters
    double  computePriorDensity();

    //! Computes the expected Fisher information associated with the parameters
    double computeFisherInformation();

    //! Prints the component parameters 
    void printParameters(ostream &);

    //! Returns the mean direction
    vector<double> getMeanDirection();

    //! Return kappa 
    double getKappa();

    //! Generate random samples from this component
    vector<vector<double>> generate(int);

    //! Conflates two components
    Component conflate(Component &);
};

#endif

