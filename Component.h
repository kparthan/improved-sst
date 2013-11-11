#ifndef COMPONENT_H
#define COMPONENT_H

#include "Header.h"
#include "VonMises3D.h"

class Component
{
  private:
    //! Unit mean direction (Cartesian)
    array<double,3> mean_direction;

    //! Mean direction on a unit sphere
    array<double,2> mu;

    //! Sample size
    double N;

    //! rbar = |R|/N
    double rbar,R;

    //! Flag to constrain kappa
    int constrain_kappa;

    //! ML and MML estimates of kappa
    double kappa_ml,kappa_mml;

    //! Representative von mises distribution
    VonMises3D von_mises;

  public:
    //! Null constructor
    Component();

    //! Constructor
    Component(array<double,3> &, double, int);

    //! Constructor
    Component(array<double,2> &, double, int);

    //! Minimize message length
    void minimizeMessageLength();

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
    double likelihood(array<double,2> &);

    //! Computes the likelihood value
    double likelihood(array<double,2> &, double);

    //! Computes the probability of the component parameters
    double computeParametersProbability();

    //! Computes the prior density on the parameters
    double  computePriorDensity();

    //! Computes the expected Fisher information associated with the parameters
    double computeFisherInformation();

    //! Prints the component parameters 
    void printParameters(ostream &);

    //! Return theta
    double getTheta();

    //! Return phi 
    double getPhi();

    //! Return kappa 
    double getKappa();

    //! Generate random samples from this component
    vector<array<double,3>> generate(int);
};

#endif

