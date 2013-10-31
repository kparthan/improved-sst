#ifndef MESSAGE_H
#define MESSAGE_H

#include "Header.h"

class Message
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

    //! ML and MML estimates of kappa
    double kappa_ml,kappa_mml;

  public:
    //! Null constructor
    Message();

    //! Constructor
    Message(array<double,3> &, double);

    //! Minimize message length
    void minimize();

    //! Estimate MML mean parameter
    void estimateVonMisesMean();

    //! Estimate ML kappa
    double estimateKappa_ML();
 
    //! Estimate MML kappa
    double estimateKappa_MML(double);

    //! Computes the message length
    double computeMessageLength(double);

    //! Computes the first part of the message
    double computeFirstPart(double);

    //! Computes the second part of the message
    double computeSecondPart(double);

    //! Computes the rate of change of message length
    double computeFirstDerivative(double);

    //! Computes the second derivative of message length
    double computeSecondDerivative(double);
};

#endif

