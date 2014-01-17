#ifndef NORMAL_H
#define NORMAL_H

#include "Header.h"

class Normal
{
  private:
    //! Mean of the distribution
		double mu;

    //! Standard deviation of the distribution
		double sigma;

  public:
		//! Constructor
		Normal() ;

		//! Constructor that sets value of parameters
		Normal(double, double);

    //! Assignment of an existing Normal distribution
    Normal operator=(const Normal &);

		//! Gets the mean 
		const double mean();

    //! Gets the standard deviation
    const double standardDeviation(); 

		//! Function value
		double value(double);

    //! Computes the negative log likelihood of a sample
    double negativeLogLikelihood(vector<double> &);
};

#endif

