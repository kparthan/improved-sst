#ifndef POISSON_H 
#define POISSON_H

#include "Header.h"

class Poisson
{
  private:
    //! Mean of the distribution
		double mu;

  public:
		//! Constructor
		Poisson() ;

		//! Constructor that sets value of parameters
		Poisson(double);

    //! Assignment of an existing Poisson distribution
    Poisson operator=(const Poisson &);

		//! Gets the mean 
		const double mean();

		//! Function value
		double density(int);
};

#endif

