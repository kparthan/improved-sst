#ifndef VON_MISES_3D_H
#define VON_MISES_3D_H

#include "Header.h"

class VonMises3D
{
  private:
    //! Mean of the distribution
		array<double,2> mu;

    //! Scale parameter of the distribution
		double kappa;

    //! Precomputed constants
    double c3k;
		array<double,3> kmu;

    //!
    void computeConstants();

  public:
		//! Constructor
		VonMises3D() ;

		//! Constructor that sets value of parameters
		VonMises3D(array<double,2> &, double);

    //! Assignment of an existing VonMises3D distribution
    VonMises3D operator=(const VonMises3D &);

		//! Gets the mean 
		const array<double,2> mean();

    //! Gets the standard deviation
    const double scale(); 

		//! Function value
		double density(double, double);

};

#endif

