#ifndef VON_MISES_3D_H
#define VON_MISES_3D_H

#include "Header.h"

class VonMises3D
{
  private:
    //! Cartesian mean of the distribution
    vector<double> unit_mean;

    //! Scale parameter of the distribution
		double kappa;

    //! Precomputed normalization constant
    double norm_constant;

    //! Precompute normalization constant
    void computeNormalizationConstant();

  public:
		//! Constructor
		VonMises3D() ;

		//! Constructor that sets value of parameters
		VonMises3D(vector<double> &, double);

    //! Assignment of an existing VonMises3D distribution
    VonMises3D operator=(const VonMises3D &);

		//! Gets the mean 
		const vector<double> mean();

    //! Gets the standard deviation
    const double scale(); 

		//! Function value
		double density(vector<double> &);

    //! Generate random samples
    vector<vector<double>> generateCanonical(int);

    //! Generate random samples
    vector<vector<double>> generateAngles(int);

    //! Generate random samples
    vector<vector<double>> generateCoordinates(int);

    //! Generate transformation matrix
    vector<vector<double>> constructRotationMatrix();
};

#endif

