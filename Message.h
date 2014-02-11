#ifndef MESSAGE_H
#define MESSAGE_H

#include "Normal.h"
#include "Poisson.h"
#include "Mixture.h"

class Message
{
  private:

  public:
    //! Null constructor
    Message();

    //! Computes the message length of stating an integer using the
    //! log star distribution
    double encodeUsingLogStarModel(double);

    //! Poisson distribution
    double encodeUsingPoissonModel(int, Poisson &);

    //! Computes the message length of stating the length using the
    //! normal model
    double encodeUsingNormalModel(double, Normal &);
 
    //! Computes the message length of stating the length using the
    //! sphere model
    double encodeUsingSphereModel(double, Normal &);

    //! Computes the message length of stating the direction using a
    //! mixture model
    double encodeUsingMixtureModel(vector<double> &, Mixture &);

    //! Computes the message length of stating the direction using a
    //! single component
    double encodeUsingComponent(vector<double> &, Component &);

};

#endif

