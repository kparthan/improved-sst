#ifndef MIXTURE_H
#define MIXTURE_H

#include "Component.h"

class Mixture
{
  private:
    //! number of components
    int M;

    //! list of components
    vector<Component> components;
    
    //! weights of the components
    vector<double> weights;

  public:
    //! Null constructor
    Mixture();

    //! Constructor
    Mixture(int);
};

#endif

