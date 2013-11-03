#ifndef MIXTURE_H
#define MIXTURE_H

#include "Component.h"

class Mixture
{
  private:
    //! Number of components
    int K;

    //! List of components
    vector<Component> components;
    
    //! Sample (x_i)
    vector<array<double,2>> angles;

    //! Sample size
    int N;

    //! Responsibility matrix (r_ik)
    vector<vector<double>> responsibility;

    //! Effective sample size for each component (n_k)
    vector<double> sample_size;

    //! Weights of the components (a_k)
    vector<double> weights;

    //! Flag to use modified weight update rule
    int update_weights_new;

    //! Alphas (alpha_k)
    vector<double> alphas;

  public:
    //! Null constructor
    Mixture();

    //! Constructor
    Mixture(int, vector<array<double,2>> &, int);

    //! Initialize parameters
    void initialize();

    //! Updates the effective sample size
    void updateEffectiveSampleSize();

    //! Estimate mixture parameters
    void estimateParameters();
};

#endif

