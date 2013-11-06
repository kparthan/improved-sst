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
    
    //! Sample (x_i) -- spherical coordinates with unit radius
    vector<array<double,2>> angles;

    //! Sample (x_i) -- Cartesian coordinates
    vector<array<double,3>> data;

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

    //! List of message lengths over several iterations
    vector<double> msglens;

  public:
    //! Null constructor
    Mixture();

    //! Constructor
    Mixture(int, vector<array<double,2>> &, int);

    //! Initialize parameters
    void initialize();

    //! Initialize parameters
    void initialize2();

    //! Initializes the parameters of the components
    void initializeComponentParameters();

    //! Updates the effective sample size
    void updateEffectiveSampleSize();

    //! Update alphas
    void updateAlphas();

    //! Update the component weights
    void updateWeights();

    //! Update components
    void updateComponents();

    //! Update the responsibility matrix
    void updateResponsibilityMatrix();

    //! Probability of a datum
    double probability(array<double,2> &);

    //! Computes the minimum message length
    double computeMinimumMessageLength();

    //! Estimate mixture parameters
    double estimateParameters();

    //! Prints the model parameters
    void printParameters(ostream &, int, double);

    //! Plot the variation in message length
    void plotMessageLengthEM();

    //! Loads the mixture file
    void load(string &);
};

#endif

