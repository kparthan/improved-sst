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
    //! (measured in degrees)
    vector<array<double,2>> angles;

    //! Sample (x_i) -- Cartesian coordinates
    //! (on a sphere of unit radius)
    vector<array<double,3>> data;

    //! Sample size
    int N;

    //! Flag to bound kappa
    int constrain_kappa;

    //! Responsibility matrix (r_ik)
    vector<vector<double>> responsibility;

    //! Effective sample size for each component (n_k)
    vector<double> sample_size;

    //! Weights of the components (a_k)
    vector<double> weights;

    //! Flag to use modified weight update rule
    int update_weights_new;

    //! Simulation flag
    int simulation;

    //! Alphas (alpha_k)
    vector<double> alphas;

    //! List of message lengths over several iterations
    vector<double> msglens;

    //! Null model message length
    double null_msglen;

  public:
    //! Null constructor
    Mixture();

    //! Constructor
    Mixture(int, vector<array<double,2>> &, int, int, int);

    //! Constructor
    Mixture(int, vector<array<double,3>> &, int, int, int);

    //! Constructor
    Mixture(int , vector<double> &, vector<Component> &);

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

    //! Computes the null model message length
    double computeNullModelMessageLength();

    //! Prints the model parameters
    void printParameters(ostream &, int, double);

    //! Prints the model parameters used in simulation
    void printParameters();

    //! Plot the variation in message length
    void plotMessageLengthEM();

    //! Loads the mixture file
    void load(string &);

    //! Randomly choose a component
    int randomComponent();

    //! Saves the data generated from a component
    void saveComponentData(int, vector<array<double,3>> &);

    //! Generate random data using arbitrary sample size
    vector<array<double,3>> generateRandomSampleSize(bool);

    //! Generate random data from the distribution using mixture proportions
    vector<array<double,3>> generateProportionally(int, bool);

    //! Generate heat map data
    void generateHeatmapData(double);
};

#endif

