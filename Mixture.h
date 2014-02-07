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
    
    //! Sample (x_i) -- Cartesian coordinates
    //! (on a sphere of unit radius)
    vector<vector<double>> data;

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

    //! Simulation flag
    int simulation;

    //! DSSP flag
    int dssp;

    //! Type of sst model
    string dssp_sst_type;

    //! Alphas (alpha_k)
    vector<double> alphas;

    //! List of message lengths over several iterations
    vector<double> msglens;

    //! Null model message length
    double null_msglen;

  public:
    //! Flag to bound kappa
    int constrain_kappa;

    //! Null constructor
    Mixture();

    //! Constructor
    Mixture(int, vector<Component> &, vector<double> &);

    //! Constructor
    Mixture(int, vector<vector<double>> &, int, int, int);

    //! Sets the DSSP flag
    void setDSSPFlag(string &);

    //! Gets the DSSP type
    string getDSSPType();

    //! Gets the list of weights
    vector<double> getWeights();

    //! Gets the list of components
    vector<Component> getComponents();

    //! Gets the sample size
    double getSampleSize();

    //! Sets the simulation flag 
    void setSimulationFlag();

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
    double probability(vector<double> &);

    //! Computes the negative log likelihood
    double negativeLogLikelihood(vector<vector<double>> &);

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
    void saveComponentData(int, vector<vector<double>> &);

    //! Generate random data using arbitrary sample size
    vector<vector<double>> generateRandomSampleSize(bool);

    //! Generate random data from the distribution using mixture proportions
    vector<vector<double>> generateProportionally(int, bool);

    //! Generate heat map data
    void generateHeatmapData(double);

    //! Conflates a component with the current mixture
    Mixture conflate(Component &);

    //! Gets the residual mixture
    Mixture getResidualMixture(vector<int> &, double &);
};

#endif

