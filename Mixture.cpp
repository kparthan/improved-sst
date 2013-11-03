#include "Mixture.h"

/*!
 *  \brief Null constructor module
 */
Mixture::Mixture()
{}

/*!
 *  \brief This is a constructor function.
 *  \param num_components an integer
 *  \param angles a reference to a vector<array<double,2>>
 *  \param update_weights_new an integer
 */
Mixture::Mixture(int num_components, vector<array<double,2>> &angles,
                 int update_weights_new) : K(num_components), angles(angles),
                 update_weights_new(update_weights_new)
{}

/*!
 *  \brief This function initializes the parameters of the model.
 */
void Mixture::initialize()
{
  N = angles.size();
  cout << "sample size: " << N << endl;
  // initialize weights of components
  /*double w = 1 / (double) K;
  for (int i=0; i<K; i++) {
    weights.push_back(w);
  }*/

  // initialize responsibility matrix
  srand(time(NULL));
  for (int i=0; i<N; i++) {
    vector<double> tmp(K,0);
    responsibility.push_back(tmp);
    int index = rand() % K;
    //cout << index << endl;
    responsibility[i][index] = 1;
  }
  sample_size = vector<double>(K,0);
  updateEffectiveSampleSize();
  weights = vector<double>(K,0);
  alphas = vector<double>(K,0);
  updateWeights();

  // initialize parameters of each component
  double theta,phi,kappa;
  for (int i=0; i<K; i++) {
    // generate theta
    theta = (rand() / (double) RAND_MAX) * PI;
    // generate phi
    phi = (rand() / (double) RAND_MAX) * 2 * PI;
    // generate kappa
    kappa = tan((rand() / (double) RAND_MAX) * (PI/2));
    array<double,3> coords = convertToCartesian(1,theta,phi);
    Component component(coords);
  }
}

/*!
 *  \brief This function updates the effective sample size of each component.
 */
void Mixture::updateEffectiveSampleSize()
{
  for (int i=0; i<K; i++) {
    double count = 0;
    for (int j=0; j<N; j++) {
      count += responsibility[j][i];
    }
    sample_size[i] = count;
  }
  /*for (int i=0; i<M; i++) {
    cout << sample_size[i] << " ";
  }cout << endl;*/
}

/*!
 *  \brief This function is used to update the weights of the components.
 */
void Mixture::updateWeights()
{
  if (update_weights_new == UNSET) {
    double normalization_constant = N + (K/2.0);
    for (int i=0; i<K; i++) {
      weights[i] = (sample_size[i] + 0.5) / normalization_constant;
    }
  } else if (update_weights_new == SET) {
    updateAlphas();
    double A = 0;
    for (int i=0; i<K; i++) {
      A += alphas[i];
    }
    double normalization_constant = N + A - (K/2.0);
    for (int i=0; i<K; i++) {
      weights[i] = alphas[i] / normalization_constant;
    }
  }
}

/*!
 *  \brief This function is used to update the alpha parameters used in the
 *  modified update rule for the weights.
 */
void Mixture::updateAlphas()
{
  for (int i=0; i<K; i++) {
    double count = 0;
    for (int j=0; j<N; j++) {
      count += responsibility[j][i] * responsibility[j][i]
    }
    alphas[i] = count / (double) sample_size[i];
  }
}

/*!
 *  \brief This function is used to estimate the model .parameters.
 */
void Mixture::estimateParameters()
{
  initialize();
  //while (1) {
  //}
}

