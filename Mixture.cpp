#include "Mixture.h"
#include "Support.h"

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
{
  for (int i=0; i<angles.size(); i++) {
    array<double,3> x = convertToCartesian(1,angles[i][0],angles[i][1]);
    data.push_back(x);
  }
}

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
  initializeComponentParameters();
  updateResponsibilityMatrix();
}

/*!
 *  \brief This function initializes the parameters of the model.
 */
void Mixture::initialize2()
{
  N = angles.size();
  cout << "sample size: " << N << endl;
  srand(time(NULL));
  double w = 1 / (double) K;
  for (int i=0; i<K; i++) {
    // initialize weights of components
    weights.push_back(w);
    // initialize component parameters
    array<double,2> mu;
    mu[0] = (rand()/(double)RAND_MAX)*180;
    mu[1] = (rand()/(double)RAND_MAX)*360;
    double kappa = (rand() / (double) RAND_MAX) * 100;
    Component component(mu,kappa);
    components.push_back(component);
  }

  // set the dimensions of r_{ik}, n_k and alpha_k matrices
  for (int i=0; i<N; i++) {
    vector<double> tmp(K,0);
    responsibility.push_back(tmp);
  }
  sample_size = vector<double>(K,0);
  alphas = vector<double>(K,0);
}

/*!
 *  \brief This function initializes the mean direction and the kappa values
 *  of the individual components.
 */
void Mixture::initializeComponentParameters()
{
  for (int i=0; i<K; i++) {
    array<double,3> mean({0,0,0});
    for (int j=0; j<N; j++) {
      if (responsibility[j][i] == 1) {
        array<double,3> x = data[j]; 
        for (int k=0; k<3; k++) {
          mean[k] += x[k];
        }
      }
    }
    Component component(mean,sample_size[i]);
    component.minimizeMessageLength();
    components.push_back(component);
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
      count += responsibility[j][i] * responsibility[j][i];
    }
    alphas[i] = count / (double) sample_size[i];
  }
}

/*!
 *  \brief This function is used to update the components.
 */
void Mixture::updateComponents()
{
  components.clear();
  for (int i=0; i<K; i++) {
    array<double,3> sum({0,0,0});
    for (int j=0; j<N; j++) {
      for (int k=0; k<3; k++) {
        sum[k] += responsibility[j][i] * data[j][k];
      }
    }
    /*for (int k=0; k<3; k++) {
      sum[k] /= sample_size[i];
    }*/
    Component component(sum,sample_size[i]);
    component.minimizeMessageLength();
    components.push_back(component);
  }
}

/*!
 *  \brief This function updates the terms in the responsibility matrix.
 */
void Mixture::updateResponsibilityMatrix()
{
  double px;
  vector<double> density(K,0);
  for (int i=0; i<N; i++) {
    px = 0;
    for (int j=0; j<K; j++) {
      density[j] = components[j].likelihood(angles[i]);
      px += weights[j] * density[j]; 
    }
    for (int j=0; j<K; j++) {
      responsibility[i][j] = weights[j] * density[j] / px;
    }
  }
}

/*!
 *  \brief This function computes the probability of a datum from the mixture
 *  model.
 *  \param x a reference to an array<double,2>
 *  \return the probability value
 */
double Mixture::probability(array<double,2> &x)
{
  double px = 0,density;
  for (int i=0; i<K; i++) {
    density = components[i].likelihood(x);
    px += weights[i] * density;
  }
  return px;
}

/*!
 *  \brief This function computes the minimum message length using the current
 *  model parameters.
 *  \return the minimum message length
 */
double Mixture::computeMinimumMessageLength()
{
  // encode the number of components
  // assume uniform priors
  double Ik = log(MAX_COMPONENTS);
  cout << "Ik: " << Ik << endl;

  // enocde the weights
  double Iw = ((K-1)/2.0) * log(N);
  Iw -= log(boost::math::factorial<double>(K-1));
  for (int i=0; i<K; i++) {
    Iw -= 0.5 * log(weights[i]);
  }
  cout << "Iw: " << Iw << endl;

  // encode the likelihood of the sample
  double Il = 0;
  for (int i=0; i<N; i++) {
    double px = probability(angles[i]);
    Il -= log(px);
  }
  Il -= 2 * N * log(AOM);
  cout << "Il: " << Il << endl;

  // encode the parameters of the components
  double It = 0,p;
  for (int i=0; i<K; i++) {
    p = components[i].computeParametersProbability();
    It -= log(p);
  }
  cout << "It: " << It << endl;

  // the constant term
  // # of continuous parameters d = 4K-1
  double cd = computeConstantTerm(4*K-1);
  cout << "cd: " << cd << endl;

  return (Ik + Iw + Il + It + cd)/(log(2));
}

/*!
 *  \brief This function is used to estimate the model parameters by running
 *  an EM algorithm.
 */
void Mixture::estimateParameters()
{
  initialize2();
  double current;
  int iter = 1;
  string file_name = string(CURRENT_DIRECTORY) + "mixture/logs/";
  file_name += boost::lexical_cast<string>(K) + ".log";
  ofstream log(file_name.c_str());
  printParameters(log,0,current);
  while (1) {
    // Expectation (E-step)
    updateResponsibilityMatrix();
    updateEffectiveSampleSize();
    // Maximization (M-step)
    updateWeights();
    updateComponents();
    current = computeMinimumMessageLength();
    msglens.push_back(current);
    printParameters(log,iter,current);
    if (iter++ == 3) break;
  }
  log.close();
  //plotMessageLength();
}

/*!
 *  \brief This function prints the parameters to a log file.
 *  \param os a reference to a ostream
 *  \param iter an integer
 *  \param msglen a double
 */
void Mixture::printParameters(ostream &os, int iter, double msglen)
{
  os << "Iteration #: " << iter << endl;
  for (int k=0; k<K; k++) {
    os << "\t" << fixed << setw(10) << setprecision(3) << sample_size[k];
    os << "\t" << fixed << setw(10) << setprecision(3) << weights[k];
    os << "\t";
    components[k].printParameters(os);
  }
  os << "\t\t\tmsglen: " << msglen << " bits." << endl;
}

/*!
 *  \brief This function is used to plot the variation in message length
 *  with each iteration.
 */
void Mixture::plotMessageLength()
{
  // output the data to a file
  string data_file = string(CURRENT_DIRECTORY) + "mixture/msglens/";
  data_file += boost::lexical_cast<string>(K) + ".dat";
  ofstream file(data_file.c_str());
  for (int i=0; i<msglens.size(); i++) {
    file << i << "\t" << msglens[i] << endl;
  }
  file.close();

  // prepare gnuplot script file
  ofstream script("script.p");
	script << "# Gnuplot script file for plotting data in file \"data\"\n\n" ;
	script << "set terminal post eps" << endl ;
	script << "set autoscale\t" ;
	script << "# scale axes automatically" << endl ;
	script << "set xtic auto\t" ;
	script << "# set xtics automatically" << endl ;
	script << "set ytic auto\t" ;
	script << "# set ytics automatically" << endl ;
	script << "set title \"# of components: " << K << "\"" << endl ;
	script << "set xlabel \"# of iterations\"" << endl ;
	script << "set ylabel \"message length (in bits)\"" << endl ;
	script << "set output \"mixture/plots/" << K << ".eps\"" << endl ;
	script << "plot \"mixture/msglens/" << K << ".dat\" using 1:2 notitle " 
         << "with linespoints lc rgb \"red\"" << endl ;
  script.close();
  system("gnuplot -persist script.p") ;	
}

