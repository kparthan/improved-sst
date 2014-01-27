#include "Mixture.h"
#include "Support.h"

extern int initialize_components_from_file;
extern string CURRENT_DIRECTORY;

/*!
 *  \brief Null constructor module
 */
Mixture::Mixture()
{}

/*!
 *  \brief This is a constructor function which instantiates a Mixture
 *  \param K an integer
 *  \param components a reference to a vector<Component>
 *  \param weights a reference to a vector<double> 
 */
Mixture::Mixture(int K, vector<Component> &components, vector<double> &weights):
                 K(K), components(components), weights(weights)
{
  assert(components.size() == K);
  assert(weights.size() == K);
}

/*!
 *  \brief This is a constructor function.
 *  \param num_components an integer
 *  \param data a reference to a vector<vector<double>>
 *  \param update_weights_new an integer
 *  \param constrain_kappa an integer
 *  \param simulation an integer
 */
Mixture::Mixture(int num_components, vector<vector<double>> &data,
                 int update_weights_new, int constrain_kappa, int simulation) : 
                 K(num_components), data(data), simulation(simulation),
                 update_weights_new(update_weights_new),
                 constrain_kappa(constrain_kappa)
{}

/*!
 *  \brief This function sets the simulation flag.
 */
void Mixture::setSimulationFlag()
{
  simulation = SET;
}

/*!
 *  \brief This function initializes the parameters of the model.
 */
void Mixture::initialize()
{
  N = data.size();
  cout << "sample size: " << N << endl;
  // initialize weights of components
  /*double w = 1 / (double) K;
  for (int i=0; i<K; i++) {
    weights.push_back(w);
  }*/

  if (initialize_components_from_file == UNSET) {
    // initialize responsibility matrix
    auto ts = high_resolution_clock::now();
    usleep(1);
    auto te = high_resolution_clock::now();
    double t = duration_cast<nanoseconds>(ts-te).count();
    srand(t);
    vector<double> tmp(K,0);
    for (int i=0; i<N; i++) {
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
  } else if (initialize_components_from_file == SET) {
    vector<double> tmp(K,0);
    for (int i=0; i<N; i++) {
      responsibility.push_back(tmp);
      //int index = rand() % K;
      //cout << index << endl;
      //responsibility[i][index] = 1;
    }
    int count = 0;
    double w = 0;
    ifstream file("mixture/initialize_components_file");
    string line;
    vector<double> numbers;
    vector<double> spherical_mean(3,1);
    vector<double> unit_mean(3,0);
    while (getline(file,line)) {
      count++;
      boost::char_separator<char> sep("mukap,:()[] \t");
      boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
      BOOST_FOREACH (const string& t, tokens) {
        istringstream iss(t);
        double x;
        iss >> x;
        numbers.push_back(x);
      }
      sample_size.push_back(numbers[0]);
      weights.push_back(numbers[1]);
      w += numbers[1];
      array<double,2> mu({numbers[2],numbers[3]});
      angleInRadians(mu[0]);
      angleInRadians(mu[1]);
      spherical_mean[1] = mu[0]; spherical_mean[2] = mu[1];
      spherical2cartesian(spherical_mean,unit_mean);
      double kappa = numbers[4];
      Component component(unit_mean,kappa);
      components.push_back(component);
      numbers.clear();
    }
    file.close();
    if (count < K) {
      int n = K - count;
      vector<double> ws = generateRandomWeights(n,1-w);
      vector<Component> cs = generateRandomComponents(n,constrain_kappa);
      for (int i=0; i<n; i++) {
        weights.push_back(ws[i]);
        components.push_back(cs[i]);
        sample_size.push_back(0);
      }
    }
  }
}

/*!
 *  \brief This function initializes the parameters of the model.
 */
void Mixture::initialize2()
{
  N = data.size();
  cout << "sample size: " << N << endl;

  // initialize weights of components
  weights = generateRandomWeights(K,1);

  // initialize components
  components = generateRandomComponents(K,constrain_kappa);

  // set the dimensions of r_{ik}, n_k and alpha_k matrices
  vector<double> tmp(K,0);
  for (int i=0; i<N; i++) {
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
    vector<double> mean(3,0);
    for (int j=0; j<N; j++) {
      if (responsibility[j][i] > 0.5) {
        vector<double> x = data[j]; 
        for (int k=0; k<3; k++) {
          mean[k] += x[k];
        }
      }
    }
    Component component(mean,sample_size[i],constrain_kappa);
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
      weights[i] = (sample_size[i] + alphas[i] - 0.5) / normalization_constant;
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
    vector<double> sum(3,0);
    for (int j=0; j<N; j++) {
      for (int k=0; k<3; k++) {
        sum[k] += responsibility[j][i] * data[j][k];
      }
    }
    /*for (int k=0; k<3; k++) {
      sum[k] /= sample_size[i];
    }*/
    Component component(sum,sample_size[i],constrain_kappa);
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
      density[j] = components[j].likelihood(data[i]);
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
double Mixture::probability(vector<double> &x)
{
  double px = 0;
  vector<double> density(K,0);
  for (int i=0; i<K; i++) {
    density[i] = components[i].likelihood(x);
    px += weights[i] * density[i];
  }
  if (!(px > 0)) {
    cout << "px: " << px << endl;
    for (int i=0; i<K; i++) {
      cout << weights[i] << "\t" << density[i] << endl;
    }
  }
  fflush(stdout);
  assert(px >= 0);
  return px;
}

/*!
 *  \brief This function computes the negative log likelihood of a data
 *  sample.
 *  \param a reference to a vector<array<double,2>>
 *  \return the negative log likelihood (base e)
 */
double Mixture::negativeLogLikelihood(vector<vector<double>> &sample)
{
  double value = 0,density;
  for (int i=0; i<sample.size(); i++) {
    density = probability(sample[i]);
    //cout << i+1 << " " << -log2(density)-2*log2(AOM) << endl;
    value -= log(density);
  }
  return value;
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
  assert(Ik > 0);

  // enocde the weights
  double Iw = ((K-1)/2.0) * log(N);
  Iw -= log(boost::math::factorial<long double>(K-1));
  for (int i=0; i<K; i++) {
    Iw -= 0.5 * log(weights[i]);
  }
  cout << "Iw: " << Iw << endl;
  assert(Iw >= 0);

  // encode the likelihood of the sample
  double Il = 0,px;
  for (int i=0; i<N; i++) {
    px = probability(data[i]);
    Il -= log(px);
  }
  Il -= 2 * N * log(AOM);
  cout << "Il: " << Il << endl;
  assert(Il > 0);

  // encode the parameters of the components
  double It = 0,p;
  for (int i=0; i<K; i++) {
    p = components[i].computeParametersProbability();
    It -= log(p);
  }
  cout << "It: " << It << endl;
  assert(It > 0);

  // the constant term
  // # of continuous parameters d = 4K-1
  double cd = computeConstantTerm(4*K-1);
  cout << "cd: " << cd << endl;

  return (Ik + Iw + Il + It + cd)/(log(2));
}

/*!
 *  \brief This function is used to estimate the model parameters by running
 *  an EM algorithm.
 *  \return the stable message length
 */
double Mixture::estimateParameters()
{

  clock_t c_start = clock();
  auto t_start = high_resolution_clock::now();

  /* prepare log file */
  string file_name = string(CURRENT_DIRECTORY) + "/mixture/";
  if (simulation == SET) {
    file_name += "simulation/";
  }
  file_name += "logs/";
  if (simulation == UNSET) {
    if (update_weights_new == UNSET) {
      file_name += "normal_weights_update/";
    } else if (update_weights_new == SET){
      file_name += "new_weights_update/";
    }
  }
  file_name += boost::lexical_cast<string>(K) + ".log";
  ofstream log(file_name.c_str());
  

  /* set the max allowed diff in msglen */
  double MAX_ALLOWED_DIFF_MSGLEN;
  if (simulation == SET) {
    // terminate EM when diff in msglen is < 1 bit
    MAX_ALLOWED_DIFF_MSGLEN = 1; 
  } else {
    // terminate EM when diff in msglen is < 10 bits
    MAX_ALLOWED_DIFF_MSGLEN = 10;
  }

  /* EM loop */
  double prev=0,current;
  int iter = 1;
  initialize();
  printParameters(log,0,current);
  null_msglen = computeNullModelMessageLength();
  cout << "null_msglen: " << null_msglen << endl;
  //MAX_ALLOWED_DIFF_MSGLEN = 0.001 * null_msglen;
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
    //if (iter == 30) break;
    if (iter != 1) {
      assert(current > 0);
      //assert(current <= prev);
      // because EM has to consistently produce lower 
      // message lengths otherwise something wrong!
      if (iter > 20 && current <= prev && (prev - current) <= 0.0001 * prev) {
      //if (iter > 20 && current <= prev && prev - current < MAX_ALLOWED_DIFF_MSGLEN) {
      //if (prev - current < 0.005 * prev) {  // if decrement is less than 0.5 %
                                              // terminates prematurely
        log << "\t\t" << current/N << " bpr" << endl;
        log << "Null model msglen: " << null_msglen << " bits.";
        log << "\t(" << null_msglen/N << " bpr)" << endl;
        break;
      }
    }
    prev = current;
    iter++;
  }
  log.close();
  plotMessageLengthEM();

  clock_t c_end = clock();
  auto t_end = high_resolution_clock::now();
  double cpu_time = double(c_end-c_start)/(double)(CLOCKS_PER_SEC);
  double wall_time = duration_cast<seconds>(t_end-t_start).count();
  // update summary file
  ofstream summary("summary",ios::app);
  summary << fixed << setw(5) << K;
  summary << fixed << setw(10) << iter;
  summary << fixed << setw(20) << setprecision(3) << cpu_time/60; // in mins
  summary << fixed << setw(20) << setprecision(3) << wall_time/60; // in mins
  summary << fixed << setw(20) << setprecision(3) << current << endl;
  summary.close();
  return current;
}

/*!
 *  \brief This function computes the null model message length.
 *  \return the null model message length
 */
double Mixture::computeNullModelMessageLength()
{
  null_msglen = 0;
  null_msglen += N * (log(4*PI)-(2*log(AOM)));
  /*for (int i=0; i<data.size(); i++) {
    double theta = angleInRadians(angles[i][0]);
    null_msglen -= log(sin(theta));
  }*/
  return null_msglen / log(2);
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
    os << "\t" << fixed << setw(10) << setprecision(5) << sample_size[k];
    os << "\t" << fixed << setw(10) << setprecision(5) << weights[k];
    os << "\t";
    components[k].printParameters(os);
  }
  os << "\t\t\tmsglen: " << msglen << " bits." << endl;
}

/*!
 *  \brief This function prints the parameters used to simulate the mixture
 *  modelling process.
 */
void Mixture::printParameters()
{
  ofstream file("mixture/simulation/simulate.dat");
  for (int k=0; k<K; k++) {
    file << "\t" << fixed << setw(10) << (int)sample_size[k];
    file << "\t" << fixed << setw(10) << setprecision(5) << weights[k];
    file << "\t";
    components[k].printParameters(file);
  }
  file.close();
}

/*!
 *  \brief This function is used to plot the variation in message length
 *  with each iteration.
 */
void Mixture::plotMessageLengthEM()
{
  // output the data to a file
  string data_file = string(CURRENT_DIRECTORY) + "/mixture/";
  string output_file = string(CURRENT_DIRECTORY) + "/mixture/";
  string script_file = string(CURRENT_DIRECTORY) + "/mixture/";
  if (simulation == SET) {
    data_file += "simulation/";
    output_file += "simulation/";
    script_file += "simulation/";
  }
  data_file += "msglens/";
  output_file += "plots/";
  script_file += "plots/";
  if (simulation == UNSET) {
    if (update_weights_new == UNSET) {
      data_file += "normal_weights_update/";
      output_file += "normal_weights_update/";
      script_file += "normal_weights_update/";
    } else if (update_weights_new == SET){
      data_file += "new_weights_update/";
      output_file += "new_weights_update/";
      script_file += "new_weights_update/";
    }
  }
  string num_comp = boost::lexical_cast<string>(K);
  data_file += num_comp + ".dat";
  output_file += num_comp + ".eps";
  script_file += num_comp + "_script.p";
  ofstream file(data_file.c_str());
  for (int i=0; i<msglens.size(); i++) {
    file << i << "\t" << msglens[i] << endl;
  }
  file.close();

  // prepare gnuplot script file
  ofstream script(script_file.c_str());
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
	script << "set output \"" << output_file << "\"" << endl ;
	script << "plot \"" << data_file << "\" using 1:2 notitle " 
         << "with linespoints lc rgb \"red\"" << endl ;
  script.close();
  string cmd = "gnuplot -persist " + script_file;
  system(cmd.c_str());
  cmd = "rm " + script_file;
  system(cmd.c_str());
}

/*!
 *  \brief This function is used to read the mixture details to aid in
 *  visualization.
 *  \param file_name a reference to a string
 */
void Mixture::load(string &file_name)
{
  sample_size.clear();
  weights.clear();
  components.clear();
  K = 0;
  ifstream file(file_name.c_str());
  string line;
  vector<double> numbers;
  vector<double> spherical_mean(3,1);
  vector<double> unit_mean(3,0);
  while (getline(file,line)) {
    K++;
    boost::char_separator<char> sep("mukap,:()[] \t");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    BOOST_FOREACH (const string& t, tokens) {
      istringstream iss(t);
      double x;
      iss >> x;
      numbers.push_back(x);
    }
    sample_size.push_back(numbers[0]);
    weights.push_back(numbers[1]);
    array<double,2> mu({numbers[2],numbers[3]});
    angleInRadians(mu[0]);
    angleInRadians(mu[1]);
    spherical_mean[1] = mu[0]; spherical_mean[2] = mu[1];
    spherical2cartesian(spherical_mean,unit_mean);
    double kappa = numbers[4];
    Component component(unit_mean,kappa);
    components.push_back(component);
    numbers.clear();
  }
  file.close();
}

/*!
 *  \brief This function is used to randomly choose a component.
 *  \return the component index
 */
int Mixture::randomComponent()
{
  auto ts = high_resolution_clock::now();
  usleep(10);
  auto te = high_resolution_clock::now();
  double t = duration_cast<nanoseconds>(ts-te).count();
  srand(t);
  double random = rand() / (double) RAND_MAX;
  //cout << random << endl;
  double previous = 0;
  for (int i=0; i<weights.size(); i++) {
    if (random <= weights[i] + previous) {
      return i;
    }
    previous += weights[i];
  }
}

/*!
 *  \brief This function saves the data generated from a component to a file.
 *  \param index an integer
 *  \param data a reference to a vector<array<double,3>>
 */
void Mixture::saveComponentData(int index, vector<vector<double>> &data)
{
  string data_file = string(CURRENT_DIRECTORY) + "/mixture/visualize/comp";
  data_file += boost::lexical_cast<string>(index+1) + ".dat";
  components[index].printParameters(cout);
  ofstream file(data_file.c_str());
  for (int j=0; j<data.size(); j++) {
    for (int k=0; k<3; k++) {
      file << fixed << setw(10) << setprecision(3) << data[j][k];
    }
    file << endl;
  }
  file.close();
}

/*!
 *  \brief This function is used to generate samples of arbitrary size from 
 *  each component.
 *  \param save_data a boolean variable
 */
vector<vector<double>> Mixture::generateRandomSampleSize(bool save_data)
{
  vector<vector<double>> sample;
  for (int i=0; i<K; i++) {
    vector<vector<double>> x = components[i].generate((int)sample_size[i]);
    if (save_data) {
      saveComponentData(i,x);
    }
    for (int j=0; j<x.size(); j++) {
      sample.push_back(x[j]);
    }
  }
  for (int i=0; i<sample_size.size(); i++) {
    cout << sample_size[i] << endl;
  }
  return sample;
}

/*!
 *  \brief This function is used to randomly sample from the mixture
 *  distribution.
 *  \param num_samples an integer
 *  \param save_data a boolean variable
 *  \return the random sample
 */
vector<vector<double>>
Mixture::generateProportionally(int num_samples, bool save_data) 
{
  sample_size = vector<double>(K,0);
  for (int i=0; i<num_samples; i++) {
    // randomly choose a component
    int k = randomComponent();
    sample_size[k]++;
  }
  vector<vector<double>> sample;
  for (int i=0; i<K; i++) {
    vector<vector<double>> x = components[i].generate((int)sample_size[i]);
    if (save_data) {
      saveComponentData(i,x);
    }
    for (int j=0; j<x.size(); j++) {
      sample.push_back(x[j]);
    }
  }
  return sample;
}


/*!
 *  \brief This function generates data to visualize the 2D/3D heat maps.
 *  \param res a double
 */
void Mixture::generateHeatmapData(double res)
{
  string data_fbins2D = string(CURRENT_DIRECTORY) + "/mixture/visualize/bins_2D.dat";
  string data_fbins3D = string(CURRENT_DIRECTORY) + "/mixture/visualize/bins_3D.dat";
  ofstream fbins2D(data_fbins2D.c_str());
  ofstream fbins3D(data_fbins3D.c_str());
  vector<double> x(3,1);
  vector<double> point(3,0);
  for (double theta=0; theta<180; theta+=res) {
    x[1] = theta * PI/180;
    for (double phi=0; phi<360; phi+=res) {
      x[2] = phi * PI/180;
      double pr = probability(x);
      // 2D bins
      fbins2D << fixed << setw(10) << setprecision(4) << pr;
      // 3D bins
      spherical2cartesian(x,point);
      for (int k=0; k<3; k++) {
        fbins3D << fixed << setw(10) << setprecision(4) << point[k];
      }
      fbins3D << fixed << setw(10) << setprecision(4) << pr << endl;
    }
    fbins2D << endl;
  }
  fbins2D.close();
  fbins3D.close();
}

/*!
 *  \brief This function conflates a component with the current mixture.
 *  \param component a reference to a Component 
 *  \return the conflated  Mixture
 */
Mixture Mixture::conflate(Component &component)
{
  vector<Component> new_components;
  for (int i=0; i<K; i++) {
    Component conflated_component = components[i].conflate(component);
    new_components.push_back(conflated_component);
  }
  return Mixture(K,new_components,weights);
}

