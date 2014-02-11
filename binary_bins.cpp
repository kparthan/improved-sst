#include <iostream>
#include <memory>
#include <cstdlib>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <iomanip>

using namespace std;
using namespace std::chrono;
using namespace boost::program_options;

struct Parameters
{
  string bins_file;
  //double res;
};

void Usage(const char *exe, options_description &desc)
{
  cout << "Usage: " << exe << " [options]" << endl;
  cout << desc << endl;
  exit(1);
}

struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;

  options_description desc("Allowed options");
  desc.add_options()
       ("help","produce help component")
       ("bins",value<string>(&parameters.bins_file),"path to the file")
       //("res",value<double>(&parameters.res),"resoultion")
  ;
  variables_map vm;
  store(parse_command_line(argc,argv,desc),vm);
  notify(vm);

  return parameters;
}

void spherical2cartesian(vector<double> &spherical, vector<double> &cartesian)
{
  cartesian[0] = spherical[0] * sin(spherical[1]) * cos(spherical[2]);
  cartesian[1] = spherical[0] * sin(spherical[1]) * sin(spherical[2]);
  cartesian[2] = spherical[0] * cos(spherical[1]);
}

void collectData(double res, vector<vector<int>> &binary_bins)
{
  vector<double> spherical(3,1);
  vector<double> cartesian(3,0);
  double theta = 0;
  ofstream file("coil.data");
  for (int i=0; i<binary_bins.size(); i++) {
    double phi = 0;
    for (int j=0; j<binary_bins[i].size(); j++) {
      assert(binary_bins[i][j] == 0 || binary_bins[i][j] == 1);
      if (binary_bins[i][j] == 1) {
        spherical[1] = theta * M_PI / 180;  // in radians
        spherical[2] = phi * M_PI / 180;    // in radians
        spherical2cartesian(spherical,cartesian);
        file << "A\t";
        for (int k=0; k<3; k++) {
          file << fixed << setw(10) << setprecision(5) << spherical[k] << "\t";
        }
        file << endl;
      }
      phi += res;
    }
    theta += res;
  }
  file.close();
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);
  ifstream bins(parameters.bins_file.c_str());
  string line;
  vector<vector<int>> binary_bins;
  vector<int> numbers;
  while (getline(bins,line)) {
    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    BOOST_FOREACH (const string& t, tokens) {
      istringstream iss(t);
      int x;
      iss >> x;
      numbers.push_back(x);
    }
    binary_bins.push_back(numbers);
    numbers.clear();
  }
  bins.close();
  ofstream binary_bins_file("binary_bins_2D");
  int total_bins = 0;
  int pop_bins = 0;
  for (int i=0; i<binary_bins.size(); i++) {
    for (int j=0; j<binary_bins[i].size(); j++) {
      assert(binary_bins[i][j] >= 0);
      total_bins++;
      if (binary_bins[i][j] > 0) {
        binary_bins[i][j] = 1;
        pop_bins++;
      }
      binary_bins_file << binary_bins[i][j] << "\t";
    }
    binary_bins_file << endl;
  }
  binary_bins_file.close();
  double res = 360 / (double)binary_bins[0].size();
  cout << "res: " << res << endl;
  cout << "total bins: " << total_bins << endl;
  cout << "pop bins: " << pop_bins << endl;
  collectData(res,binary_bins);
  return 0;
}

