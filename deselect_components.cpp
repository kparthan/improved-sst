#include <iostream>
#include <memory>
#include <cstdlib>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

using namespace std;
using namespace boost::program_options;

struct Parameters
{
  string components_file;
  int remove;
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

  cout << "Checking command-line input ..." << endl;
  options_description desc("Allowed options");
  desc.add_options()
       ("help","produce help component")
       ("components_file",value<string>(&parameters.components_file),"path to the file")
       ("remove",value<int>(&parameters.remove),"# of components to remove")
  ;
  variables_map vm;
  store(parse_command_line(argc,argv,desc),vm);
  notify(vm);
}

vector<string> removeComponents(struct Parameters &parameters)
{
  vector<string> all_lines;
  string line;
  ifstream file(parametersc.components_file.c_str());
  // read all lines
  while (getline(file,line)) {
    all_lines.push_back(line);
  } 
  file.close();
  // randomly select components to remove
  int num_lines = all_lines.size();
  vector<int> index(num_lines,0);
  for (int i=0; i<parameters.remove; i++) {
    
  }
  return lines;
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  vector<string> lines = removeComponents(parameters);
  return 0;
}
/*
void Mixture::load(string &file_name)
{
  sample_size.clear();
  weights.clear();
  components.clear();
  K = 0;
  ifstream file(file_name.c_str());
  string line;
  vector<double> numbers;
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
    double kappa = numbers[4];
    Component component(mu,kappa,constrain_kappa);
    components.push_back(component);
    numbers.clear();
  }
  file.close();
}
*/
