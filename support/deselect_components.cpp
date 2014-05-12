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

using namespace std;
using namespace std::chrono;
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

  options_description desc("Allowed options");
  desc.add_options()
       ("help","produce help component")
       ("components_file",value<string>(&parameters.components_file),"path to the file")
       ("remove",value<int>(&parameters.remove),"# of components to remove")
  ;
  variables_map vm;
  store(parse_command_line(argc,argv,desc),vm);
  notify(vm);

  return parameters;
}

vector<string> removeComponents(struct Parameters &parameters)
{
  vector<string> all_lines;
  string line;
  ifstream file(parameters.components_file.c_str());
  // read all lines
  while (getline(file,line)) {
    all_lines.push_back(line);
  } 
  file.close();
  // randomly select components to remove
  int num_lines = all_lines.size();
  vector<int> index(num_lines,0);
  auto ts = high_resolution_clock::now();
  usleep(10);
  auto te = high_resolution_clock::now();
  double t = duration_cast<nanoseconds>(ts-te).count();
  srand(t);
  for (int i=0; i<parameters.remove; i++) {
    int random_index = rand() % num_lines;
    if (index[random_index] == 0) {
      index[random_index] = 1;
    } else {
      i--;
    }
  }
  // retain components not selected to remove
  cout << "Components removed are:\n";
  vector<string> lines;
  for (int i=0; i<all_lines.size(); i++) {
    if (index[i] == 0) {
      lines.push_back(all_lines[i]);
    } else {
      cout << "[" << i+1 << "]\t";
      cout << all_lines[i] << endl;
    }
  }
  return lines;
}

void createFile(vector<string> &lines)
{
  ofstream initial("mixture/initialize_components_file");
  for (int i=0; i<lines.size(); i++) {
    initial << lines[i] << endl;
  }
  initial.close();
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  vector<string> lines = removeComponents(parameters);

  createFile(lines);

  return 0;
}

