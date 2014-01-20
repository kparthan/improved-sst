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
  int profiles;
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
       ("profiles",value<int>(&parameters.profiles),
          "# of random profiles to be selected from each class")
  ;
  variables_map vm;
  store(parse_command_line(argc,argv,desc),vm);
  notify(vm);

  return parameters;
}

vector<string> readProfiles(string path)
{
  ifstream file(path.c_str());
  vector<string> all_lines;
  string line;
  while (getline(file,line)) {
    all_lines.push_back(line);
  } 
  file.close();
  return all_lines;
}

vector<string> selectProfiles(vector<string> &all_profiles, int num_profiles)
{
  int num_all_profiles = all_profiles.size();
  vector<int> select(num_all_profiles,0);
  vector<string> select_profiles;
  int count = 0;
  srand(time(NULL));
  while (count < num_profiles) {
    int index = rand() % num_profiles;
    if (select[index] == 0) {
      select[index] = 1;
      count++;
      select_profiles.push_back(all_profiles[index]);
    }
  }
  return select_profiles;
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  vector<string> paths;
  string cl = "./spherical_system/class-a-profiles";
  paths.push_back(cl);
  cl = "./spherical_system/class-b-profiles";
  paths.push_back(cl);
  cl = "./spherical_system/class-c-profiles";
  paths.push_back(cl);
   
  vector<string> profiles;
  for (int i=0; i<paths.size(); i++) {
    // read profiles of a class
    vector<string> all_profiles = readProfiles(paths[i]);
    cout << paths[i] << "\t" << all_profiles.size() << "\t";
    // select random number of profiles
    vector<string> 
    select_profiles = selectProfiles(all_profiles,parameters.profiles);
    cout << select_profiles.size() << endl;
    assert(select_profiles.size() == parameters.profiles);
    // store the selected profiles
    for (int j=0; j<parameters.profiles; j++) {
      profiles.push_back(select_profiles[j]);
    }
  }

  string cmd = "cp ./spherical_system/profiles/";
  for (int i=0; i<profiles.size(); i++) {
    string exec = cmd + profiles[i] + " ./spherical_system/random_profiles/";
    system(exec.c_str());
  }


  return 0;
}

