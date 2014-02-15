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
#include <boost/filesystem.hpp>
#include <iomanip>

using namespace std;
using namespace std::chrono;
using namespace boost::program_options;
using namespace boost::filesystem;

#define ALPHA 2 

int getIntegerIndex(char symbol)
{
  switch(symbol) {
    case 'C':
      return 0;

    case 'E':
      return 1;

    case 'G':
      return 2;

    case 'H':
      return 3;

    case 'I':
      return 4;
  }
}

void updateSSTLengths(vector<vector<int>> &sst_lengths, string &file_name)
{
  ifstream file(file_name.c_str());
  string line;
  char symbol;
  int index;
  int length;
  while(getline(file,line)) {
    boost::char_separator<char> sep("\t ");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    int i = 0;
    BOOST_FOREACH(const string& t, tokens) {
      if (i == 1) {
        symbol = boost::lexical_cast<char>(t);
      } else if (i == 2) {
        length = boost::lexical_cast<int>(t);
      }
      i++;
    }
    //cout << symbol << "\t" << length << endl;
    int index = getIntegerIndex(symbol);
    sst_lengths[index].push_back(length);
  }
  file.close();
}

double getPoissonMean(vector<int> &lengths)
{
  double sum_lengths = 0;
  for (int i=0; i<lengths.size(); i++) {
    //cout << lengths[i] << " ";
    sum_lengths += lengths[i];
  }
  //cout << endl;
  double num = sum_lengths + 0.5;
  double denom = lengths.size() + 1.0/ALPHA;
  //cout << "sum_lengths: " << sum_lengths << "\tN: " << lengths.size() << endl;
  return num/denom;
}

int main(int argc, char **argv)
{
  vector<vector<int>> sst_lengths;
  vector<int> lengths;
  for (int i=0; i<5; i++) {
    sst_lengths.push_back(lengths);
  }
  string lengths_dir = "./dssp/sst_lengths/";
  path p(lengths_dir);
  cout << "path: " << p.string() << endl;
  if (exists(p)) { 
    if (is_directory(p)) { 
      vector<path> files; // store paths,
      copy(directory_iterator(p), directory_iterator(), back_inserter(files));
      cout << "# of profiles: " << files.size() << endl;
      for (int i=0; i<files.size(); i++) {
        string file_name = files[i].string();
        cout << i+1 << ". Reading " << file_name << endl;
        updateSSTLengths(sst_lengths,file_name);
      }
    }
  } 
  vector<double> means(sst_lengths.size(),0);
  for (int i=0; i<sst_lengths.size(); i++) {
    means[i] = getPoissonMean(sst_lengths[i]);
    cout << sst_lengths[i].size() << "\t" << means[i] << endl;
  }
}

