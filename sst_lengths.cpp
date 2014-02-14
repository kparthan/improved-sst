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
  string dssp_file;
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
       ("dssp",value<string>(&parameters.dssp_file),"path to the file")
  ;
  variables_map vm;
  store(parse_command_line(argc,argv,desc),vm);
  notify(vm);

  return parameters;
}

vector<vector<string>> parseDSSP(string &dssp_file)
{
  ifstream file(dssp_file.c_str());
  int count = 0;
  string line;
  vector<string> chain_lines;
  vector<vector<string>> all_lines;

  int HASH = 3 - 1;
  char prev_chain_id,current_chain_id;
  bool begin = 0;
  bool ignore_current_chain = 1;
  bool new_chain = 1;
  int CHAIN_ID = 12 - 1;
  int CHAIN_BREAK = 14 - 1;
  int CHAIN_END = 15 - 1;

  vector<char> suitable_chains;
  int num_chains = 0;
  int num_suitable_chains = 0;
  ofstream log("dssp_parser.log",ios::app);
  log << dssp_file;
  while (getline(file,line)) {
    count++;
    if (line[HASH] == '#') {
      begin = 1;
      log << "\tParsing begins at line # " << count+1 << endl;
    }
    if (begin == 1 && line[HASH] != '#') {
      if (line[CHAIN_BREAK] == '!' && line[CHAIN_END] == '*') { // chain ends
        if (!ignore_current_chain) {
          all_lines.push_back(chain_lines);
          num_suitable_chains++;
          suitable_chains.push_back(current_chain_id);
        }
        new_chain = 1;
        chain_lines.clear();
      } else if (line[CHAIN_BREAK] == '!') {  // ignore the current chain
        ignore_current_chain = 1;
        log << "\tChain break in chain " << current_chain_id << " at line " 
            << count << endl;
      } else {
        current_chain_id = line[CHAIN_ID];
        if (!new_chain) {
          assert(current_chain_id == prev_chain_id);
        } else {
          new_chain = 0;
          ignore_current_chain = 0;
          num_chains++;
        }
        if (!ignore_current_chain) {
          chain_lines.push_back(line);
        }
        prev_chain_id = current_chain_id;
      }
    }
  }
  if (!ignore_current_chain) {
    all_lines.push_back(chain_lines);
    num_suitable_chains++;
    suitable_chains.push_back(current_chain_id);
  }
  log << "\tnum_chains: " << num_chains << ";"
      << "\tnum_suitable: " << num_suitable_chains << " | ";
  for (int i=0; i<all_lines.size(); i++) {
    log << "\t [" << suitable_chains[i] << "," << all_lines[i].size() << "]";
  }
  log << endl;
    
  file.close();

  ofstream check("check");
  for (int i=0; i<all_lines.size(); i++) {
    for (int j=0; j<all_lines[i].size(); j++) {
      check << all_lines[i][j] << endl;
    }
  }
  check.close();

  return all_lines;
}

string extractName(string &file)
{
  unsigned pos1 = file.find_last_of("/");
  unsigned pos2 = file.find(".");
  int length = pos2 - pos1 - 1;
  string sub = file.substr(pos1+1,length);
  return sub;
}

vector<char> getSequence(vector<string> &lines)
{
  int ASSIGN_ID = 17 - 1;
  vector<char> sequence;
  char symbol;
  for (int i=0; i<lines.size(); i++) {
    symbol = lines[i][ASSIGN_ID];
    if (!(symbol == 'E' || symbol == 'G' || symbol == 'H' || symbol == 'I')) {
      symbol = 'C';
    }
    sequence.push_back(symbol);
  }
  assert(lines.size() == sequence.size());
  /*for (int i=0; i<sequence.size(); i++) {
    cout << sequence[i];
  } cout << endl;*/
  return sequence;
}

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

void getLengths(vector<vector<string>> &all_lines, string output_file)
{
  ofstream out(output_file.c_str());
  int CHAIN_ID = 12 - 1;
  char prev,current;
  int iprev,icurrent;
  int c,e,g,h,i;
  c = e = g = h = i = 0;
  int frequency[5] = {0};
  for (int i=0; i<all_lines.size(); i++) {
    char chain = all_lines[i][0][CHAIN_ID];
    vector<char> sequence = getSequence(all_lines[i]);
    prev = sequence[0];
    iprev = getIntegerIndex(prev);
    frequency[iprev] = 1;
    for (int j=1; j<sequence.size(); j++) {
      current = sequence[j];
      if (current == prev) {
        frequency[iprev]++;
      } else if (current != prev) {
        out << chain << "\t\t" << prev << "\t\t" << frequency[iprev] << endl;
        frequency[iprev] = 0;
        prev = current;
        iprev = getIntegerIndex(prev);
        frequency[iprev] = 1;
      }
    }
    out << chain << "\t\t" << prev << "\t\t" << frequency[iprev] << endl;
  }
  out.close();
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);
  string name = extractName(parameters.dssp_file);
  cout << "name: " << name << endl;
  vector<vector<string>> all_lines = parseDSSP(parameters.dssp_file);
  string output = "./dssp/sst_lengths/" + name;
  cout << "output file: " << output << endl;
  getLengths(all_lines,output);
}

