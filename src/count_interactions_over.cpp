/* count_interactions_over.cpp
 * Seungsoo Kim
 * January 15, 2018
 * 
 */

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>

using namespace std;

int main(int argc, char** argv) {

  int MINLEN = 500;
  int MAXLEN = 10000;
 
  // read and process chromosome sizes 
  ifstream chrfile(argv[1]);
  string line;
  vector<string> chrs;
  map<string,long> chrlens;
  while (getline(chrfile, line)) {
    stringstream ss(line);
    string chr;
    long chrlen;
    ss >> chr;
    ss >> chrlen;
    chrlens[chr] = chrlen;
    chrs.push_back(chr);
  }
  chrfile.close();

  string currentchr;
  int* counts = NULL;
  
  // process pairs
  ifstream pairsfile(argv[2]);
  while (getline(pairsfile, line)) {
    if (line.compare(0,1,"#") == 0) {
      continue;
    }
    stringstream ss(line);
    string id;
    string chr1;
    long pos1;
    string chr2;
    long pos2;
    char strand1;
    char strand2;
    int insert;
    string orient;
    ss >> id;
    ss >> chr1;
    ss >> pos1;
    ss >> chr2;
    ss >> pos2;
    ss >> strand1;
    ss >> strand2;
    ss >> insert;
    ss >> orient;

    if (chr1 != chr2) continue;

    if (chr1 != currentchr) {
      for (int i = 0; i < chrlens[currentchr]; i++) {
        cout << currentchr << '\t' << i + 1 << '\t' << counts[i] << '\n';
      }
      delete [] counts;
      counts = NULL;

      counts = new int[chrlens[chr1]];
      for (int i = 0; i < chrlens[chr1]; i++) {
        counts[i] = 0;
      }
      currentchr = chr1;
    }

    if ((insert > MINLEN) && (insert < MAXLEN)) {
      for (int i = pos1; i < pos2 - 1; i++) {
        counts[i]++;
      }
    }

  }
  for (int i = 0; i < chrlens[currentchr]; i++) {
    cout << currentchr << '\t' << i + 1 << '\t' << counts[i] << '\n';
  }
  delete [] counts;

  pairsfile.close();
}
