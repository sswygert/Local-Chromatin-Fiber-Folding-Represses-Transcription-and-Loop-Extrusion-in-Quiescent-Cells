/* dist_decay.cpp
 * takes pairs file and creates histogram counts
 * Seungsoo Kim
 */
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char** argv) {

  // arg 1: pairs file
  ifstream pairsf(argv[1]);

  // arg 2: min bin size
  long minb;
  stringstream arg2(argv[2]);
  arg2 >> minb;
  
  // arg 3: max bin size
  long maxb;
  stringstream arg3(argv[3]);
  arg3 >> maxb;
  
  // arg 4: increment (multiplicative)
  double incr;
  stringstream arg4(argv[4]);
  arg4 >> incr;
  double logincr = log(incr);
  
  vector<double> binmins;
  int nbins = ceil((log(maxb)-log(minb))/logincr);
  vector<long> countin;
  vector<long> countout;
  vector<long> countsame;
  for (int i = 0; i < nbins; i++) {
    countin.push_back(0);
    countout.push_back(0);
    countsame.push_back(0);
    binmins.push_back(ceil((double)minb*pow(incr,i)));
  } 

  string line;
  string id;
  string chr1;
  long pos1;
  string chr2;
  long pos2;
  char strand1;
  char strand2;
  int insert;
  string orient;
  while (getline(pairsf, line)) {
    stringstream ss(line);
    ss >> id;
    ss >> chr1;
    ss >> pos1;
    ss >> chr2;
    ss >> pos2;
    ss >> strand1;
    ss >> strand2;
    ss >> insert;
    ss >> orient;

    if (insert < minb) continue;
    if (chr1.compare(chr2) != 0) continue;
    if (insert > maxb) insert = maxb;
    int binno = floor((log(insert)-log(minb))/logincr);
    if (orient.compare("IN") == 0) countin[binno]++;
    else if (orient.compare("OUT") == 0) countout[binno]++;
    else if (orient.compare("SAME") == 0) countsame[binno]++;
  }
  pairsf.close();

  cout << "binno" << '\t' << "min" << '\t' << "IN" << '\t' << "OUT" << '\t' << "SAME" << '\n';
  for (int i = 0; i < nbins; i++) {
    cout << i << '\t' << binmins[i] << '\t' << countin[i] << '\t' << countout[i] << '\t' << countsame[i] << '\n';
  }
}
