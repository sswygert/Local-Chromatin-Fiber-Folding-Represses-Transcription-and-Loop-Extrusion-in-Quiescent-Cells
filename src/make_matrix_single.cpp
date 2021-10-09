/* make_matrix.cpp
 * takes assigned read pairs and creates matrix of interaction counts (long form)
 * Seungsoo Kim
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <map>

using namespace std;

int main(int argc, char** argv) {

  // arg 1: chr
  stringstream chrarg(argv[1]);
  string selectchr;
  chrarg >> selectchr;

  // arg 2: chrlen
  stringstream chrlenarg(argv[2]);
  long selectchrlen;
  chrlenarg >> selectchrlen;

  // arg 3: bin size
  stringstream bsizearg(argv[3]);
  long bsize;
  bsizearg >> bsize;

  // arg 4: pairs file
  ifstream countfile(argv[4]);

  // arg 5: minimum average value across row
  stringstream minrowavg_arg(argv[5]);
  double minrowavg;
  minrowavg_arg >> minrowavg;

  // arg 6: output row sum file
  ofstream outrowsumf(argv[6]);

  // arg 7: output matrix
//  ofstream outmatrf(argv[7]);

  // parse chromosome lengths  
  string line;
  
  // calculate total number of bins, and cumulative number of bins for each chromosome
  long nbin = selectchrlen/bsize+1;
  long minrowsum = (long) (minrowavg * nbin);

  // initialize matrices
  long* outmatr = new long[nbin*nbin];
  double* normmatr = new double[nbin*nbin];
  for (long i = 0; i < nbin*nbin; i++) {
    outmatr[i]=0;
    normmatr[i]=nan("");
  }
  
  string id;
  string chr1, chr2;
  long pos1, pos2;
  char dir1, dir2;
  int insert;
  string orient;
  long bin1, bin2;
  while (getline(countfile,line)) {
    if (line.compare(0,1,"#")==0) {
      continue;
    }

    // parse one read pair
    stringstream ss(line);
    ss >> id;
    ss >> chr1;
    ss >> pos1;
    ss >> chr2;
    ss >> pos2;
    ss >> dir1;
    ss >> dir2;
    ss >> insert;
    ss >> orient;
    // add to matrix count (symmetrically) if both reads are from chromosomes included in the matrix and not in same bin
    if ((chr1.compare(selectchr) == 0) && (chr2.compare(selectchr) == 0)) {
      bin1 = pos1/bsize;
      bin2 = pos2/bsize;
      if (bin1 != bin2) {
        outmatr[nbin*bin2 + bin1]++; 
        outmatr[nbin*bin1 + bin2]++; 
      }
    }  
  }

  // calculate row sums
  long* rowsums = new long[nbin];
  for (int i = 0; i < nbin; i++) {
    rowsums[i] = 0;
    for (int j = 0; j < nbin; j++) {
      rowsums[i] += outmatr[nbin*i + j];
    }
    outrowsumf << i << '\t' << rowsums[i] << '\n';
  }

  // calculate total sum, excluding rows that don't meet minimum
  long allsum = 0;
  for (int i = 0; i < nbin; i++) {
    if (rowsums[i] >= minrowsum) {
      allsum += rowsums[i];
    }
  }

  // calculate coverage-normalized matrix
  for (int i = 0; i < nbin; i++) {
    if (rowsums[i] >= minrowsum) {
      for (int j = 0; j < nbin; j++) {
        if (rowsums[j] >= minrowsum) {
          normmatr[nbin*i + j] = ((double)outmatr[nbin*i + j]*allsum)/((double)rowsums[i]*rowsums[j]);
        }
      }
    }
  }

  // print matrix, one line per entry of matrix
  for (int i = 0; i < nbin; i++) {
    for (int j = 0; j < nbin; j++) {
      cout << i << '\t' << j << '\t' << outmatr[nbin*i + j] << '\t' << normmatr[nbin*i + j] << '\n';
    }
  }

  delete[] outmatr;
  delete[] normmatr;
  delete[] rowsums;
  countfile.close();
  outrowsumf.close();
//  outmatrf.close();
}
