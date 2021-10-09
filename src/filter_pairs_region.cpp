/* filter_pairs_region.cpp
 * Seungsoo Kim
 * April 7, 2018
 * 
 * Takes two read-name-sorted BAM files and creates a pairs file for pairs where one read maps to a region e.g. rDNA.
 */

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

class Seqread {
  public:
    string id;
    int flag;
    string chr;
    long pos;
    int mapq;
    string seq;
    Seqread (string);
    long endpos (void);
    string strand (void);
};

Seqread::Seqread(string line) {
  stringstream ss(line);
  string tempid;
  int tempflag;
  string tempchr;
  long temppos;
  int tempmapq;
  string tempseq;
  string blank;

  if (!line.empty()) {
    ss >> tempid;
    ss >> tempflag;
    ss >> tempchr;
    ss >> temppos;
    ss >> tempmapq;
    ss >> blank;
    ss >> blank;
    ss >> blank;
    ss >> blank;
    ss >> tempseq;
    id = tempid;
    flag = tempflag;
    chr = tempchr;
    pos = temppos;
    mapq = tempmapq;
    seq = tempseq;
  }
}


long Seqread::endpos (void) {
  long p = pos;
  if (flag & 0x0010) {
    p = p + seq.length() - 1;
  }
  return p;
}

string Seqread::strand (void) {
  string out;
  if (flag & 0x0010) {
    out = "-";
  }
  else {
    out = "+";
  }
  return out;
}

long pairdist(Seqread r1, Seqread r2) {
  long dist = -1;
  if ((r1.chr.compare(r2.chr)==0) && (r1.flag != 4) && (r2.flag != 4)) {
    dist = r1.pos - r2.pos;
    if (dist < 0) dist = -dist;
  }
  return dist;
}

int main(int argc, char** argv) {
  // print pairs file header
  cout << "## pairs format v1.0\n";

  // process reads
  ifstream file1(argv[1]);
  ifstream file2(argv[2]);
  
  int minmapq;
  stringstream mapqss(argv[3]);
  mapqss >> minmapq;
  
  string selchr;
  stringstream selchrss(argv[4]);
  selchrss >> selchr;

  long selst;
  stringstream selstss(argv[5]);
  selstss >> selst;
  long selend;
  stringstream selendss(argv[6]);
  selendss >> selend;

  // skip SAM headers, if present
  string line1;
  string line2;
  getline(file1,line1);
  getline(file2,line2);
  while (line1.compare(0,1,"@") == 0) {
    getline(file1,line1);
  }
  while (line2.compare(0,1,"@") == 0) {
    getline(file2,line2);
  }
 
  while (!file1.eof() && !file2.eof()) {
    // read line from file 1
    Seqread read1(line1);
    getline(file1,line1);
    // read line from file 2
    Seqread read2(line2);    
    getline(file2,line2);

    // skip pair if both reads map to region
    if (((read1.chr.compare(selchr) == 0) && (read1.endpos() > selst) && (read1.endpos() < selend)) && ((read2.chr.compare(selchr) == 0) && (read2.endpos() > selst) && (read2.endpos() < selend))) {
      continue;
    }
    // skip pair if neither read maps to region
    if (!((read1.chr.compare(selchr) == 0) && (read1.endpos() > selst) && (read1.endpos() < selend)) && !((read2.chr.compare(selchr) == 0) && (read2.endpos() > selst) && (read2.endpos() < selend))) {
      continue;
    }

    // print pair data: 1) read ID, 2) chr1, 3) pos1, 4) chr2, 5) pos2, 6) strand1, 7) strand2, 8) insertsize, 9) orientation
    // print read mapping to region second
    if ((read1.chr.compare(selchr) == 0) && (read1.endpos() > selst) && (read1.endpos() < selend)) {
      cout << read2.chr << '\t';
      cout << read2.endpos() << '\t';
      cout << read1.id << '\t';
      cout << read1.chr << '\t';
      cout << read1.endpos() << '\t';
      cout << read2.strand() << '\t';
      cout << read1.strand() << '\t';
    }
    else {
      cout << read1.id << '\t';
      cout << read1.chr << '\t';
      cout << read1.endpos() << '\t';
      cout << read2.chr << '\t';
      cout << read2.endpos() << '\t';
      cout << read1.strand() << '\t';
      cout << read2.strand() << '\t';
    }
    if (read1.chr.compare(read2.chr) == 0) {
      if (read1.endpos() < read2.endpos()) {
        cout << read2.endpos() - read1.endpos() << '\t';
      }
      else {
        cout << read1.endpos() - read2.endpos() << '\t';
      }
    }
    else {
      cout << "-1\t";
    }
    if (read1.strand().compare(read2.strand()) == 0) {
      cout << "SAME\n";
    }
    else if (read1.endpos() < read2.endpos()) {
      if (read1.strand().compare(0,1,"+")==0) {
        cout << "IN\n";
      }
      else {
        cout << "OUT\n";
      }
    }
    else {
      if (read1.strand().compare(0,1,"+")==0) {
        cout << "OUT\n";
      }
      else {
        cout << "IN\n";
      }
    }
  } 
  
  file1.close();
  file2.close();
}
