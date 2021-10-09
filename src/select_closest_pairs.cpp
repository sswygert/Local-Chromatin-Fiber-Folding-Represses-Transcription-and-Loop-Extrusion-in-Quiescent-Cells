/* select_closest_pair.cpp
 * Seungsoo Kim
 * January 7, 2018
 * 
 * Takes two read-name-sorted BAM files and selects the mapping for each read pair that minimizes the insert size. 
 */

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

string firstfield(string line) {
  stringstream ss(line);
  string field;
  ss >> field;
  return field;
};

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
  
  // skip SAM headers, if present
  string line1,line1next;
  string line2,line2next;
  getline(file1,line1next);
  getline(file2,line2next);
  while (line1.compare(0,1,"@") == 0) {
    getline(file1,line1next);
  }
  while (line2.compare(0,1,"@") == 0) {
    getline(file2,line2next);
  }
 
  string read1id, read2id;
  while (!file1.eof() && !file2.eof()) {
    vector<Seqread> f1mappings;
    vector<Seqread> f2mappings;
    // read lines from file 1 until read ID changes
    do {
      line1 = line1next;
      Seqread tempseqread(line1);
      read1id = tempseqread.id;
      f1mappings.push_back(tempseqread);
      getline(file1,line1next);
    } while (read1id.compare(firstfield(line1next)) == 0);

    // read lines from file 2 until read ID changes
    do {
      line2 = line2next;
      Seqread tempseqread(line2);
      read2id = tempseqread.id;
      f2mappings.push_back(tempseqread);
      getline(file2,line2next);
    } while (read2id.compare(firstfield(line2next)) == 0);
    
    // select pair
    int read1best = 0;
    int read2best = 0;
    int mindist = -1;
    for (int i = 0; i < f1mappings.size(); i++) {
      for (int j = 0; j < f2mappings.size(); j++) {
        long dist = pairdist(f1mappings[i],f2mappings[j]);
        if (dist > 0) {
          if (dist < mindist) {
            mindist = dist;
            read1best = i;
            read2best = j;
          }
        }
      }
    }    

    if ((f1mappings[read1best].flag == 4) || (f2mappings[read2best].flag == 4)) {
      continue;
    }
    
    // print pair data: 1) read ID, 2) chr1, 3) pos1, 4) chr2, 5) pos2, 6) strand1, 7) strand2, 8) insertsize, 9) orientation
    cout << f1mappings[read1best].id << '\t';
    cout << f1mappings[read1best].chr << '\t';
    cout << f1mappings[read1best].endpos() << '\t';
    cout << f2mappings[read2best].chr << '\t';
    cout << f2mappings[read2best].endpos() << '\t';
    cout << f1mappings[read1best].strand() << '\t';
    cout << f2mappings[read2best].strand() << '\t';
    if (f1mappings[read1best].chr.compare(f2mappings[read2best].chr) == 0) {
      if (f1mappings[read1best].endpos() < f2mappings[read2best].endpos()) {
        cout << f2mappings[read2best].endpos() - f1mappings[read1best].endpos() << '\t';
      }
      else {
        cout << f1mappings[read1best].endpos() - f2mappings[read2best].endpos() << '\t';
      }
    }
    else {
      cout << "-1\t";
    }
    if (f1mappings[read1best].strand().compare(f2mappings[read2best].strand()) == 0) {
      cout << "SAME\n";
    }
    else if (f1mappings[read1best].endpos() < f2mappings[read2best].endpos()) {
      if (f1mappings[read1best].strand().compare(0,1,"+")==0) {
        cout << "IN\n";
      }
      else {
        cout << "OUT\n";
      }
    }
    else {
      if (f1mappings[read1best].strand().compare(0,1,"+")==0) {
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
