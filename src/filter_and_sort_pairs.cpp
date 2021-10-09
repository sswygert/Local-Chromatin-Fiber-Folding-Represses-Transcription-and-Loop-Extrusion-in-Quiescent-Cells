/* filter_pairs.cpp
 * Seungsoo Kim
 * January 11, 2018
 * 
 * Takes two read-name-sorted BAM files and a MAPQ filter (>=) creates pairs file with pairs where both reads pass MAPQ filter. 
 */

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <tuple>
#include <string>
#include <numeric>
#include <stdlib.h>
#include <type_traits>
#include <stdexcept>

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
    char strand (void);
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

char Seqread::strand (void) {
  char out;
  if (flag & 0x0010) {
    out = '-';
  }
  else {
    out = '+';
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

//c++14 can compare tuple field-by-field 
//so that the built record is sorted accordingly
//chr1 chr2 pos1 pos2 strand1 strand2 insertSize orientation
using Record = tuple<string,string,long,long,char,char,int,string>; 

void outputRecord(const string& fname, const map<Record, string>& records, const string& header) {
  fstream fsout(fname, fstream::out);
  fsout << header << "\n";
  //insertSize and orientation are added by Seungsoo
  fsout << "#readID chr1 position1 chr2 position2 strand1 strand2 insertSize orientation\n";
  for(const auto& record : records) {
    const auto& other = record.first;
    fsout << record.second; //qname
    fsout << "\t";
    fsout << get<0>(other); //chr1
    fsout << "\t";
    fsout << get<2>(other); //pos1
    fsout << "\t";
    fsout << get<1>(other); //chr2
    fsout << "\t";
    fsout << get<3>(other); //pos2
    fsout << "\t";
    fsout << get<4>(other); //strand1
    fsout << "\t";
    fsout << get<5>(other); //strand2
    fsout << "\t";
    fsout << get<6>(other); //insertSize
    fsout << "\t";
    fsout << get<7>(other); //orientation
    fsout << "\n";
  }
  fsout.close();
} 

template < class T >
double percent(const T& n, const T& N) {
  return static_cast<double>(n) / N * double(100);
}

int main(int argc, char** argv) {
  const string args("#filter_and_sort_pairs samfile1 samefile2 fchrSize MinMapQ minInsertSize bSepOrient Outprefix\n"); 
  if(argc != 8) { cerr << args; exit(-1); }

  //output a header as required by 4DN-DCIC
  const string header = "## pairs format v1.0\n" +
    args + accumulate(argv, argv+argc, string("#"),
    [&](const string& a, const string& b) {
      return string(a)+" "+string(b);
    }) + string("\n");

  int k = 1;
  ifstream file1(argv[k++]);
  ifstream file2(argv[k++]);
  const std::string fChrSize(argv[k++]);
  ifstream fsChrSize(fChrSize);
  int minmapq, minInsertSize;
  stringstream mapqss(argv[k++]);
  mapqss >> minmapq;
  stringstream minInsertSizeSS(argv[k++]);
  minInsertSizeSS >> minInsertSize;
  const bool bSepOrient = stoi(argv[k++]);
  const string outputPrefix(argv[k++]);
  
  //read and build a map chrname -> index of chromosome
  //This contains only the genome chromosomes
  if(!fsChrSize.is_open()) {
    throw std::runtime_error("Can't open "+fChrSize+" for read");
  }
  string line;
  vector<string> chrs;
  map<string,int> chrnum;
  int chrn = 0;
  while (getline(fsChrSize, line)) {
    stringstream ss(line);
    string chr;
    ss >> chr;
    chrs.push_back(chr);
    chrnum[chr] = chrn;
    chrn++;
  }
  fsChrSize.close();
  if(chrnum.empty()) {
    throw std::runtime_error("Chromosome name -> size map is empty");
  }
  
  //map: Record -> qname
  //always keep the records separately and only output them
  //differently depending on bSepOrient
  map<Record,string> recordsIn, recordsOut, recordsSame; 

  // process reads
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
 
  std::size_t nFMapQ = 0, nFNoChr = 0, nFin = 0, nFDup = 0, nRemain = 0;
  std::size_t nIN = 0, nOUT = 0, nSAME = 0;
  std::size_t nTotal = 0;
  //parse the sam files and filter and build pair records
  while (!file1.eof() && !file2.eof()) {
    ++nTotal;
    // read line from file 1
    Seqread read1(line1);
    getline(file1,line1);
    // read line from file 2
    Seqread read2(line2);    
    getline(file2,line2);

    if ((read1.mapq < minmapq) || (read2.mapq < minmapq)) {
      ++nFMapQ;
      continue;
    }

    //exclude chromosomes not listed, e.g., mitochrondrial DNA
    if(chrnum.find(read1.chr) == chrnum.end() || 
       chrnum.find(read2.chr) == chrnum.end()) {
      ++nFNoChr;
      continue;
    }

    int insertSize = -1;
    if (read1.chr.compare(read2.chr) == 0) {
      if (read1.endpos() < read2.endpos()) {
        insertSize = read2.endpos() - read1.endpos();
      }
      else {
        insertSize = read1.endpos() - read2.endpos();
      }
    }

    string orientation;
    if (read1.strand() == read2.strand()) {
      orientation = "SAME";
    }
    else if (read1.endpos() < read2.endpos()) {
      if (read1.strand() == '+') {
        orientation = "IN";
      }
      else {
        orientation = "OUT";
      }
    }
    else {
      if (read1.strand() == '+') {
        orientation = "OUT";
      }
      else {
        orientation = "IN";
      }
    }

    if(insertSize != -1 && orientation == "IN" && insertSize < minInsertSize) {
      ++nFin;
      continue;
    }

    //canonicalize the read order so that smaller chromosome and position 
    //is always at front
    const bool bSwap = (chrnum[read1.chr] > chrnum[read2.chr]) ||
      ((chrnum[read1.chr]==chrnum[read2.chr]) && (read1.endpos() > read2.endpos()));
    Record record = bSwap ? 
      make_tuple(read2.chr,read1.chr,read2.endpos(),read1.endpos(),read2.strand(),read1.strand(),insertSize,orientation) : 
      make_tuple(read1.chr,read2.chr,read1.endpos(),read2.endpos(),read1.strand(),read2.strand(),insertSize,orientation);

    auto& records = 
      orientation == "IN" ? recordsIn : 
      orientation == "OUT" ? recordsOut : 
      recordsSame;
    //exclude duplicate reads
    if(records.find(record) != records.end()) {
      ++nFDup;
      continue;
    }

    ++nRemain;
    if(orientation == "IN") { ++nIN; }
    else if(orientation == "OUT") { ++nOUT; }
    else { ++nSAME; }

    records[record] = read1.id;
  }
  
  file1.close();
  file2.close();

  if(bSepOrient) {
    outputRecord(outputPrefix + ".IN.pair", recordsIn, header);
    outputRecord(outputPrefix + ".OUT.pair", recordsOut, header);
    outputRecord(outputPrefix + ".SAME.pair", recordsSame, header);
  } else {
    //Merge the records
    auto records = recordsIn;
    for(const auto& recordOut : recordsOut) {
      records.emplace(recordOut.first, recordOut.second);
    }
    for(const auto& recordSame : recordsSame) {
      records.emplace(recordSame.first, recordSame.second);
    }
    outputRecord(outputPrefix + ".pair", records, header);
  }

  //Output the total number of reads processed and filtered
  const std::size_t nF = nFMapQ + nFNoChr + nFDup + nFin;
  const std::size_t nFpnR = nF + nRemain;
  fstream fsout(outputPrefix + ".filtersort.log", fstream::out);
  fsout << header << "\n";
  fsout << "# Total of number of input read pairs: " << nTotal << " ( "  << percent(nTotal, nTotal) << "% )\n";
  fsout << "# number of filtered read pairs with mapQ < " << minmapq << ": " << nFMapQ << " ( "  << percent(nFMapQ, nTotal) << "% )\n";
  fsout << "# number of filtered read pairs from unknown chromosome: " << nFNoChr << " ( "  << percent(nFNoChr, nTotal) << "% )\n";
  fsout << "# number of filtered read pairs cis IN-IN within " <<  minInsertSize << ": " << nFin << " ( "  << percent(nFin, nTotal) << "% )\n";
  fsout << "# number of filtered read pairs duplication: " << nFDup << " ( "  << percent(nFDup, nTotal) << "% )\n";
  fsout << "# number of read pairs cis+trans IN-IN output: " << nIN  << " ( "  << percent(nIN, nTotal) << "% )\n";
  fsout << "# number of read pairs cis+trans OUT-OUT output: " << nOUT << " ( "  << percent(nOUT, nTotal) << "% )\n"; 
  fsout << "# number of read pairs cis+trans SAME-SAME output: " << nSAME << " ( "  << percent(nSAME, nTotal) << "% )\n";
  fsout << "# Total number of read pairs output: " << nRemain << " ( "  << percent(nRemain, nTotal) << "% )\n";
  fsout << "# Total number of filtered + output read pairs: " << nFpnR << " ( "  << percent(nFpnR, nTotal) << "% )\n";
  fsout.close();
}
