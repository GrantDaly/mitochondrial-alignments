#pragma once
#include "htslib/hts.h"
#include <map>
#include <string>
#include <iostream>
class SampleBasePileup
{
 public:
  hts_pos_t pos = -1;
  std::string contig = "";
  char ref = 'N';
  // F1R2 Read 1's
  std::map<std::string, long int> F1R2R1s = {
    {"A", 0},
    {"T", 0},
    {"C", 0},
    {"G", 0}
  };

  // F1R2 Read 2's
    std::map<std::string, long int> F1R2R2s = {
    {"A", 0},
    {"T", 0},
    {"C", 0},
    {"G", 0}
  };
  
  // F2R1 Read 1's
    std::map<std::string, long int> F2R1R1s = {
    {"A", 0},
    {"T", 0},
    {"C", 0},
    {"G", 0}
  };

  // F2R1 Read 2's
    std::map<std::string, long int> F2R1R2s = {
    {"A", 0},
    {"T", 0},
    {"C", 0},
    {"G", 0}
  };
  //SampleBasePileup(hts_pos_t in_pos,  std::string in_contig);
  SampleBasePileup(hts_pos_t in_pos, std::string in_contig);
  SampleBasePileup(hts_pos_t in_pos, std::string in_contig, char in_ref);
  
  friend std::ostream& operator<<(std::ostream& os, const SampleBasePileup& pl);

  static void writeHeader(std::ostream & os){
    os << "Chrom" << "\t" << "Start" << "\t" << "ref" << "\t" << "base" << "\t" <<
      "F1R2R1" << "\t" << "F1R2R2" << "\t" << "F2R1R1" << "\t" << "F2R2R2" << std::endl;

  };

 
};
// SampleBasePileup(hts_pos_t in_pos,  std::string in_contig){
//     pos = in_pos;
//     contig = in_contig;
//   }
// SampleBasePileup::SampleBasePileup(hts_pos_t in_pos,  std::string in_contig): pos(in_pos), contig(in_contig) {}
SampleBasePileup::SampleBasePileup(hts_pos_t in_pos, std::string in_contig)
    : pos(in_pos), contig(in_contig) {}

SampleBasePileup::SampleBasePileup(hts_pos_t in_pos, std::string in_contig, char in_ref): pos(in_pos), contig(in_contig), ref(in_ref) {}
std::ostream& operator<<(std::ostream& os, const SampleBasePileup& pl)
  {
    os << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t" << "A" << "\t" << 
      pl.F1R2R1s.at("A") << "\t" << pl.F1R2R2s.at("A") << "\t" <<
      pl.F2R1R1s.at("A") << "\t" << pl.F2R1R2s.at("A") << std::endl  <<
      pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t" << "T" << "\t" <<
      pl.F1R2R1s.at("T") << "\t" << pl.F1R2R2s.at("T") << "\t" <<
      pl.F2R1R1s.at("T") << "\t" << pl.F2R1R2s.at("T") << std::endl <<
      pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t" << "C" << "\t" << 
      pl.F1R2R1s.at("C") << "\t" << pl.F1R2R2s.at("C") << "\t" <<
      pl.F2R1R1s.at("C") << "\t" << pl.F2R1R2s.at("C") << std::endl <<
      pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t" << "G" << "\t" <<
      pl.F1R2R1s.at("G") << "\t" << pl.F1R2R2s.at("G") << "\t" <<
      pl.F2R1R1s.at("G") << "\t" << pl.F2R1R2s.at("G") << std::endl;
                        
    return os;

  }
