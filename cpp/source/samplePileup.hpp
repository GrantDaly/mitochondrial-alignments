#pragma once
#include "htslib/hts.h"
#include <map>
#include <string>
#include <iostream>
#include <memory>
#include "sampleBasePileup.hpp"
class SamplePileup
{
 public:
  hts_pos_t pos = -1;
  std::string contig = "";
  char ref = 'N';

  //std::unique_ptr<SampleBasePileup> operator[](std::string baseName);
  std::unique_ptr<SampleBasePileup>& getSampleBasePileup(std::string baseName);
  std::map<std::string,  std::unique_ptr<SampleBasePileup>> sampleBases;
  
  //SampleBasePileup(hts_pos_t in_pos,  std::string in_contig);
  SamplePileup(hts_pos_t in_pos, std::string in_contig);
  SamplePileup(hts_pos_t in_pos, std::string in_contig, char in_ref);
  
  friend std::ostream& operator<<(std::ostream& os, SamplePileup& pl);

  static void writeHeader(std::ostream & os){
    os << "Chrom" << "\t" << "Start" << "\t" << "ref" << "\t" << "base" << "\t" <<
      "F1R2R1" << "\t" << "F1R2R2" << "\t" << "F2R1R1" << "\t" << "F2R1R2" "\t" <<
      "Total" << "\t" << "VAF" << "\t"
       << "FR-Pass" << "\t" << "Artifact-Score" << std::endl;
  };
 
};

SamplePileup::SamplePileup(hts_pos_t in_pos, std::string in_contig, char in_ref)
    : pos(in_pos), contig(in_contig), ref(in_ref) {
  sampleBases["A"] = std::make_unique<SampleBasePileup>();
  sampleBases["T"] = std::make_unique<SampleBasePileup>();
  sampleBases["C"] = std::make_unique<SampleBasePileup>();
  sampleBases["G"] = std::make_unique<SampleBasePileup>();
  sampleBases["ins"] = std::make_unique<SampleBasePileup>();
  sampleBases["del"] = std::make_unique<SampleBasePileup>();  
}

std::unique_ptr<SampleBasePileup>&  SamplePileup::getSampleBasePileup(std::string baseName){
  return sampleBases[baseName];
}


std::ostream& operator<<(std::ostream& os, SamplePileup& pl)
  {
    long int totalBases = pl.getSampleBasePileup("A")->totalBases() +
       pl.getSampleBasePileup("T")->totalBases() +
       pl.getSampleBasePileup("C")->totalBases() +
       pl.getSampleBasePileup("G")->totalBases() +
       pl.getSampleBasePileup("ins")->totalBases() +
       pl.getSampleBasePileup("del")->totalBases();
    
    os << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "A" << "\t" << pl.getSampleBasePileup("A")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("A")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("A")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("A")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("A")->totalBases() << "\t"
       << pl.getSampleBasePileup("A")->totalBases() / (double) totalBases << "\t"
       << pl.getSampleBasePileup("A")->strandTest() << "\t"
       <<pl.getSampleBasePileup("A")->artifactScore() << std::endl

       << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "T" << "\t" << pl.getSampleBasePileup("T")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("T")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("T")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("T")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("T")->totalBases() << "\t"
       << pl.getSampleBasePileup("T")->totalBases() / (double) totalBases << "\t"
       << pl.getSampleBasePileup("T")->strandTest() << "\t"
       << pl.getSampleBasePileup("T")->artifactScore() << std::endl
      
       << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "C" << "\t" << pl.getSampleBasePileup("C")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("C")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("C")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("C")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("C")->totalBases() << "\t"
       << pl.getSampleBasePileup("C")->totalBases() / (double) totalBases << "\t"
       << pl.getSampleBasePileup("C")->strandTest() << "\t"
       << pl.getSampleBasePileup("C")->artifactScore() << std::endl

       << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "G" << "\t" << pl.getSampleBasePileup("G")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("G")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("G")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("G")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("G")->totalBases() << "\t"
       << pl.getSampleBasePileup("G")->totalBases() / (double) totalBases << "\t"
       << pl.getSampleBasePileup("G")->strandTest() << "\t"
       << pl.getSampleBasePileup("G")->artifactScore() << std::endl

       << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "ins" << "\t" << pl.getSampleBasePileup("ins")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("ins")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("ins")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("ins")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("ins")->totalBases() << "\t"
       << pl.getSampleBasePileup("ins")->totalBases() / (double) totalBases << "\t"
       << pl.getSampleBasePileup("ins")->strandTest() << "\t"
       << pl.getSampleBasePileup("ins")->artifactScore() << std::endl

       << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "del" << "\t" << pl.getSampleBasePileup("del")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("del")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("del")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("del")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("del")->totalBases() << "\t"
       << pl.getSampleBasePileup("del")->totalBases() / (double) totalBases << "\t"
       << pl.getSampleBasePileup("del")->strandTest() << "\t"
       << pl.getSampleBasePileup("del")->artifactScore() << std::endl;
                        
   return os;

 }
