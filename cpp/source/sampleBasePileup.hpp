#pragma once
#include "htslib/hts.h"
#include <map>
#include <string>
#include <iostream>
#include <math.h>

class SampleBasePileup
{
 public:
  
  // giving SampleBasePileup knowledge of which base this is
  // F1R2 Read 1's
  long int F1R2R1s = 0;
  long int F1R2R2s = 0;
  long int F2R1R1s = 0;
  long int F2R1R2s = 0;
  long int properly_paired = 0;
  long int improperly_paired = 0;
  SampleBasePileup & operator=(SampleBasePileup &&) = default;

  double artifactScore();
  double strandTest();
  long int totalBases();
  double fractionUnpaired();
  double getVAF(auto allSiteBases);
};

/*
  F1R2 R1 | F1R2 R2
  -----------------
  F2R1 R1 | F2R1 R2
  Strand artifacts are diagonals (F1R2 R1 and F2R1 R2 for +) and (F1R2 R2 and F2R1 R1 for +)
  Read orientation (pre-adapter) are (F1R2 R1 and R2 for +) and (F2R1 R1 and R2 for -)
 */
double SampleBasePileup::artifactScore() {
  // log2 ratio of F1R2 to F2R1 to detect orientation artifacts (pre-adapter)
  return log2( (double)(F1R2R1s + F1R2R2s + 1) / (F2R1R1s + F2R1R2s + 1) );
}
double SampleBasePileup::strandTest() {
  // flag for if mimimum number of reads go to both strands to avoid sequencing artifacts
  return log2( (double)(F1R2R1s + F2R1R2s + 1) / (F1R2R2s + F2R1R1s + 1));
}

long int SampleBasePileup::totalBases() {
  return F1R2R1s + F1R2R2s + F2R1R1s + F2R1R2s;
}

double SampleBasePileup::fractionUnpaired() {
  if(improperly_paired == 0)
    return 0;
  if(properly_paired == 0)
    return 1;
  
  return (double)(improperly_paired)  / ( properly_paired + improperly_paired);
}

double SampleBasePileup::getVAF(auto allSiteBases) {
  if(allSiteBases > 0) {
    return (double) totalBases() / allSiteBases;
  }
  else {
    return -1.0;
  }
}
