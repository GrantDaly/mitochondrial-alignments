#pragma once
#include "htslib/hts.h"
#include <map>
#include <string>
#include <iostream>
#include <memory>
#include <algorithm>
#include <random>

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
    os << "Chrom" << "\t" << "Start" << "\t" << "ref" << "\t" << "alt" << "\t" <<
      "F1R2R1" << "\t" << "F1R2R2" << "\t" << "F2R1R1" << "\t" << "F2R1R2" "\t" <<
      "Total-Alt-HQ" << "\t"<< "Site-Depth-HQ" << "\t" << "VAF" << "\t"
       << "FR-Pass" << "\t" << "Artifact-Score" << "\t" <<  "Fraction-Unpaired"  << std::endl;
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

// taking total site depth and assumption that worst case error is the alignment
// error (30Q -> 1/1000), what's the probability that number of bases observed
// by chance?
struct doubleBaseHolder_t {
  // double a{0}, t{0}, c{0}, g{0}, ins{0}, del{0};
  double a=0, t=0, c=0, g=0, ins=0, del=0;
};

struct longBaseHolder_t {
  // double a{0}, t{0}, c{0}, g{0}, ins{0}, del{0};
  long int a=0, t=0, c=0, g=0, ins=0, del=0;
};

  
doubleBaseHolder_t depthTest(long int totalBases, longBaseHolder_t alts, doubleBaseHolder_t unpairedFractions, const int numSamples = 100) {
  // int a[numSamples], t[numSamples], c[numSamples], g[numSamples], ins[numSamples], del[numSamples];
  if( totalBases < 500 ){
    doubleBaseHolder_t pValues = {1,1, 1, 1, 1, 1};
    return pValues;
  }
    std::vector<int> a_diffs(numSamples);
    std::vector<int> t_diffs(numSamples);
    std::vector<int> c_diffs(numSamples);
    std::vector<int> g_diffs(numSamples);
    std::vector<int> ins_diffs(numSamples);
    std::vector<int> del_diffs(numSamples);

    // effective error rate is maximum of base/mapping quality and % unpaired
    doubleBaseHolder_t errors{
      std::max(0.01, unpairedFractions.a),
      std::max(0.01, unpairedFractions.t),
      std::max(0.01, unpairedFractions.c),
      std::max(0.01, unpairedFractions.g),
      std::max(0.01, unpairedFractions.ins),
      std::max(0.01, unpairedFractions.del)
    };
    //std::cout << "A error % " << unpairedFractions.a << std::endl;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::binomial_distribution<typeof(totalBases)> a_errors(totalBases, errors.a);
  std::binomial_distribution<typeof(totalBases)> t_errors(totalBases, errors.t);
  std::binomial_distribution<typeof(totalBases)> c_errors(totalBases, errors.c);
  std::binomial_distribution<typeof(totalBases)> g_errors(totalBases, errors.g);
  std::binomial_distribution<typeof(totalBases)> ins_errors(totalBases, errors.ins);
  std::binomial_distribution<typeof(totalBases)> del_errors(totalBases, errors.del);

  int  a_draw;
  int  t_draw;
  int  c_draw;
  int  g_draw;
  int  ins_draw;
  int  del_draw;
  // std::cout << "Error Rates A " << errors.a << " T " << errors.t <<
  //   " C " << errors.c << " G "  << errors.g << std::endl;
  for (int i=0; i < numSamples; i++){
    a_draw = a_errors(gen);
    t_draw = t_errors(gen);
    c_draw = c_errors(gen);
    g_draw = g_errors(gen);
    ins_draw = ins_errors(gen);
    del_draw = del_errors(gen);

    // null hypothesis alt !>= error rate. 
    a_diffs[i] = alts.a - a_draw;
    t_diffs[i] = alts.t - t_draw;
    c_diffs[i] = alts.c - c_draw;
    g_diffs[i] = alts.g - g_draw;
    ins_diffs[i] = alts.ins - ins_draw;
    del_diffs[i] = alts.del - del_draw;
    // if(a_draw > 0){
    // std::cout << "Depth : " << totalBases
    // << " A errors : " << a_draw << "A observed " << alts.a <<
    //   " diff obs - draw " << alts.a - a_draw << std::endl;

    // }
  }
    auto gt_zero_a = std::count_if(a_diffs.begin(), a_diffs.end(),
				   [](long int i) { return i < 0; });
    auto pval_a = (double) gt_zero_a / numSamples;
    std::cout << "Bases " << totalBases
	      << " A " << alts.a << " Error " << errors.a << " Number > 0 " << gt_zero_a << " p-value " << pval_a << std::endl;
    
    auto gt_zero_t = std::count_if(t_diffs.begin(), t_diffs.end(),
				   [](long int i) { return i <= 0; });
    auto pval_t = (double) gt_zero_t / numSamples;
        std::cout << "Bases " << totalBases
		  << " T " << alts.t << " Error " << errors.t
		  << " Number > 0 " << gt_zero_t << " p-value " << pval_t << std::endl;
    auto gt_zero_c = std::count_if(c_diffs.begin(), c_diffs.end(),
				   [](long int i) { return i <= 0; });
    auto pval_c = (double) gt_zero_c / numSamples;
        std::cout << "Bases " << totalBases
		  << " C " << alts.c << " Error " << errors.c
		  << " Number > 0 " << gt_zero_c << " p-value " << pval_c << std::endl;
    auto gt_zero_g = std::count_if(g_diffs.begin(), g_diffs.end(),
				   [](long int i) { return i <= 0; });
    auto pval_g = (double) gt_zero_g / numSamples;
        std::cout << "Bases " << totalBases
		  << " G " << alts.g << " Error " << errors.g
	      << " Number > 0 " << gt_zero_g << " p-value " << pval_g << std::endl;

    auto gt_zero_ins = std::count_if(ins_diffs.begin(), ins_diffs.end(),
				   [](long int i) { return i <= 0; });
    auto pval_ins = (double) gt_zero_ins / numSamples;
    auto gt_zero_del = std::count_if(del_diffs.begin(), del_diffs.end(),
				   [](long int i) { return i <= 0; });
    auto pval_del = (double) gt_zero_del / numSamples;


    doubleBaseHolder_t pValues = {pval_a, pval_t, pval_c, pval_g, pval_ins, pval_del};
  return pValues;


}
std::ostream& operator<<(std::ostream& os, SamplePileup& pl)
  {
    long int totalBases = pl.getSampleBasePileup("A")->totalBases() +
       pl.getSampleBasePileup("T")->totalBases() +
       pl.getSampleBasePileup("C")->totalBases() +
       pl.getSampleBasePileup("G")->totalBases() +
       pl.getSampleBasePileup("ins")->totalBases() +
       pl.getSampleBasePileup("del")->totalBases();

    
    doubleBaseHolder_t unpairedFractions = {
      pl.getSampleBasePileup("A")->fractionUnpaired(),
      pl.getSampleBasePileup("T")->fractionUnpaired(),
      pl.getSampleBasePileup("C")->fractionUnpaired(),
      pl.getSampleBasePileup("G")->fractionUnpaired(),
      pl.getSampleBasePileup("ins")->fractionUnpaired(),
      pl.getSampleBasePileup("del")->fractionUnpaired()
    };
    
    longBaseHolder_t totalAlt = {
      pl.getSampleBasePileup("A")->totalBases(),
      pl.getSampleBasePileup("T")->totalBases(),
      pl.getSampleBasePileup("C")->totalBases(),
      pl.getSampleBasePileup("G")->totalBases(),
      pl.getSampleBasePileup("ins")->totalBases(),
      pl.getSampleBasePileup("del")->totalBases()
    };
    
      doubleBaseHolder_t basePvals = depthTest(totalBases,totalAlt, unpairedFractions, 10000);


      if((pl.getSampleBasePileup("A")->totalBases() > 0) && (pl.ref != 'A')) {
    os << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "A" << "\t" << pl.getSampleBasePileup("A")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("A")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("A")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("A")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("A")->totalBases() << "\t"
       << totalBases << "\t"
       << pl.getSampleBasePileup("A")->getVAF(totalBases) << "\t"
       << pl.getSampleBasePileup("A")->strandTest() << "\t"
       <<pl.getSampleBasePileup("A")->artifactScore() << "\t"
       <<pl.getSampleBasePileup("A")->fractionUnpaired() << "\t"
       << basePvals.a << std::endl;
    }
    
    if((pl.getSampleBasePileup("T")->totalBases() > 0) && (pl.ref != 'T')) {
       os << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "T" << "\t" << pl.getSampleBasePileup("T")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("T")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("T")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("T")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("T")->totalBases() << "\t"
       << totalBases << "\t"
       << pl.getSampleBasePileup("T")->getVAF(totalBases) << "\t"
       << pl.getSampleBasePileup("T")->strandTest() << "\t"
       << pl.getSampleBasePileup("T")->artifactScore() << "\t"
	  <<pl.getSampleBasePileup("T")->fractionUnpaired() << "\t"
	 << basePvals.t << std::endl;
    }
          if((pl.getSampleBasePileup("C")->totalBases() > 0) && (pl.ref != 'C')) {
       os << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "C" << "\t" << pl.getSampleBasePileup("C")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("C")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("C")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("C")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("C")->totalBases() << "\t"
       << totalBases << "\t"
       << pl.getSampleBasePileup("C")->getVAF(totalBases) << "\t"
       << pl.getSampleBasePileup("C")->strandTest() << "\t"
       << pl.getSampleBasePileup("C")->artifactScore() << "\t"
	  <<pl.getSampleBasePileup("C")->fractionUnpaired() << "\t"
		 << basePvals.c << std::endl;
          }
      if((pl.getSampleBasePileup("G")->totalBases() > 0) && (pl.ref != 'G')) {
       os << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "G" << "\t" << pl.getSampleBasePileup("G")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("G")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("G")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("G")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("G")->totalBases() << "\t"
       << totalBases << "\t"
       << pl.getSampleBasePileup("G")->getVAF(totalBases) << "\t"
       << pl.getSampleBasePileup("G")->strandTest() << "\t"
       << pl.getSampleBasePileup("G")->artifactScore() << "\t"
	  <<pl.getSampleBasePileup("G")->fractionUnpaired() << "\t"
	 << basePvals.g << std::endl;
    }
    if((pl.getSampleBasePileup("ins")->totalBases() > 0)) {
       os << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "ins" << "\t" << pl.getSampleBasePileup("ins")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("ins")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("ins")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("ins")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("ins")->totalBases() << "\t"
       << pl.getSampleBasePileup("ins")->getVAF(totalBases) << "\t"
       << pl.getSampleBasePileup("ins")->strandTest() << "\t"
       << pl.getSampleBasePileup("ins")->artifactScore() << "\t"
	  <<pl.getSampleBasePileup("ins")->fractionUnpaired() << "\t"
	 << basePvals.ins << std::endl;
    }
    if((pl.getSampleBasePileup("del")->totalBases() > 0)) {
       os << pl.contig << "\t" << pl.pos << "\t" << pl.ref << "\t"
       << "del" << "\t" << pl.getSampleBasePileup("del")->F1R2R1s
       << "\t" << pl.getSampleBasePileup("del")->F1R2R2s
       << "\t" << pl.getSampleBasePileup("del")->F2R1R1s << "\t"
       << pl.getSampleBasePileup("del")->F2R1R2s << "\t"
       << pl.getSampleBasePileup("del")->totalBases() << "\t"
       << totalBases << "\t"
       << pl.getSampleBasePileup("del")->getVAF(totalBases) << "\t"
       << pl.getSampleBasePileup("del")->strandTest() << "\t"
       << pl.getSampleBasePileup("del")->artifactScore() << "\t"
       <<pl.getSampleBasePileup("del")->fractionUnpaired() << "\t"
	  << basePvals.del << std::endl;
	}
                        
   return os;

 }
