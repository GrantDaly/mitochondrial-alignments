#include "seqan3/argument_parser/argument_parser.hpp"
#include "seqan3/core/debug_stream/debug_stream_type.hpp"
#include <cctype>
#include <fstream>
#include <iostream>
#include <filesystem>

#include <iterator>

#include <hts-wrapper/faidx-cpp.h>
#include <hts-wrapper/bam-header.h>
#include <hts-wrapper/bam-record.h>
#include <hts-wrapper/cigar.h>
#include <hts-wrapper/sam-file.h>
#include <hts-wrapper/sam-iterator.h>

#include <ostream>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <utility>
//#include <htslib/faidx.h>

// struct BaseCount {
//   long int A{0};
//   long int T{0};
//   long int C{0};
//   long int G{0};
// };

struct BaseCount {
  long int A{0}, T{0}, C{0}, G{0};
};

struct PileupCount {
  char ref;
  BaseCount bases;
};

// std::vector<base_count> makePileupCount(const std::string & inSequence)
std::vector<PileupCount>  makePileupCount(const std::string & inSequence)
{
  std::vector<PileupCount> out_vector(inSequence.length());

  for (long unsigned int i{0}; i < inSequence.length(); i++)
    {
      out_vector[i].ref = toupper(inSequence[i]);
    }
  return out_vector;
}

void updatePileups(const auto & read,  std::vector<PileupCount> & pileups)
{
  int current_query_offset{0};
  int next_query_offset{0};
  
  int current_subject_offset{0};
  int next_subject_offset{0};


  const auto start_pos = read.getPos();
  const auto ref_bases = read.getReferenceLength();
  const auto pos = read.getPos();
  
  const auto seq = read.getSequence();
  auto seq_view = std::string_view(seq.data(), seq.size());
  for(const auto cig: read.getCigar()){
    
    const auto query_consumed =  cig.queryBasesConsumed();

    auto subject_consumed{0};

    // for now don't handle soft clipping well under the hood so correcting
    if((cig.getCigar() == CigarOperationType::SoftClip) || (cig.getCigar() == CigarOperationType::HardClip))
      {
	subject_consumed = 0;
      }
    else {
    subject_consumed =  cig.subjectBasesConsumed();
    }
    

    auto match_count{0}, mismatch_count{0};
    next_subject_offset = current_subject_offset + subject_consumed;
    
    next_query_offset = current_query_offset + query_consumed;
    // check if cigar match and check mismatches
    // std::cout << " " << current_query_offset << "," <<  next_query_offset;
    if(cig.getCigar() == CigarOperationType::Match){
        
      auto query_match_string = seq_view.substr(current_query_offset, query_consumed);

      for(auto i{0}; i < query_match_string.length(); i++)
	{
	  auto & subject_pileup = pileups[pos + i + current_subject_offset];
	  auto ref_base = subject_pileup.ref;
	  auto query_base = toupper(query_match_string[i]);
	  switch(query_base){
	  case 'A':
	    subject_pileup.bases.A += 1;
	    break;
	  case 'T':
	    subject_pileup.bases.T += 1;
	    break;
	  case 'C':
	    subject_pileup.bases.C += 1;
	    break;
	  case 'G':
	    subject_pileup.bases.G += 1;
	    break;
	  default:
	    break;
	  }
	  // std::cout << ref_base << ">" << query_base;
	}
      // std::cout << "matches " << match_count << " mismatches " << mismatch_count << "|" << std::endl;
    }
    
  current_query_offset = next_query_offset;
  current_subject_offset = next_subject_offset;
  }
  
  // std::cout << std::endl;
}

void outputPileups(const std::vector<PileupCount> & mito, const std::vector<PileupCount> & numt, std::ostream & out) {

  if (mito.size() != numt.size()){
    std::cerr << "mito and numt pileup vectors not the same length" << std::endl;
    std::exit(1);
  }

  out << "Start" << "\t" << "ref" << "\t" << "base" << "\t"
       <<"mito_count" << "\t" << "numt_count" << "\t"
      << "numt_percent" << "\t" << "numt_vaf" <<  std::endl;
  
  for(auto i{0}; i < mito.size(); i++)
    {
      long int mito_depth = mito[i].bases.A + mito[i].bases.T + mito[i].bases.C + mito[i].bases.G;
      long int numt_depth = numt[i].bases.A + numt[i].bases.T + numt[i].bases.C + numt[i].bases.G;

      // of any given allele wat is the % numt vs total of just that allele
      double a_numt_percent{-1.0}, t_numt_percent{-1.0},
	c_numt_percent{-1.0}, g_numt_percent{-1.0};

      // of any given allele what is the % numt vs total of all bases, so comparable to VAF
      double a_numt_vaf_percent{-1.0}, t_numt_vaf_percent{-1.0},
	c_numt_vaf_percent{-1.0}, g_numt_vaf_percent{-1.0};

      if((mito[i].bases.A > 0) || (numt[i].bases.A > 0)) {
      a_numt_percent = ((double) (numt[i].bases.A) / (mito[i].bases.A + numt[i].bases.A)) *100;
      a_numt_vaf_percent = ((double) (numt[i].bases.A) / (mito_depth + numt_depth)) *100;
      }
      

      if((mito[i].bases.T > 0) || (numt[i].bases.T > 0))
	{
      t_numt_percent = ((double) (numt[i].bases.T) / (mito[i].bases.T + numt[i].bases.T)) *100;
      t_numt_vaf_percent = ((double) (numt[i].bases.T) / (mito_depth + numt_depth)) *100;
	}

      if((mito[i].bases.C > 0) || (numt[i].bases.C > 0))
	{
      c_numt_percent = ((double) (numt[i].bases.C) / (mito[i].bases.C + numt[i].bases.C)) *100;
      c_numt_vaf_percent = ((double) (numt[i].bases.C) / (mito_depth + numt_depth)) *100;
	}

      if((mito[i].bases.G > 0) || (numt[i].bases.G > 0))
	{
      g_numt_percent = ((double) (numt[i].bases.G) / (mito[i].bases.G + numt[i].bases.G)) *100;
      g_numt_vaf_percent = ((double) (numt[i].bases.G) / (mito_depth + numt_depth)) *100;
      
	}


      // A
      out << i+1 << "\t" << mito[i].ref << "\tA\t" << mito[i].bases.A
	  << "\t" << numt[i].bases.A << "\t" << a_numt_percent
	  << "\t" << a_numt_vaf_percent << std::endl;

      // T
      out << i+1 << "\t" << mito[i].ref << "\tT\t" << mito[i].bases.T
	  << "\t" << numt[i].bases.T << "\t" << t_numt_percent
	  << "\t" << t_numt_vaf_percent << std::endl;

     // C
      out << i+1 << "\t" << mito[i].ref << "\tC\t" << mito[i].bases.C
	  << "\t" << numt[i].bases.C << "\t" << c_numt_percent
       	<< "\t" << c_numt_vaf_percent << std::endl;

     // G
      out << i+1 << "\t" << mito[i].ref << "\tG\t" << mito[i].bases.G
	  << "\t" << numt[i].bases.G << "\t" << g_numt_percent
	<< "\t" << g_numt_vaf_percent << std::endl;
       }
    
}
int main(int argc, char ** argv) {

  
  seqan3::argument_parser myparser{"NUMT_Identifier", argc, argv, seqan3::update_notifications::off};

  std::filesystem::path inBamName{};
  std::string outPrefix{};
  
  myparser.add_option(inBamName,'b',"bam","Input bam file sorted by position",
                        seqan3::option_spec::standard, seqan3::input_file_validator{{"bam"}});

  std::filesystem::path inFastaName{}; 
  myparser.add_option(inFastaName,'r',"ref","Reference Genome",
		      seqan3::option_spec::standard, seqan3::input_file_validator{{"fa","fasta"}});

    myparser.add_option(outPrefix, 'p', "prefix", "Prefix");
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    if(outPrefix.length() > 0){
    std::cout << outPrefix << std::endl;
    }
    else{
      outPrefix = "out";
      std::cout << "Default output out.tsv" << std::endl;
    }
    auto inFasta = FastaWrapper{inFastaName};
    auto mito_bam = SamFile(inBamName, "rb");
    // numt bam is the same file which I make seperate iterators for the nuclear regions
    auto numt_bam = SamFile(inBamName, "rb");

    auto bam_header = mito_bam.read_header();

    if (mito_bam.index == nullptr) {
    std::cerr << "No Index for Bam File " << std::endl;
    std::exit(1);
  }

  std::string main_chrom{"chrM"};
  int main_chrom_tid = bam_header.contigNameToTid(main_chrom);


    std::map<std::string, std::pair<std::unique_ptr<BamRecord>,
				  std::unique_ptr<BamRecord>
			  >
	   >
      numt_read_pairs;

  auto counter{0};
  auto mito_iter_wrapper = BamIterator(&mito_bam, &bam_header, main_chrom);

  // vector to store var counts
  auto pileup_vec_mito = makePileupCount(inFasta.getRegion(main_chrom));
  auto pileup_vec_numt = makePileupCount(inFasta.getRegion(main_chrom));

  
  while (const auto & mito_record = mito_iter_wrapper.getNextRecord()) {
    auto mito_read = mito_record.value();
    auto const mito_read_name = mito_read.getQName();

    // needs to be primary alignment
    if(BAM_FSECONDARY == (mito_read.getFlag() & BAM_FSECONDARY))
      continue;

    auto const read_tid = mito_read.getTidNumber();
    auto const mate_tid = mito_read.getMateTidNumber();
    auto const read_contig_name = bam_header.contigNameToTid(read_tid);
    auto const mate_contig_name = bam_header.contigNameToTid(mate_tid);
    

    if (read_tid == mate_tid) {
      updatePileups(mito_read, pileup_vec_mito);
    }
    else if (read_tid != mate_tid){
      updatePileups(mito_read, pileup_vec_numt);
    }
    /* call */

    // not doing anything with the mate aligned to nuclear at this time, but later plan t
    // if (read_tid != mate_tid) {
  //           auto mate_iter_wrapper = BamIterator(&numt_bam, mate_tid,
  // 					   mate_pos, mate_pos+1);
  //     const auto & numt_record = mate_iter_wrapper.getRecordByName(mito_read_name);
  //     if(numt_record) {

  // 	std::cout << mito_read_name << std::endl;
  // 	auto numt_read = numt_record.value();
  //       auto const numt_read_flag = numt_read.getFlag();
  //     }
  //     else {
  // 	std::cerr << "Could not retreive mate" << std::endl;
  //     }
  // }
    counter++;
    if((counter % 10000000) == 0)
      {
	std::cerr << counter / 10000000 << "0M reads processed" << std::endl;
      }

  }

  std::ofstream outFile;

  outFile.open(outPrefix + ".tsv");

  outputPileups(pileup_vec_mito, pileup_vec_numt,outFile);
  outFile.close();
  return 0;
}
