#include "seqan3/argument_parser/argument_parser.hpp"
#include "seqan3/core/debug_stream/debug_stream_type.hpp"
#include <iostream>
#include <filesystem>

#include <iterator>

#include <hts-wrapper/faidx-cpp.h>
#include <hts-wrapper/bam-header.h>
#include <hts-wrapper/bam-record.h>
#include <hts-wrapper/cigar.h>
#include <hts-wrapper/sam-file.h>
#include <hts-wrapper/sam-iterator.h>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sequence_file/input.hpp>
//#include <htslib/faidx.h>


int main(int argc, char ** argv) {

  int test = fai_build("test.fasta");
  seqan3::argument_parser myparser{"NUMT_Identifier", argc, argv, seqan3::update_notifications::off};

  std::filesystem::path inBamName{};
  
  myparser.add_option(inBamName,'b',"bam","Input bam file sorted by position",
                        seqan3::option_spec::standard, seqan3::input_file_validator{{"bam"}});

  std::filesystem::path inFastaName{}; 
  myparser.add_option(inFastaName,'r',"ref","Reference Genome",
		      seqan3::option_spec::standard, seqan3::input_file_validator{{"fa","fasta"}});
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }


    auto inFasta = FastaWrapper{inFastaName};
    auto numt_bam = SamFile(inBamName, "rb");

    auto bam_header = numt_bam.read_header();

    if (numt_bam.index == nullptr) {
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
  auto mito_iter_wrapper = BamIterator(&numt_bam, &bam_header, main_chrom);
    
}
