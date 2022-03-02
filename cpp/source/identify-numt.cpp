#include "seqan3/argument_parser/argument_parser.hpp"
#include "seqan3/core/debug_stream/debug_stream_type.hpp"
#include <iostream>
#include <filesystem>

#include <iterator>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <htslib/faidx.h>

struct reference_storage_t
{
    std::vector<std::string> ids;
    std::vector<std::vector<seqan3::dna5>> seqs;
};

void read_reference(std::filesystem::path const & reference_path,
                    reference_storage_t & storage)
{
    seqan3::sequence_file_input reference_in{reference_path};
    for (auto && record : reference_in)
    {
        seqan3::debug_stream << "Seq: " << record.id() << "\n";
        storage.ids.push_back(record.id());
        storage.seqs.push_back(record.sequence());
    }
}

int main(int argc, char ** argv) {

  int test = fai_build("test.fasta");
  seqan3::argument_parser myparser{"NUMT_Identifier", argc, argv, seqan3::update_notifications::off};

  std::filesystem::path inBamName{};
  
  myparser.add_option(inBamName,'f',"file","The input file containing the sequences.",
                        seqan3::option_spec::standard, seqan3::input_file_validator{{"bam"}});

  //std::filesystem::path inFastaName{}; 
  //myparser.add_option(inFastaName,'r',"ref","Reference Genome",
  //		      seqan3::option_spec::standard, seqan3::input_file_validator{{"fa","fasta"}});
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    // takes too long to read in reference. Seeing if I can avoid this
    //reference_storage_t ref_storage{};
    //read_reference(inFastaName, ref_storage);
 
    seqan3::sam_file_input inBam{inBamName};
    //seqan3::sam_file_input inBam{inBamName, ref_storage.ids, ref_storage.seqs};
    //seqan3::sequence_file_input inFasta{inFastaName};
    auto sortType = inBam.header().sorting;
    seqan3::debug_stream << "Sorting type " << sortType << std::endl;

    auto mitoName = "chrM";
    
    if(sortType != "queryname"){
      std::cerr << "Bam File Not Sorted by Query Name" << std::endl;
    }

    auto it = inBam.begin();
    
    do{
    auto first = *it;
    it++;
    //auto nextIt =  std::ranges::next(it);
    auto second = *it;
    //seqan3::debug_stream << firstRecord << std::endl;
    // seqan3::debug_stream << "id:  " << first.id() << '\n';
    // seqan3::debug_stream << "read sequence: " << first.sequence() << '\n';
    // seqan3::debug_stream << "mapping position: " << first.reference_position() << '\n';
    // seqan3::debug_stream << "mapping quality: " << first.mapping_quality() << '\n';

    // could later only generate these for the first read. I'm having trouble finding access to the header without accessing individual records, as opposed to getting a pointer at initialization of the iterator
    auto tempHeader = first.header_ptr();
    auto mitoRef_ids = tempHeader->ref_ids();

    bool firstMito = false;
    std::string refNameFirst;
    if(auto refIDFirst = first.reference_id()){
      //	seqan3::debug_stream << "Ref ID " << refIDFirst << std::endl;
	//	seqan3::debug_stream << "Ref Name " << tempHeader->ref_ids()[22] << std::endl;
	//	seqan3::debug_stream << "Ref Name " << tempHeader->ref_ids()[*refIDFirst] << std::endl;
	refNameFirst = tempHeader->ref_ids()[*refIDFirst];
	if(refNameFirst == mitoName){
	  firstMito = true;
	    }
	//seqan3::debug_stream << "Ref Name " << refNameFirst << std::endl;
	}
    else{ continue;}
    bool secondMito = false;
    std::string refNameSecond;
    if(auto refIDSecond = second.reference_id()){
    // seqan3::debug_stream << std::endl;
    // seqan3::debug_stream << "random " <<  << std::endl ;
    //std::cout << " mito code " << mitoCode << std::endl;
       refNameSecond = tempHeader->ref_ids()[*refIDSecond];
      if(refNameSecond == mitoName){
      secondMito = true;
      }
    }
    
    else{continue;}

    // find possible numt
    if((firstMito && !secondMito) || (secondMito && !firstMito)){

        seqan3::debug_stream << "Possible NUMT " << std::endl;
	//	seqan3::debug_stream << "Read ID' " << first.id() << std::endl;
	seqan3::debug_stream << "Read One Chromosome  " << refNameFirst << std::endl;
	seqan3::debug_stream << "Read One Position  " << first.reference_position() << std::endl;
	seqan3::debug_stream << "Read Two Chromosome  " << refNameSecond << std::endl;
	seqan3::debug_stream << "Read Two Position  " << second.reference_position() << std::endl;
	//seqan3::debug_stream << "Read ID's Match " << first.sequence() << std::endl;
	
	
	//auto firstChr = tempHeaderFirst->ref_dict[first.reference_id()];
	// for(const auto [key, val] : tempHeader ->ref_dict){
	//   //std::cout << i << std::endl;
	//   seqan3::debug_stream <<"key: " <<  key << std::endl;
	//   seqan3::debug_stream << val << std::endl;
	// }
	
      }
    //first = std::ranges::next(it);
    }
    while(it != inBam.end());

    // for (auto & rec : inBam)
    // {
    //     seqan3::debug_stream << "id:  " << rec.id() << '\n';
    //     seqan3::debug_stream << "read sequence: " << rec.sequence() << '\n';
    //     seqan3::debug_stream << "mapping position: " << rec.reference_position() << '\n';
    //     seqan3::debug_stream << "mapping quality: " << rec.mapping_quality() << '\n';
 
    //     // there are more fields read on default
    // }
    return 0;

    
}
