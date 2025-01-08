#include "htslib/hts.h"
#include "seqan3/argument_parser/argument_parser.hpp"
#include "seqan3/core/debug_stream/debug_stream_type.hpp"
#include <iostream>
#include <filesystem>
#include <random>
#include <iterator>
#include <ranges>
#include <fstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sam_file/output.hpp>

#include <seqan3/alphabet/views/complement.hpp> 

//#include <htslib/faidx.h>

#include "hts-wrapper/faidx-cpp.h"
#include "discrete-distributions.hpp"

using namespace seqan3::literals;

// struct reference_storage_t
// {
//     std::vector<std::string> ids;
//     std::vector<std::vector<seqan3::dna5>> seqs;
// };

// void read_reference(std::filesystem::path const & reference_path,
//                     reference_storage_t & storage)
// {
//     seqan3::sequence_file_input reference_in{reference_path};
//     for (auto && record : reference_in)
//     {
//         seqan3::debug_stream << "Seq: " << record.id() << "\n";
//         storage.ids.push_back(record.id());
//         storage.seqs.push_back(record.sequence());
//     }
// }
struct Bed{
  std::string chrom;
  hts_pos_t start;
  hts_pos_t end;
  std::string name;
  double score;
  std::string strand;

  friend std::istream& operator>>(std::istream& in, Bed& bed);
  friend std::ostream& operator<<(std::ostream& out, Bed& bed);
};

std::istream& operator>>(std::istream& in, Bed& bed){
  in >> bed.chrom >> bed.start >> bed.end >> bed.name >> bed.score >> bed.strand;
  return in;
} 
std::ostream& operator<<(std::ostream& out, Bed& bed){
  out << bed.chrom << "\t" << bed.start << "\t"
  << bed.end << "\t" << bed.name << "\t"
  << bed.score << "\t" << bed.strand;
  return out;
}

seqan3::dna5_vector stringToDNA(std::string inString) {
  //took literal overload from dna5 source code
     seqan3::dna5_vector outVector;
     auto n = inString.length();
     outVector.resize(n);

     
     for (size_t i = 0; i < n; ++i)
     {
     //    outVector[i].assign_char(s[i]);
       outVector[i].assign_char(inString.at(i));
      
     }
       
       //seqan3::debug_stream << outVector << std::endl;
   return outVector;
}

// std::pair<seqan3::dna5_vector,seqan3::dna5_vector>
// createPairedEndReads(const seqan3::dna5_vector & insertVector, int readLength = 150) {
//   int insert_size = insertVector.size();

//   auto insertRCVector = insertVector | seqan3::views::complement | std::views::reverse;			 
//     // seqan3::views::to<seqan3::dna5_vector>;
//   auto correctedReadLength = std::min(readLength, insert_size);
//   seqan3::dna5_vector read_one(insertVector.begin(), insertVector.begin() + correctedReadLength);
//   seqan3::dna5_vector read_two(insertRCVector.begin(), insertRCVector.begin() + correctedReadLength);

//   std::pair<seqan3::dna5_vector,seqan3::dna5_vector> return_pair{read_one, read_two};
//   return  return_pair;

// 				  }
//std::pair<T, T>
template <std::ranges::bidirectional_range T>
std::pair<T,T>
createPairedEndReads(const T &  insertVector, int readLength = 150) {
  int insert_size = insertVector.size();

  auto insertRCVector = insertVector | seqan3::views::complement | std::views::reverse
    | std::ranges::to<std::vector>();			 
    // seqan3::views::to<seqan3::dna5_vector>;
  auto correctedReadLength = std::min(readLength, insert_size);
  T read_one(insertVector.begin(), insertVector.begin() + correctedReadLength);
  T read_two(insertRCVector.begin(), insertRCVector.begin() + correctedReadLength);

  std::pair<T,T> return_pair{read_one, read_two};
  return  return_pair;

				  }
// hts_pos_t nuclearFragmentDistribution(void) { 
// hts_pos_t nuclearFragmentDistribution(void) {

//   // hard coding for now. eventually will parameterize 
// 			  return 165;
// }

int main(int argc, char ** argv) {

  // int test = fai_build("test.fasta");
  seqan3::argument_parser arg_parser{"NUMT_Simulator", argc, argv, seqan3::update_notifications::off};

  // std::filesystem::path inBamName{};
  
  // myparser.add_option(inBamName,'f',"file","The input file containing the sequences.",
  //                       seqan3::option_spec::standard, seqan3::input_file_validator{{"bam"}});

  std::filesystem::path inFastaName{}; 
  arg_parser.add_option(inFastaName,'r',"ref","Reference Genome",
		      seqan3::option_spec::standard, seqan3::input_file_validator{{"fa","fasta"}});

  std::filesystem::path bedName{};
  arg_parser.add_option(bedName, 'b', "bed", "Bed file name (bed must be 6 column format",
			seqan3::option_spec::required, seqan3::input_file_validator{{"bed"}});

  std::filesystem::path outPrefix{};
  arg_parser.add_option(outPrefix, 'o', "out-prefix", "Output Fastq Prefix");

  // add read length, target coverage, overhang
  int overhang_length = 300;
  int read_length = 150;
  int target_coverage = 10000;
  std::string distString{"uniform"};

  // appears to be a bug where 0 overhang to very short interval throws std::out_of_range. Shouldn't be encountered for meaningful intervals b/c either they're long or default 300bp overhang ocurrs
  // later can figure out where this is happening.
  arg_parser.add_option(overhang_length, 'x', "overhang", "Number of Bases Overhang");
  arg_parser.add_option(read_length, 'l', "read-length", "PE read length");
  arg_parser.add_option(target_coverage, 'c', "coverage", "Target Coverage");
  arg_parser.add_option(distString, 'd', "dist", "Distribution (mito,numt,uniform)");

  
    try
    {
        arg_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }


    if(outPrefix.empty()) {
    outPrefix = "simulated";
    }
    
    std::filesystem::path  outNameOne = outPrefix;
    outNameOne += "_R1.fq";
    std::filesystem::path  outNameTwo = outPrefix;
    outNameTwo += "_R2.fq";

    auto reference = FastaWrapper{inFastaName};

    std::vector<std::string> ref_ids;
    std::vector<size_t> ref_lengths;
    std::map<std::string, int> contig_tid_map;
    for( int i=0; i< reference.getNumberContigs(); i++){
      std::string tempString{reference.getContigNameByIndex(i)};
      ref_ids.push_back(tempString);
      size_t tempLength = reference.getLengthByContig(tempString);
      ref_lengths.push_back(tempLength);
      contig_tid_map[tempString] = i;
    }
        
    std::filesystem::path  outBamName = outPrefix;
    outBamName += ".bam";


    seqan3::sam_file_output outBam{outBamName, ref_ids, ref_lengths};

    //out fastqs
    seqan3::sequence_file_output outFastqOne{outNameOne};
    seqan3::sequence_file_output outFastqTwo{outNameTwo};
    // std::ofstream outFastqOne(outNameOne, std::ios::binary);
    // std::ofstream outFastqTwo(outNameTwo, std::ios::binary);
    

    
    long long int read_pair_id_counter = 1;
  //   // takes too long to read in reference. Seeing if I can avoid this
  //   //reference_storage_t ref_storage{};
  //   //read_reference(inFastaName, ref_storage);

    
    
    
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    // std::uniform_int_distribution<> distrib(10,500);
    discreteDistributions::Distribution distParam;
    
    if(distString == "mito"){
      distParam =     discreteDistributions::Distribution::mito;
    }
    else if(distString ==  "numt") {
     distParam =     discreteDistributions::Distribution::numt;
    }
    else if( distString ==  "uniform") {
      distParam =     discreteDistributions::Distribution::uniform;
    }
    else {
      std::cerr << "Invalid distribution parameter" << std::endl;
      exit(1);
    }
  
    
    auto distrib = discreteDistributions(distParam);
    std::uniform_int_distribution<> coinflip(0, 1);
    // iterate through bed file
    std::ifstream inBed;
    inBed.open(bedName);
    Bed bedLine;

    while(inBed >> bedLine){
    std::cout << bedLine << std::endl;
      std::string contig = bedLine.chrom;
     hts_pos_t start = bedLine.start;
     // bed is a start 0 based inclusive and end exclusive, so needs to be one less
     hts_pos_t end = bedLine.end -1;

     //std::cout << bedLine << std::endl;
     hts_pos_t contigLength = reference.getLengthByContig(contig);

     hts_pos_t correctedStart = std::max(start - overhang_length, (hts_pos_t) 0);
     hts_pos_t correctedEnd = std::min(end + overhang_length, contigLength);



    std::string region_string = reference.getRegion(contig,correctedStart, correctedEnd);
    int region_length = region_string.length();
    // std::cout << "string length" << region_length << "calculated length " <<
    // 	      correctedEnd - correctedStart + 1 << std::endl;

    // running total index is added by fragment size each iteration to set the index
    int running_total_index = 0;
    // running total bases is only added to if we're in the region so copy number is corrected
    // aka so copy number only reflects number of actual simulated reads
    int running_total_bases = 0;
    //int current_start=0, current_end=0;
    // integer division for copy number, modulo for index.
    // std::cout << "Draws: " << std::endl;
     int copy_count=0;
    do {

      // start at the beginning of the contig
      // draw random fragment size from distribution
      int draw = distrib.draw(gen);
      int current_start = running_total_index % region_length;
      // could theoretically go past the end of string depending on parameters
      int current_end = current_start + draw;

      running_total_index += draw;
      copy_count = running_total_bases / region_length;

      if(current_end >= region_length) {
	continue;
      }
      else {
	running_total_bases += draw;
      }

      // std::cout << "Running total: " << running_total_bases
      // 		<< " Length: " << region_length
      // 		<< "Start " << current_start
      // 		<< " End " << current_end
      // 		<< " Diff " << current_end - current_start
      // 		<< " Copy Number: " << copy_count << std::endl;
      hts_pos_t insert_length = current_end - current_start + 1 ;
      

      //write out a read of this size. end becomes next iteration's beginning
      seqan3::dna5_vector insertVector =
	stringToDNA(region_string.substr(current_start, insert_length));

      // randomly choose sense or antisense fragment
      auto corrected_start_full_contig = correctedStart + current_start;
      auto corrected_end_full_contig = corrected_start_full_contig + insert_length;
      int sense_draw = coinflip(gen);
      
      std::pair<seqan3::dna5_vector,seqan3::dna5_vector> read_pair;
      std::stringstream read_name;
      read_name << contig << ":" << corrected_start_full_contig << "-"
		<< corrected_end_full_contig;
      // initialize read flags
      seqan3::sam_flag read_one_flag = seqan3::sam_flag::paired
	| seqan3::sam_flag::proper_pair
	| seqan3::sam_flag::first_in_pair;

      seqan3::sam_flag read_two_flag = seqan3::sam_flag::paired
	| seqan3::sam_flag::proper_pair
	| seqan3::sam_flag::second_in_pair;

            // output bam file
      using bam_types = seqan3::type_list<std::vector<seqan3::dna5>,
					  std::string,
					  std::string,
					  std::string,
					  size_t,
					  std::vector<seqan3::cigar>,
					  size_t,
					  seqan3::sam_flag,
					  std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>>;



      using bam_fields = seqan3::fields<seqan3::field::seq,
					seqan3::field::id,
					seqan3::field::qual,
					seqan3::field::ref_id,
					seqan3::field::ref_offset,
					seqan3::field::cigar,
					seqan3::field::mapq,
					seqan3::field::flag,
					seqan3::field::mate>;
    using bam_record_type = seqan3::sam_record<bam_types, bam_fields>;

    // declare bam specific vals
	size_t ref_offset_first;
	size_t ref_offset_second;
	int32_t mate_insert_first;
	int32_t mate_insert_second;
	size_t map_qual_first, map_qual_second;
	seqan3::dna5_vector read_one_seq_bam, read_two_seq_bam;
      switch(sense_draw){
      case 0:
	// sense 
	read_pair = createPairedEndReads(insertVector, read_length);
	read_one_flag |= seqan3::sam_flag::mate_on_reverse_strand;
	read_two_flag |= seqan3::sam_flag::on_reverse_strand;
	read_name  << "_F1R2" << "_" << read_pair_id_counter;

	// bam fields
	ref_offset_first = corrected_start_full_contig;
	ref_offset_second = corrected_end_full_contig - read_pair.second.size();
	mate_insert_first = insert_length;
	mate_insert_second = -1*insert_length;
	map_qual_first = 60;
	map_qual_second = 60;

	// read pair to this point is the "off-the sequencer true R1 & R2,
	// bam format wants them both in 5'-3' from perspective of the reference genome,
	// creating views for the bam
	 // = insertVector | std::views::reverse | seqan3::views::complement| std::ranges::to<std::vector>();;
	read_one_seq_bam = read_pair.first;
	read_two_seq_bam = read_pair.second  | std::views::reverse | seqan3::views::complement| std::ranges::to<std::vector>();

	break;
	
      case 1:
	// antisense
	// auto insertRCView = insertVector | seqan3::views::complement | std::views::reverse;
	auto insertRCView = insertVector | std::views::reverse | seqan3::views::complement| std::ranges::to<std::vector>();;
	read_pair = createPairedEndReads(insertRCView ,
	 				      read_length);
	read_one_flag |= seqan3::sam_flag::on_reverse_strand;
	read_two_flag |= seqan3::sam_flag::mate_on_reverse_strand;
	read_name  << "_F2R1" << "_" << read_pair_id_counter;

	// bam fields
	ref_offset_second = corrected_start_full_contig;
	ref_offset_first = corrected_end_full_contig - read_pair.first.size();
	mate_insert_second = insert_length;
	mate_insert_first = -1*insert_length;
	map_qual_first = 60;
	map_qual_second = 60;

	read_one_seq_bam = read_pair.first | std::views::reverse | seqan3::views::complement| std::ranges::to<std::vector>();
	read_two_seq_bam = read_pair.second;

	
	break;
      }
      // increment read name counter
      read_pair_id_counter++;
      
      
      // output fastq files
      using fastq_types = seqan3::type_list<std::vector<seqan3::dna5>, std::string, std::string>;
      using fastq_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::qual>;
      using fastq_record_type = seqan3::sequence_record<fastq_types, fastq_fields>;

      std::string first_read_qualities(read_pair.first.size(), ']');
      fastq_record_type outFastqOneRecord{read_pair.first, read_name.str(), first_read_qualities};
      outFastqOne.push_back(outFastqOneRecord);

      std::string second_read_qualities(read_pair.second.size(), ']');
      fastq_record_type  outFastqTwoRecord{read_pair.second, read_name.str(), second_read_qualities};
      outFastqTwo.push_back(outFastqTwoRecord);

      // output bam file
      using bam_types = seqan3::type_list<std::vector<seqan3::dna5>,
					  std::string,
					  std::string,
					  std::string,
					  size_t,
					  std::vector<seqan3::cigar>,
					  size_t,
					  seqan3::sam_flag,
					  std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>>;



      using bam_fields = seqan3::fields<seqan3::field::seq,
					seqan3::field::id,
					seqan3::field::qual,
					seqan3::field::ref_id,
					seqan3::field::ref_offset,
					seqan3::field::cigar,
					seqan3::field::mapq,
					seqan3::field::flag,
					seqan3::field::mate>;
    using bam_record_type = seqan3::sam_record<bam_types, bam_fields>;

    	// make bam elements
	seqan3::cigar cigar_match_first{(unsigned int) read_pair.first.size(),'M'_cigar_operation};
	std::vector<seqan3::cigar> cigar_string_first{cigar_match_first};

	std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> read_one_next;
	read_one_next = std::make_tuple(contig_tid_map[contig], ref_offset_second, mate_insert_first);

	seqan3::cigar cigar_match_second{(unsigned int) read_pair.second.size(),'M'_cigar_operation};
	std::vector<seqan3::cigar> cigar_string_second{cigar_match_second};

	std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> read_two_next;
	read_two_next = std::make_tuple(contig_tid_map[contig], ref_offset_first, mate_insert_second);

    
    bam_record_type outBamRecordR1{
		     read_one_seq_bam,
				   read_name.str(),
				   first_read_qualities,
				   contig,
		                   ref_offset_first,
				   cigar_string_first,
				   map_qual_first,
				   read_one_flag,
				   read_one_next};
    outBam.push_back(outBamRecordR1);

    bam_record_type outBamRecordR2{
		     read_two_seq_bam,
				   read_name.str(),
				   second_read_qualities,
				   contig,
		                   ref_offset_second,
				   cigar_string_second,
				   map_qual_second,
				   read_two_flag,
				   read_two_next};
    outBam.push_back(outBamRecordR2);
    
     // check if we are past the end, if so advance "copy_count" and set at beginning + how much you went over
      //
    }
    while(copy_count <= target_coverage - 1);
    // std::cout << std::endl;
    }
    return 0;
    
}
