
//#include <bits/getopt_ext.h>
#include "seqan3/argument_parser/argument_parser.hpp"
#include "seqan3/argument_parser/validators.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <stdlib.h>
//#include <getopt.h>

#include <bits/stdc++.h>
//#include <bits/getopt_ext.h>
//#include <getopt-gnu.h

#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

#include <gsl/gsl_statistics_int.h>

#include <algorithm>
#include <cmath>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
namespace fs = std::filesystem;
struct Bed{
  std::string chrom;
  long long int start;
  long long int end;
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

int cmpfunc (const void * a, const void * b) {
  //return ( *(double*)a - *(double*)b );
   if (*(double*)a > *(double*)b) return 1;
   else if (*(double*)a < *(double*)b) return -1;
   else return 0;
}
// void PrintHelp(){
// std::cout <<
// 	   "--bed or -b Input Bed file\n"
// 	   "--prefix or -p Prefix for output\n"
//            "--bam or -s (for sam) for alignments\n"
//            "--outDir or -d for output directory\n"
// -           "--name or -n for sample name\n"
//            "--help or -h Print help\n";
//  exit(1);

// }

int modeValue(const std::vector<int> & inVector) {
  std::map<int,int> count;
  //for(int i = 0; i < inVector.size(); i++) {
  //}
  // put values into a map of intert size and count
  for( const auto & value : inVector ){
    if(count.contains(value) == true) {
      count[value] += 1;
    }
    else {
      count[value] = 0;
    }
  }
  // now that counts have been completed, create vector of pairs of insert size/count
  std::vector<std::pair<int,int>> countPairVec;
  for( const auto & [key, value] : count ) {
    // std::cout << "[" << key << "] = " << value << std::endl;
    countPairVec.push_back(std::make_pair(key, value));
  }
  
  // sort vector of pairs by second value, which is count
  
  std::sort(countPairVec.begin(), countPairVec.end(), [](auto &left, auto &right) {
    return left.second > right.second ;
  });
  
  return countPairVec[0].first;
}

int readFilter(bam1_t *inRead){
  //picking reads that are leftmost so pos + insert
  //gives the fragment span
  //originally did 99 and 163, but BWA puts reads with
  //inserts ~ 300bp as not being properly paired
  //I will remove the "2" bit and check that the tid
  //and mate tid match. Leaves 97 and 161
  // read one of F1R2
  bool isRead1 = ((inRead->core.flag & 97) == 97);
  // read 2 of F2R1
  bool isRead2 = ((inRead->core.flag & 161) == 161);
  if(inRead->core.tid != inRead->core.mtid){
	  return 0;
  }
  else if((inRead->core.qual >= 20) && isRead1){
    return 1;
      }
  else if((inRead->core.qual >= 20) && isRead2){
    return 2;
  }
  else{
  return 0;
  }
}

void updateCovArray(int inArray[], long long int start, long long int end){
  for(long long int i=start; i <= end; i++){
    inArray[i] += 1;
  }
}

//void initialize_arg_parser(
int main(int argc, char* argv[]){
  std::string filePrefix;
  std::string bedName;
  std::string bamName;
  std::string sampleName;
  std::string outDirectoryName;

  seqan3::argument_parser arg_parser{"Inserts-and-Coverage", argc, argv, seqan3::update_notifications::off};

  arg_parser.add_option(filePrefix, 'p', "prefix", "Prefix",
			seqan3::option_spec::required);
  arg_parser.add_option(bedName, 'b', "bed", "Bed file name (bed must be 6 column format",
			seqan3::option_spec::required, seqan3::input_file_validator{{"bed"}});
  arg_parser.add_option(bamName, 's', "bam", "Bam file name",
			seqan3::option_spec::required, seqan3::input_file_validator{{"bam"}});
  arg_parser.add_option(sampleName, 'n', "name", "Sample name",
			seqan3::option_spec::required);
  arg_parser.add_option(outDirectoryName, 'd', "outDir", "Output directory name",
			seqan3::option_spec::required, seqan3::output_directory_validator{});

  try
    {
         arg_parser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "[Error parsing arguments] " << ext.what() << "\n"; // customise your error message
        return -1;
    }
  // const char* const short_opts = "p:b:s:d:h";
  // const option long_opts[] = {	
  //   {"prefix", required_argument, nullptr, 'p'},
  //   {"bed", required_argument, nullptr, 'b'},
  //   {"bam", required_argument, nullptr, 's'},
  //   {"name", required_argument, nullptr, 'n'},
  //   {"outDir", required_argument, nullptr, 'd'},
  //   {"help", no_argument, nullptr, 'h'},
  //   {0,0,0,0}
  // };


  //  while(true){
  //    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
  //    if(opt == -1){
  //      break;
  //    }
  //    switch(opt)
  //      {
  //      case 'p':
  // 	 filePrefix = std::string(optarg);
  // 	 break;
  //      case 'b':
  // 	 bedName = std::string(optarg);
  // 	 break;
  //      case 's':
  // 	 bamName = std::string(optarg);
  // 	 break;
  //      case 'n':
  // 	 sampleName = std::string(optarg);
  // 	 break;
  //      case 'd':
  // 	 outDirectoryName = std::string(optarg);
  // 	 break;
  //      case 'h':
  // 	 PrintHelp();
  // 	 break;
  //      default:
  // 	 PrintHelp();
  // 	 break;
  //      }
  //   }
   // check for required parameters
   // if(filePrefix == ""){
   //   std::cerr << "no output prefix provided" << std::endl;
   //   exit(1);
   // }
   // if(bedName == ""){
   //   std::cerr << "no bed file provided" << std::endl;
   //   exit(1);
   // }
   // if(bamName == ""){
   //   std::cerr << "no bam file provided" << std::endl;
   //   exit(1);
   // }
   // if(sampleName == ""){
   //   std::cerr << "no sample name provided" << std::endl;
   //   exit(1);
   // }

   // iterate through bed file
   std::ifstream inBed;
   inBed.open(bedName);
   Bed bedLine;

   if(outDirectoryName == ""){
     outDirectoryName = ".";
   }
   fs::path outDirectory{outDirectoryName};
   if(! fs::exists(outDirectory)){
     fs::create_directory(outDirectory);   
   }
   fs::path coverageDir{"coverage-outputs"};
   if(! fs::exists(outDirectory / coverageDir))
      fs::create_directory(outDirectory / coverageDir);
   
   fs::path insertDir{"insert-outputs"};
      if(! fs::exists(outDirectory / insertDir))
	fs::create_directory(outDirectory / insertDir);
   // open output file
   std::ofstream outCovFile;
   
   fs::path outCovName{sampleName + "." + filePrefix + ".coverage.tsv"};
   outCovFile.open(outDirectory / coverageDir / outCovName  , std::ios::out);
   // todo: add columns for number of start and number of ends (forward and reverse)
   outCovFile << "Chromosome\tStart\tEnd\tName\tScore\tStrand\tSample\tOffset\tForward_Depth\t" <<
     "Reverse_Depth\tDepth\tForward_Starts\tForward_Ends\tReverse_Starts\tReverse_Ends" << std::endl;
   std::vector<int> covVec;

   std::ofstream outInsertHistFile;
   
   fs::path outInsertHistName{sampleName + "." + filePrefix + ".insert.hist.tsv"};
   outInsertHistFile.open(outDirectory / insertDir /  outInsertHistName, std::ios::out);
   outInsertHistFile << "Sample\tIntervals\tRaw Density\t% Density\tSmoothed % Density" << std::endl;
   // read in bam

   samFile * inBam = sam_open(bamName.c_str(), "r");
   bam_hdr_t * bamHeader = sam_hdr_read(inBam);
   hts_idx_t * inBamIndex = bam_index_load(bamName.c_str());

    long long int runningInsertTotal = 0;
       long long int insertCount = 0;
     auto maxInsert = 1000;
     long long int rawInsertHist[maxInsert];
     double percentDensityInsert[maxInsert];
     double smoothedPercentDensity[maxInsert];
     for(int i=0; i < maxInsert; i++){
       rawInsertHist[i] = 0;
       percentDensityInsert[i] = 0;
       smoothedPercentDensity[i] = 0;
     }
     std::vector<int> insertVec;

     auto insertVecBinSize = 100;
     std::ofstream outBinInsertFile;
   
   fs::path outBinInsertName{sampleName + "." + filePrefix + ".insert.bin.tsv"};
   outBinInsertFile.open(outDirectory / insertDir / outBinInsertName  , std::ios::out);
   outBinInsertFile << "Chromosome\tStart\tEnd\tName\tScore\tStrand\tBin\tSample\tMean\tStd-Dev\tMin\t25%\t50%\t75%\tMax\tMode" << std::endl;
     
   while(inBed >> bedLine){
     //std::cout << bedLine << std::endl;
     long int start = bedLine.start;
     // bed is a start 0 based inclusive and end exclusive, so needs to be one less
     long int end = bedLine.end -1;
     
     auto intervalLength = end - start + 1;
     // implemented these as primitive arrays out of convenience, wouldn't have to keep reallocating memory and could easily use gnu scientific library (although should work with vector as well)
     int FoneRtwo[intervalLength];
     int FtwoRone[intervalLength];

     int FoneRtwoStarts[intervalLength];
     int FoneRtwoEnds[intervalLength];

     int FtwoRoneStarts[intervalLength];
     int FtwoRoneEnds[intervalLength];


     
     for(auto i=0; i < intervalLength; i++){
       FoneRtwo[i] = 0;
       FtwoRone[i] = 0;

       FoneRtwoStarts[i] = 0;
       FoneRtwoEnds[i] = 0;

       FtwoRoneStarts[i] = 0;
       FtwoRoneEnds[i] = 0;
     }
     
     // 2d array indexed by start, end and value is count. So each Bed interval will have the count of start/stop
     // these being 2d arrays appears to blow stack.
     // std::vector<std::vector<int>> fragmentPosCounts(intervalLength);
     // for(auto i=0; i < intervalLength; i++){
     //   //auto tempVector = fragmentPosCounts[i];
     //   //fragmentPosCounts[i].resize(intervalLength);
     //   fragmentPosCounts[i] = std::vector<int>(intervalLength,0);
     //   // for(auto j=0; j < intervalLength; j++){
     //   // 	 // tempVector[j] = 0;
     //   // 	 std::cout << fragmentPosCounts[i][j] << "\t";
     //   // }
     // }
     // std::cout << std::endl;
     
     // make a vector of vectors
     auto numberInsertBins = ((intervalLength - 1) / insertVecBinSize) + 1;
     // this is ugly b/c primitive array of vector<int>'s. Eventually should make vector of vectors.
     // probably not as pressing a need as "fragmentPosCounts" because this will only be # bin dimension array
     std::vector<int> binInsertVecs[numberInsertBins];
     
     
     
  // iterate over region
  // I'm retreiving reads from upstream of the region to make sure
     //I get more for the F2R1 orientations. Has the advantage that
     // I can make the minimum 0 whereas the converse I wouldn't
     //know the range of contig
     auto upstreamStart = std::max(bedLine.start - 300, (long long int) 0);
  int contigCode = bam_name2id (bamHeader, bedLine.chrom.c_str());
  hts_itr_t* bamItr =  sam_itr_queryi(inBamIndex, contigCode, upstreamStart, bedLine.end);
  if(bamItr == nullptr){
    continue;
  }

  // to iterate over reads
  bam1_t *bam_rec = bam_init1();
  
  while(sam_itr_next(inBam, bamItr, bam_rec) > 0){
    int strandReturn = readFilter(bam_rec);
    if( strandReturn > 0){
      int insertSize = bam_rec->core.isize;
      auto absInsert = abs(insertSize);
      
      long int leftPos = bam_rec->core.pos;
      long int rightPos = bam_rec->core.pos + insertSize;

      /*if((rightPos <= start) || (leftPos >= end)){
	bam_destroy1(bam_rec);
        bam_rec = bam_init1();
	continue;
	}*/
      // checks that the read is within the interval and not 0 insert size 
      if((leftPos  < end) && (rightPos >  start) &&
	  (insertSize > 0)){
	//std::cout << "bed start " << start << " bed end " << end << std::endl;
	//std::cout << "read start " << leftPos << " insert " << insertSize << " end " << rightPos << std::endl;
	
	auto covArrayStart = std::max(leftPos, start) - start;
        auto covArrayEnd = std::min(rightPos, end) - start ;
        //std::cout << "corrected start " << covArrayStart << "corrected end " << covArrayEnd << std::endl;

	// coverages
      if(strandReturn == 1){
	updateCovArray(FoneRtwo, covArrayStart, covArrayEnd);
	FoneRtwoStarts[covArrayStart] += 1;
        FoneRtwoEnds[covArrayEnd] += 1;

	// forward orientation, so start is leftmost
	binInsertVecs[covArrayStart / insertVecBinSize].push_back(absInsert);
      }
      else if(strandReturn == 2){
	updateCovArray(FtwoRone, covArrayStart, covArrayEnd);
	// since I care about beginning of fragment vs. end, need to flip start and end.
	FtwoRoneStarts[covArrayEnd] += 1;
        FtwoRoneEnds[covArrayStart] += 1;

	// reverse orientation, so start is rightmost
	binInsertVecs[covArrayEnd / insertVecBinSize].push_back(absInsert);
      }
      }

      // for inserts I'm just taking read 1s
      if(((bam_rec->core.flag & 64) == 64) &&
	(bam_rec->core.qual >= 20) && (absInsert >0) &&
		   (absInsert <= maxInsert)){
	runningInsertTotal += absInsert;
	insertCount += 1;
     
	rawInsertHist[absInsert -1]+=1;
	insertVec.push_back(absInsert);

	// add to count of that # of inserts
	// fragmentPosCounts[leftPos][rightPos] += 1;
	
    }
      } 
	
      
      
      
    bam_destroy1(bam_rec);
    bam_rec = bam_init1();
    }
  bam_itr_destroy(bamItr);
  
  // went through the interval, so writing out to file
  for(int offset = 0; offset < intervalLength; offset++){
    auto fullCoverage = FoneRtwo[offset] + FtwoRone[offset];
    covVec.push_back(fullCoverage);
    outCovFile << bedLine.chrom << "\t" << bedLine.start << "\t" << bedLine.end << "\t" <<
      bedLine.name << "\t" << bedLine.score << "\t" << bedLine.strand << "\t" <<
      sampleName << "\t" << offset << "\t" <<     
      FoneRtwo[offset] << "\t" << FtwoRone[offset] << "\t" << fullCoverage << "\t" <<
      FoneRtwoStarts[offset] << "\t" << FoneRtwoEnds[offset] << "\t"
	       << FtwoRoneStarts[offset] << "\t" << FtwoRoneEnds[offset]  << std::endl;
  }
  // writing interval binned inserts to file
  for(int binNum=0; binNum < numberInsertBins; binNum++){
    auto tempInsertVec = binInsertVecs[binNum];   
  if(tempInsertVec.size() > 0){
   std::sort(tempInsertVec.begin(), tempInsertVec.end());
   auto insertMean = gsl_stats_int_mean(tempInsertVec.data(), 1, tempInsertVec.size());
   auto insertStdev = gsl_stats_int_sd(tempInsertVec.data(), 1, tempInsertVec.size());
   auto insertMin = gsl_stats_int_min(tempInsertVec.data(), 1, tempInsertVec.size());
   auto insertTwentyFive = gsl_stats_int_quantile_from_sorted_data(tempInsertVec.data(), 1, tempInsertVec.size(), 0.25);
   auto insertMedian = gsl_stats_int_median_from_sorted_data(tempInsertVec.data(), 1, tempInsertVec.size());
   auto insertSeventyFive = gsl_stats_int_quantile_from_sorted_data(tempInsertVec.data(), 1, tempInsertVec.size(), 0.75);
   auto insertMax = gsl_stats_int_max(tempInsertVec.data(), 1, tempInsertVec.size());
   auto insertMode = modeValue( insertVec);
   outBinInsertFile << bedLine.chrom << "\t" << bedLine.start << "\t"
		    << bedLine.end << "\t" << bedLine.name << "\t"
		    << bedLine.score << "\t" << bedLine.strand << "\t" 
		    <<  binNum + 1 << "\t" << sampleName << "\t"
		    << insertMean << "\t" << insertStdev << "\t"
		    << insertMin << "\t" << insertTwentyFive
		    << "\t" << insertMedian << "\t" <<
     insertSeventyFive << "\t" << insertMax << "\t" << insertMode <<std::endl;
   }
   else{
     outBinInsertFile << bedLine.chrom << "\t" << bedLine.start << "\t"
                      << bedLine.end << "\t" << bedLine.name << "\t"
                      << bedLine.score << "\t" << bedLine.strand << "\t"
                      << binNum + 1 << "\t" << sampleName << "\t"
                      << "-nan\t-nan\t-nan\t-nan\t-nan\t-nan\t-nan\t-nan" << std::endl;
   }
  }

   }

   // calculate % density and smoothed % density
   for(int outInsertSize=0; outInsertSize < maxInsert; outInsertSize++){
     auto tempInsert = rawInsertHist[outInsertSize];
     percentDensityInsert[outInsertSize] = (tempInsert / (double) insertCount ) * 100.0;
     //std::cout <<  percentDensityInsert[outInsertSize] << std::endl;
   }
   int smoothBinSize = 5;
   int medianIndex = smoothBinSize / 2;
   for(int i=0; i < (maxInsert - smoothBinSize + 1); i++){
     //std::cout << "insert index " << i << std::endl;
     double tempArray[smoothBinSize];
     for(int j=0; j < smoothBinSize; j++){
       tempArray[j] = percentDensityInsert[i + j];
       //std::cout << tempArray[j] << " ";
     }
     //std::cout << std::endl << "After sort" << std::endl;
     // sort and find median
     qsort(tempArray, smoothBinSize, sizeof(double), cmpfunc);
     auto medianVal = tempArray[medianIndex];
     smoothedPercentDensity[i + medianIndex] = medianVal;
     // for(int k=0; k < smoothBinSize; k++){
     //   std::cout << tempArray[k] << " ";
     // }
     // std::cout << std::endl;

   }
   // write out insert file
   for(int outInsertSize=0; outInsertSize < maxInsert; outInsertSize++){
     //auto tempInsert = rawInsertHist[outInsertSize];
     //percentDensityInsert[outInsertSize] = (tempInsert / (double) insertCount ) * 100.0;
     
       outInsertHistFile << sampleName << "\t" << outInsertSize + 1 << "\t"
		     << rawInsertHist[outInsertSize] << "\t" 
		     << percentDensityInsert[outInsertSize] << "\t"
		     << smoothedPercentDensity[outInsertSize] << std::endl;
     }
   
   
   //output stats coverage
   std::ofstream outCovStatsFile;
   fs::path outCovStatsName{sampleName + "." + filePrefix + ".coverage.stats.tsv"};
   outCovStatsFile.open(outDirectory / coverageDir / outCovStatsName  , std::ios::out);

   outCovStatsFile << "Sample\tMean\tStd-Dev\tMin\t25%\t50%\t75%\tMax" << std::endl;
   if(covVec.size() > 0){
   std::sort(covVec.begin(), covVec.end());
   
   auto covMean = gsl_stats_int_mean(covVec.data(), 1, covVec.size());
   auto covStdev = gsl_stats_int_sd(covVec.data(), 1, covVec.size());
   auto covMin = gsl_stats_int_min(covVec.data(), 1, covVec.size());
   auto covTwentyFive = gsl_stats_int_quantile_from_sorted_data(covVec.data(), 1, covVec.size(), 0.25);
   auto covMedian = gsl_stats_int_median_from_sorted_data(covVec.data(), 1, covVec.size());
   auto covSeventyFive = gsl_stats_int_quantile_from_sorted_data(covVec.data(), 1, covVec.size(), 0.75);
   auto covMax = gsl_stats_int_max(covVec.data(), 1, covVec.size());
   
   
   outCovStatsFile << sampleName << "\t" << covMean << "\t" << covStdev << "\t" <<
     covMin << "\t" << covTwentyFive << "\t" << covMedian << "\t" << covSeventyFive <<
     "\t" << covMax << std::endl;
   }
   else {
     outCovStatsFile << sampleName << "\t-nan\t-nan\t-nan\t-nan\t-nan\t-nan\t-nan" << std::endl;
   }
   // output stats inserts
   std::ofstream outInsertStatsFile;
   fs::path outInsertStatsName{sampleName + "." + filePrefix + ".insert.stats.tsv"};
   outInsertStatsFile.open(outDirectory / insertDir /  outInsertStatsName, std::ios::out);

   outInsertStatsFile << "Sample\tMean\tStd-Dev\tMin\t25%\t50%\t75%\tMax\tMode" << std::endl;
   if(insertVec.size() > 0){
   std::sort(insertVec.begin(), insertVec.end());
   auto insertMean = gsl_stats_int_mean(insertVec.data(), 1, insertVec.size());
   auto insertStdev = gsl_stats_int_sd(insertVec.data(), 1, insertVec.size());
   auto insertMin = gsl_stats_int_min(insertVec.data(), 1, insertVec.size());
   auto insertTwentyFive = gsl_stats_int_quantile_from_sorted_data(insertVec.data(), 1, insertVec.size(), 0.25);
   auto insertMedian = gsl_stats_int_median_from_sorted_data(insertVec.data(), 1, insertVec.size());
   auto insertSeventyFive = gsl_stats_int_quantile_from_sorted_data(insertVec.data(), 1, insertVec.size(), 0.75);
   auto insertMax = gsl_stats_int_max(insertVec.data(), 1, insertVec.size());
   auto insertMode = modeValue( insertVec);
   outInsertStatsFile << sampleName << "\t" << insertMean << "\t" << insertStdev <<
     "\t" << insertMin << "\t" << insertTwentyFive << "\t" << insertMedian << "\t" <<
     insertSeventyFive << "\t" << insertMax << "\t" << insertMode <<std::endl;
   }
   else{
     outInsertStatsFile << sampleName << "\t-nan\t-nan\t-nan\t-nan\t-nan\t-nan\t-nan\t-nan" << std::endl;
   }
   hts_idx_destroy(inBamIndex);
   bam_hdr_destroy(bamHeader);
   inBed.close();
   outCovFile.close();
   outInsertHistFile.close();
}
