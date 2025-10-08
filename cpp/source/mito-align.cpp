#include "htslib/kstring.h"
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <stdio.h> /*for printf */
#include <getopt.h> /* for getopt_long */
#include <stdlib.h> /* for exit */
#include <regex.h>

#include <htslib/kseq.h>
#include <htslib/hts.h>
#include <htslib/sam.h>



/* ******* Read Name ******** */
#define MAX_READNAME_SIZE 256 - 2
#define MAX_READ_SIZE 65536
typedef struct {

  char name[MAX_READNAME_SIZE];
  char bases[MAX_READ_SIZE];
  unsigned base_qualities[MAX_READ_SIZE];
  // these are basically pointers. 0 means read1 is invalid. MAX_READ_SIZE means read2 is invalid. 
  unsigned short int read_1_end{0};
  unsigned short int read_2_begin{(unsigned short int)MAX_READ_SIZE};
} read_info_t;

/* ******** End Read Name ********** */

typedef struct {
  int mapping_qual=30;
  char* input_filename= NULL;
  char* read_two_filename = NULL;
} args_t;

args_t parse_args(int argc, char* argv[]){
  args_t parsed_params = args_t{};

  int c;
  int digit_optind = 0;

  while (1) {
               int this_option_optind = optind ? optind : 1;
               int option_index = 0;
               static struct option long_options[] = {
		 {"input",     required_argument, NULL,  'i' },
		 {"R2",     required_argument, NULL,  'm' },
                   {"quality",  required_argument,  0,  'q' },
                   // {"delete",  required_argument, 0,  0 },
                   {"verbose", no_argument,       0,  0 },
                   {0,         0,                 0,  0 }
               };

               c = getopt_long(argc, argv, "i:q:m:",
                               long_options, &option_index);
	       printf("c=%d, option_index=%d\n", c, option_index);
               if (c == -1)
                   break;
	       switch(c){
	       case 'i':
		 printf("input\n");
		 printf("input file name: %s\n", optarg);
		 parsed_params.input_filename = optarg;
		 break;
	       case 'm':
		 printf("read2\n");
		 printf("read2 file name: %s\n", optarg);
		 parsed_params.read_two_filename = optarg;
		 break;
	       case 'q':
		 printf("qual\n");
		 printf("qual: %d\n", atoi(optarg));
		 break;
	       case '?':
		 printf("unkown option\n");
		 break;
	       default:
		 printf("error parsing argument");
		 exit(EXIT_FAILURE);
    }
  }

  
  return parsed_params;
}

int parse_fastq_entry(htsFile * fp,
		      read_info_t* read,
		      int paired_end_flag) {
  // a properly formed fastq should be 4 lines

  // paired_end_flag, 0 for first (or only read), -1 for second read in pair 
  kstring_t line_one_string;
  ks_initialize(&line_one_string);
  int line_one_length = hts_getline(fp, KS_SEP_LINE, &line_one_string);
  size_t readname_length{0};
  
  char* readname_delimiter_p = NULL;
  // printf("%s\n", line_one_string.s);

  /* ************ 1st Fastq Line ************** */
  // first line starts with '@', and read name is everything b/w '@' & first space
  int read_name_end = 0;

  // if end of file (-1) or not starting with '@', it's an invalid first line entry
  if((-2 >= line_one_length) || ('@' != line_one_string.s[0])) {
    fprintf(stderr, "Invalid read name line from fastq '%s'", line_one_string.s);
    // other requirements for read name. regex [!-?A-~]{1,254} seems like a good guideline. They also include a more extensive regex in the paper.
    exit(EXIT_FAILURE);
  }
  else if (-1 == line_one_length){
      return 0; // -1 means it was at EOF, which should only ocurr for a 1st fastq line
    }
  else {
    readname_delimiter_p = strchr(line_one_string.s, ' ');
    if(NULL == readname_delimiter_p) {

      // actually, valid read names do not need the extra info. Just going to skip if no extras.
          // fprintf(stderr, "Malformed read name line from fastq '%s'", line_one_string.s);
	  // exit(EXIT_FAILURE);
      readname_length = (size_t) line_one_string.s - 1; // should be index of last character of read name
    }
    else {
      readname_length = readname_delimiter_p - line_one_string.s; //  -1 (to be last letter, not delimiter) + 1 (account for 0 indexing);
      // printf("Read name length: %lu, Read name: '%.*s'\n", readname_length -1,
	     // (int) readname_length -1, line_one_string.s + 1);
      // could add a check to make sure readname_length is never larger than MAX_READNAME_SIZE
      
      //printf("Read name after copying in parse function: %s\n", read->name);
    }
    strncpy(read->name, line_one_string.s + 1, readname_length);
  }

  /* ************* 2nd Fastq Line ******************* */
  kstring_t line_two_string;
  ks_initialize(&line_two_string);
  int line_two_length = hts_getline(fp, KS_SEP_LINE, &line_two_string);
  // printf("%s\n", line_two_string.s);

  strncpy(read->bases, line_two_string.s, line_two_length);
  // printf("Bases after copying %s\n", read->bases);
    /* ************* 3rd Fastq Line ******************* */
  kstring_t line_three_string;
  ks_initialize(&line_three_string);
  int line_three_length = hts_getline(fp, KS_SEP_LINE, &line_three_string);
  // printf("%s\n", line_three_string.s);

  if((-1 >= line_three_length) || ('+' != line_three_string.s[0])) {
    fprintf(stderr, "Invalid separator line from line three of fastq '%s'", line_three_string.s);
    exit(EXIT_FAILURE);
  }

      /* ************* 4th Fastq Line ******************* */
  kstring_t line_four_string;
  ks_initialize(&line_four_string);
  int line_four_length = hts_getline(fp, KS_SEP_LINE, &line_four_string);
  // printf("%s\n", line_four_string.s);

  if(-1 >= line_four_length) {
    fprintf(stderr, "Invalid read quality line from line four of fastq '%s'", line_four_string.s);
    exit(EXIT_FAILURE);
  }
  else{
    //paired_end_flag when 0, starts array at beginning. When 1, starts array at max size - read length
    int i{0};
      while( i < line_four_length){
      read->base_qualities[i + paired_end_flag*MAX_READ_SIZE -
			   paired_end_flag*line_four_length] = (unsigned)line_four_string.s[i] -33;
      // printf("%c", line_four_string.s[i]);
      // printf("%u",read->base_qualities[i + paired_end_flag*MAX_READ_SIZE -
      // 			   paired_end_flag*line_four_length]);
      i++;
    }
    //read->base_qualities[line_four_length] = "\0";
      // printf("Base Qualities after copying %.u\n", (int) readname_length -1,read->base_qualities);

  }

  // free up kstring resources. Long term not allocating and freeing these resources is probably ideal.
  ks_free(& line_one_string);
  ks_free(& line_two_string);
  ks_free(& line_three_string);
  ks_free(& line_four_string);

  // line two and four should both be the same length of the read
  return line_two_length;
}

// to debug
// DEBUGINFOD_URLS="" gdb ./mito-align
//  run --input ~/science/mitochondrial/rachek_2024/fastqs/LRH01_L0022_R1.fastq.gz

int main (int argc, char *argv[]) {
 
 // printf("\ncmdline args count=%d", argc);

 // /* First argument is executable name only */
 // printf("\nexe name=%s", argv[0]);

  // this is a read which is designed to be large enough to handle any one illumina or nanoport read. Currently I plan to only allocate one read so I don't have to keep allocating memory.
  read_info_t read;
  read_info_t read_two;

 args_t params = parse_args(argc, argv);
 // printf("mapping qual: %d", params.mapping_qual);


 // will first only allow to read in from unsorted or name sorted bams. This is to allow metadata to be added for sample name, flowcell, etc.
 htsFile * infile_p = hts_open (params.input_filename, "r");
 const htsFormat * infile_format  = hts_get_format (infile_p);
 const htsFormatCategory infile_format_cat = infile_format->category;
 const htsExactFormat infile_format_exact =  infile_format->format;
 const htsCompression infile_compression = infile_format->compression;

 htsFile * infile_two_p;
 if(NULL != params.read_two_filename) {
 infile_two_p = hts_open (params.read_two_filename, "r");
 // not currently checking if read 2 is in correct format
 }
 // from htslib, "1" for 'Sequence data -- SAM, BAM, CRAM, etc'. Appears to also recognize fasta and fastq.gz
 printf("file format %d\n", infile_format->category);
 if (1 != infile_format_cat)              
   {
     printf("Not a sequence file\n");
     exit(EXIT_FAILURE);
   }

 // if fastq
 if((htsExactFormat::fastq_format == infile_format_exact)) {
   // fastq file
   // apparently fastq support isn't well supported by htslib. Also, it really is an advantage to include metadata in an unaligned bam.
     printf("Fastq file\n");
     if(htsCompression::no_compression == infile_compression){
       printf("Uncompressed Fastq\n");
     }
     else if(htsCompression::gzip == infile_compression){
       printf("Gzipped compressed Fastq\n");
     }
     else {
       printf("Don't support other fastq compression schemes\n");
     }
     // trying to just parse with htslib functions
     int read_one_written{1};
     while(0 < read_one_written){
       read_one_written = parse_fastq_entry(infile_p, & read ,0);
       printf("Read 1 number bases written %u\n", read_one_written);
     // if(0 < read_one_written ){
     //   read.read_1_end = read_one_written;
     //   printf("Read 1 name: %s, Number of bases %u\n", read.name, read.read_1_end);
     //   }

       if(NULL != params.read_two_filename) {
	 int read_two_written{1};
	 read_two_written = parse_fastq_entry(infile_two_p, & read_two ,0);
	 if(0 < read_two_written){
	 printf("Read name 2 after returning: %s\n", read_two.name);
	 }
     }
     }

     
   }
 else if(0 != ((htsExactFormat::bam | htsExactFormat::sam) & infile_format_exact)){
   // bam or sam file. 
   printf("Bam or Sam file\n");
   // header should work for sam, bam, or cram file
   sam_hdr_t * sam_hdr_p = sam_hdr_read(infile_p);
   // sort order
   //auto sort_order_p = (kstring_t*)malloc(sizeof(kstring_t));
   kstring_t sort_order_p;
   ks_initialize(&sort_order_p);
   sam_hdr_find_tag_hd(sam_hdr_p, "SO", &sort_order_p);
   printf("sort order %s\n", sort_order_p.s);
 }
 
 // now find exact file format. Especially important for paired end reads
 exit(EXIT_SUCCESS);
 }
