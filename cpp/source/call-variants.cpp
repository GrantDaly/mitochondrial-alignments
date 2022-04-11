/*  Adapted by Grant Daly from
    bam_plcmd.c -- mpileup subcommand.

    Copyright (C) 2008-2015, 2019-2021 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>
    Modified by: Grant Daly

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

// #include <config.h>

// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <unistd.h>
#include <ctype.h>
// #include <string.h>
// #include <strings.h>
// #include <limits.h>
// #include <errno.h>
// #include <sys/stat.h>
// #include <getopt.h>
// #include <inttypes.h>
#include "htslib/hts.h"
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
// #include <htslib/klist.h>
// #include <htslib/khash_str2int.h>
// #include <htslib/cram.h>
// #include "samtools.h"
// #include "bedidx.h"
// #include "sam_opts.h"
// #include "bam_plbuf.h"

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <stdlib.h>


#include <seqan3/argument_parser/all.hpp>

#include "sampleBasePileup.hpp"

typedef struct {
    samFile *fp;
    hts_itr_t *iter;
    sam_hdr_t *header;
  //mplp_ref_t *ref;
  //const mplp_conf_t *conf;
} mplp_aux_t;

struct mpileup_params_t {
  std::string refName;
  std::string bamName;
  std::string outFileName;
  std::string regionString;
  int minBaseQ=30;
} ;

// typedef struct {
//     int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, all, rev_del;
//     int rflag_require, rflag_filter;
//     char *reg, *pl_list, *fai_fname, *output_fname;
//     faidx_t *fai;
//     void *bed, *rghash, *auxlist;
//     int argc;
//     char **argv;
//     char sep, empty, no_ins, no_ins_mods, no_del, no_ends;
//     sam_global_args ga;
// } mplp_conf_t;

// typedef struct {
//     samFile *fp;
//     hts_itr_t *iter;
//     sam_hdr_t *h;
//   //mplp_ref_t *ref;
//     const mplp_conf_t *conf;
// } mplp_aux_t;

// int bam_mplp64_auto_test(bam_mplp_t iter, int *_tid, hts_pos_t *_pos, int *n_plp, const bam_pileup1_t **plp)
// {
//     int i, ret = 0;
//     hts_pos_t new_min_pos = HTS_POS_MAX;
//     uint32_t new_min_tid = (uint32_t)-1;
//     for (i = 0; i < iter->n; ++i) {
//         if (iter->pos[i] == iter->min_pos && iter->tid[i] == iter->min_tid) {
//             int tid;
//             hts_pos_t pos;
//             iter->plp[i] = bam_plp64_auto(iter->iter[i], &tid, &pos, &iter->n_plp[i]);
//             if ( iter->iter[i]->error ) return -1;
//             if (iter->plp[i]) {
//                 iter->tid[i] = tid;
//                 iter->pos[i] = pos;
//             } else {
//                 iter->tid[i] = 0;
//                 iter->pos[i] = 0;
//             }
//         }
//         if (iter->plp[i]) {
//             if (iter->tid[i] < new_min_tid) {
//                 new_min_tid = iter->tid[i];
//                 new_min_pos = iter->pos[i];
//             } else if (iter->tid[i] == new_min_tid && iter->pos[i] < new_min_pos) {
//                 new_min_pos = iter->pos[i];
//             }
//         }
//     }
//     iter->min_pos = new_min_pos;
//     iter->min_tid = new_min_tid;
//     if (new_min_pos == HTS_POS_MAX) return 0;
//     *_tid = new_min_tid; *_pos = new_min_pos;
//     for (i = 0; i < iter->n; ++i) {
//         if (iter->pos[i] == iter->min_pos && iter->tid[i] == iter->min_tid) {
//             n_plp[i] = iter->n_plp[i], plp[i] = iter->plp[i];
//             ++ret;
//         } else n_plp[i] = 0, plp[i] = 0;
//     }
//     return ret;
// }



static int mplp_func(void *data, bam1_t *b)
{
  
  //char *ref = "test";
  mplp_aux_t *plp_data = (mplp_aux_t*)data;
  int ret, skip = 0;
    // hts_pos_t ref_len;

  do {
  	if(plp_data->iter){
	  // std::cout << "Has iterator" << std::endl;
	  // std::cout << "Bam Position " << b->core.pos << std::endl;
	  ret = sam_itr_next(plp_data->fp, plp_data->iter, b);
	  // std::cout << "Return value " << ret << std::endl;
	    }
	else
	  {
	    // std::cout << "Does Not Have Iterator" << std::endl;
	    ret = sam_read1(plp_data->fp, plp_data->header, b);
	  }

	if (ret < 0) break;
	// for some reason reads are being incorrectly marked as unmapped. flag&BMA_FUNMAP = 4
	// for now I'm going to remove this check.
	// if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
	//   std::cout << "tid " << b->core.tid << std::endl;
	  
	//   std::cout << "unmap flag " << (b->core.flag&BAM_FUNMAP) << std::endl;
        //     skip = 1;
        //     continue;
        // }
    } while (skip > 0);
    // std::cout << "Returning with this value " << ret << std::endl;
    return ret;
    // do {
    //     int has_ref;
    //     // ret = plp_data->iter? sam_itr_next(plp_data->fp, plp_data->iter, b) : sam_read1(plp_data->fp, plp_data->header, b);
    // 	// std::cout << "Bam Return Value " << ret << std::endl;
    //     if (ret < 0) break;
    // //     // The 'B' cigar operation is not part of the specification, considering as obsolete.
    // //     //  bam_remove_B(b);
    //     if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
    //         skip = 1;
    //         continue;
    //     }
    // 
}
	// todo; filter by read1, read2 if necessary. May be able to deal with this in the pileup processing function
    //     if (ma->conf->rflag_require && !(ma->conf->rflag_require&b->core.flag)) { skip = 1; continue; }
    //     if (ma->conf->rflag_filter && ma->conf->rflag_filter&b->core.flag) { skip = 1; continue; }
    //     if (ma->conf->bed && ma->conf->all == 0) { // test overlap
    //         skip = !bed_overlap(ma->conf->bed, sam_hdr_tid2name(ma->h, b->core.tid), b->core.pos, bam_endpos(b));
    //         if (skip) continue;
    //     }
    //     if (ma->conf->rghash) { // exclude read groups
    //         uint8_t *rg = bam_aux_get(b, "RG");
    //         skip = (rg && khash_str2int_get(ma->conf->rghash, (const char*)(rg+1), NULL)==0);
    //         if (skip) continue;
    //     }
    //     if (ma->conf->flag & MPLP_ILLUMINA13) {
    //         int i;
    //         uint8_t *qual = bam_get_qual(b);
    //         for (i = 0; i < b->core.l_qseq; ++i)
    //             qual[i] = qual[i] > 31? qual[i] - 31 : 0;
    //     }

    //     if (ma->conf->fai && b->core.tid >= 0) {
    //         has_ref = mplp_get_ref(ma, b->core.tid, &ref, &ref_len);
    //         if (has_ref && ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
    //             fprintf(stderr,"[%s] Skipping because %"PRIhts_pos" is outside of %"PRIhts_pos" [ref:%d]\n",
    //                     __func__, (int64_t) b->core.pos, ref_len, b->core.tid);
    //             skip = 1;
    //             continue;
    //         }
    //     } else {
    //         has_ref = 0;
    //     }

    //     skip = 0;
    //     if (has_ref && (ma->conf->flag&MPLP_REALN)) sam_prob_realn(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
    //     if (has_ref && ma->conf->capQ_thres > 10) {
    //         int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
    //         if (q < 0) skip = 1;
    //         else if (b->core.qual > q) b->core.qual = q;
    //     }
    //     if (b->core.qual < ma->conf->min_mq) skip = 1;
    //     else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) skip = 1;


/*
 * Performs pileup
 * @param conf configuration for this pileup
 * @param n number of files specified in fn
 * @param fn filenames
 * @param fn_idx index filenames
 */


int process_pileup(const bam_pileup1_t *pileup, hts_pos_t pos,  mpileup_params_t * const params, SampleBasePileup & outPileup) {

  // check if the bam is reverse bam_is_rev() or mate reverse bam_ismrev()
  // that is, separate F1R2 R1 and R2, and F2R1 R1 and R2

  //start with ignoring indels
  if( ! pileup->indel) {
    //std::cout << "not an indel" << std::endl;
    int baseQual = -1;
    if(  pileup->b->core.l_qseq)
      {
	char bamChar = -1;
	baseQual = bam_get_qual(pileup->b)[pileup->qpos];
    if(baseQual <  params->minBaseQ){
      return 0;
      }
	
    bamChar = toupper(seq_nt16_str[bam_seqi(bam_get_seq(pileup->b), pileup->qpos)]);
    std::string baseStr(1, bamChar);
	// std::cout << bamChar << std::endl;

	int bamFlag = pileup->b->core.flag;
	// if((bamFlag & BAM_FPAIRED & BAM_FPROPER_PAIR & BAM_FMREVERSE & BAM_FREAD1) ==
	// std::cout << "Flag " << bamFlag << std::endl;

	// F1R2 Read 1's
	if((bamFlag & 99) == 99){
	  //std::cout << "F1R2 Read1" << std::endl;
	  outPileup.F1R2R1s[baseStr] += 1;
	}
	// F1R2 Read 2's
	else if((bamFlag & 147) == 147){
	  //std::cout << "F1R2 Read2" << std::endl;
	  outPileup.F1R2R2s[baseStr] += 1;
	}
	// F2R1 Read 1's
	if((bamFlag & 83) == 83){
	  //std::cout << "F2R1 Read1" << std::endl;
	  outPileup.F2R1R1s[baseStr] += 1;
	}
	else if((bamFlag & 163) == 163){
	  //std::cout << "F1R2 Read2" << std::endl;
	  outPileup.F2R1R2s[baseStr] += 1;
	}
	
  }
    else
      {
	// wasn't a base here. Samtools outputs "N", but I'm just going to skip
	return 0;
      }
  }
  
  return 0;
}

void parse_arguments(mpileup_params_t & params, int argc, char * argv[]) {
  seqan3::argument_parser myparser{"Mitochondrial-Pileup-Caller", argc, argv, seqan3::update_notifications::off};

  std::filesystem::path inBamName{};
  myparser.add_option(inBamName,'b',"bam","The input bam file.",
                        seqan3::option_spec::standard, seqan3::input_file_validator{{"bam"}});

  
  std::filesystem::path inFastaName{}; 
  myparser.add_option(inFastaName,'r',"ref","Reference Genome",
		      seqan3::option_spec::standard, seqan3::input_file_validator{{"fa","fasta"}});


    std::filesystem::path outFileName{"mito.pileups.tsv"}; 
    myparser.add_option(outFileName,'o',"out","Output File Name");

  std::string mitoName{"chrM"};
  myparser.add_option(mitoName, 'n', "mitoName", "Mitochondrial Name");
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
	exit(1);
        
    }

    params.bamName = inBamName.string();
    params.refName = inFastaName.string();
    params.outFileName = outFileName.string();
    params.regionString = mitoName;
    params.minBaseQ = 30;
    

}

int main(int argc, char * argv[]){

  
  //std::string reference =  "../tests/mpileup.ref.fa";

  // filename
  //std::string bamName = "../tests/mpileup.1.sam";

  // mpileup_params_t params = {.refName="../tests/mpileup-inputs/susScr11.mito.fa", .bamName= "../tests/mpileup-inputs/sim-variant.bam",
  // 			     .regionString="chrM", .minBaseQ=30};

  mpileup_params_t params;
  parse_arguments(params, argc,argv);

  std::ofstream outFile;
  outFile.open(params.outFileName);
  SampleBasePileup::writeHeader( outFile);
  // number of samples
  int n_samples = 1;
  //initialize for pileup
  auto max_depth = INT_MAX;
    mplp_aux_t **data;
    int i, tid, *n_plp, tid0 = 0;
    hts_pos_t pos, beg0 = 0, end0 = HTS_POS_MAX, ref_len;
    int minBaseQ = 30;
    const bam_pileup1_t **plp;
    // mplp_ref_t mp_ref = MPLP_REF_INIT;
    bam_mplp_t iter;
    sam_hdr_t *h = NULL; /* header of first file in input list */
    char *ref;
    bool REGION = true;
    // faidx_t* fasta_idx = fai_load(params.refName.c_str());
    faidx_t* fasta_idx = fai_load3_format(params.refName.c_str(),
					  nullptr, nullptr, 0, FAI_FASTA);
    

    // FILE *pileup_fp = NULL;

    //bam_sample_t *sm = NULL;
    // kstring_t buf;
    //mplp_pileup_t gplp;

    //memset(&gplp, 0, sizeof(mplp_pileup_t));
    //memset(&buf, 0, sizeof(kstring_t));
    data = static_cast<mplp_aux_t**>(calloc(n_samples, sizeof(mplp_aux_t*)));
    plp = static_cast<const bam_pileup1_t**>(calloc(n_samples, sizeof(bam_pileup1_t*)));
    n_plp = static_cast<int*>(calloc(n_samples, sizeof(int)));
    //sm = bam_smpl_init();

    // initialize the bam file
    for(int i=0 ; i <n_samples; i++){
    sam_hdr_t *bam_hdr;
    data[i] = static_cast<mplp_aux_t*>(calloc(1, sizeof(mplp_aux_t)));
    data[i]->fp = sam_open(params.bamName.c_str(), "rb");
    if(! data[i]->fp)
      {
	std::cerr << "Could not open alignment file" << std::endl;
      }
    bam_hdr = sam_hdr_read(data[i]->fp);
    if(! bam_hdr){
      std::cerr << "Could not open alignment file header" << std::endl;
      exit(1);
    }

    
    hts_idx_t *idx = nullptr;
    // Loading index
    //I'm assuming the bam index is .bam.bai
    idx = sam_index_load(data[i]->fp, params.bamName.c_str());

    if (idx == nullptr){
      std::cerr << "Could Not Open Bam Index" << std::endl; 
    }

    if(REGION){
    // if ( (data[i]->iter=sam_itr_queryi(idx, 18, 1000, 6000)) == 0){
    // going to try not using an iterator, so It goes through all reads
    
    if ( (data[i]->iter=sam_itr_querys(idx, bam_hdr, params.regionString.c_str())) == 0){
      std::cerr << "Could not parse region" << std::endl;
      exit(1);
    }
    // std::cout << "Iter Start " <<  data[i]->iter->beg << " End " << data[i]->iter->end << std::endl;

    if (i == 0) beg0 = data[i]->iter->beg, end0 = data[i]->iter->end, tid0 = data[i]->iter->tid;
    }
    else {
      data[i]->iter = NULL;
    }
    // for some reason the samtools version only stores first header. I'm going to try loading each header
    data[i]->header = bam_hdr;
}

    // In samtools they have a struct that stores the current TID so they can verify which contig
    // we're on. As of now I only plan on calling mitochondrial pileups, but calling NUMT would change this
    int mtTid = sam_hdr_name2tid(data[0]->header, params.regionString.c_str());
    //int mtTid = sam_hdr_name2tid(data[0]->header, "chrM");
    hts_pos_t mtLength = sam_hdr_tid2len(data[0]->header, mtTid);
    hts_pos_t returnContigLength = -1;
    // this assumes we're only using one contig. I'm retrieving the mito sequence tid0
    // params.regionString.c_str(),
    ref = faidx_fetch_seq64(fasta_idx,
			    sam_hdr_tid2name(data[0]->header, mtTid),
                                0,
                                mtLength,
			    & returnContigLength);
  //
    // ref = fai_parse_region(fasta_idx, mtTid, 0, mtLength, int flags)
    // init pileup
  iter = bam_mplp_init(n_samples, mplp_func, (void**)data);

  // commenting these out for debugging
  // bam_mplp_init_overlaps(iter);
  // bam_mplp_set_maxcnt(iter, max_depth);

  // int last_tid = -1;
  // hts_pos_t last_pos = -1;

    // iterate through pileup
  int ret = 0;
  ret=bam_mplp64_auto(iter, &tid, &pos, n_plp, plp);
  // while ( (ret=bam_mplp64_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
  while( ret > 0) {
    //std::cout << "Current return value " << ret << std::endl;
    std::cout << "Position " << pos << std::endl;
    // todo: check if in region
    if( (tid != mtTid) || (pos < 0) || (pos >= mtLength) ){
      // continue or break / exit? Samtools continues here
      std::cerr << "Out of mitochondrial region" << std::endl;
      exit(1);
    }
    // if (conf->reg && (pos < beg0 || pos >= end0)) continue; // out of the region requested

    // samtools tries to deal with "missing portions of previous tids", which I'm not sure what they mean.

    // Samtools prints out the position and reference base here
    // iterate through the samples
    for (int i = 0; i < n_samples; i++){
      // iterate through the pileups returned for a sample
      // each sample at each position has its own sample base pileup
      char refBase = ref[pos];
      SampleBasePileup samplePileup = SampleBasePileup(pos +1, params.regionString, refBase);
      for( int j = 0; j < n_plp[i]; j++){
	// std::cout << "Current # of pileups " << n_plp[i] << std::endl;

	const bam_pileup1_t *tempPileup = plp[i] + j;
        // int baseChar = tempPileup->qpos < tempPileup->b->core.l_qseq
        //             ? bam_get_qual(tempPileup->b)[tempPileup->qpos]
        //             : 0;
	// if (baseChar - 33 >= minBaseQ) {
	    // from here I will process the pileup of this base.
    	    // Samtools passes ref and reflenth here

	// std::cout << "Processing Pileup" << std::endl;
	if(process_pileup(tempPileup ,pos, & params, samplePileup) <  0) {
	      std::cerr << "Error processing pileup" << std::endl;
	    }
	    
      }
      // finished iterating through sample, so output
      //std::cout << samplePileup;
      outFile << samplePileup;

    }
    ret=bam_mplp64_auto(iter, &tid, &pos, n_plp, plp);
  }
  std::cout << "Final return value " << ret << std::endl;
  outFile.close();
  bam_mplp_destroy(iter);
  return 0;
}
