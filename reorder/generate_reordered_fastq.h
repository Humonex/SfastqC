/*
* Copyright 2018 University of Illinois Board of Trustees and Stanford
University. All Rights Reserved.
*We have used the reordering part of spring and improved them. Please refer
*to the paper to get the specific improvement methods.
* We promise that this algorithm is not used for any commercial purpose.
* And there is no exclusive behavior
*/

#ifndef SPRING_GENERATE_REORDERED_FASTQ_H_
#define SPRING_GENERATE_REORDERED_FASTQ_H_

#include <string>
#include <vector>

#include "util.h"

namespace spring {

void generate_reordered_fastq(const std::string &temp_dir,
                              const compression_params &cp,
                              const std::vector<std::string> &infile_vector,
                              const std::vector<std::string> &outfile_vector,
                              const bool gzipped_output_flag,
                              const bool gzipped_input_flag,
                              const bool fasta_flag);



void generate_order_pe(const std::string &file_order, uint32_t *order_array,
                       const uint32_t &numreads);

void generate_order_se(const std::string &file_order, uint32_t *order_array,
                       const uint32_t &numreads);

} // namespace spring

# endif // SPRING_GENERATE_REORDERED_FASTQ_H_
