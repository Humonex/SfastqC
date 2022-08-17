/*
* Copyright 2018 University of Illinois Board of Trustees and Stanford
University. All Rights Reserved.
*We have used the reordering part of spring and improved them. Please refer
*to the paper to get the specific improvement methods.
* We promise that this algorithm is not used for any commercial purpose.
* And there is no exclusive behavior
*/

#ifndef SPRING_SPRING_H_
#define SPRING_SPRING_H_

#include <string>
#include "util.h"

namespace spring {

void spring_reorder(const std::string &temp_dir,
              const std::vector<std::string> &infile_vec,
              const std::vector<std::string> &outfile_vec, const int &num_thr,
              const bool &gzipped_input_flag, const bool &gzipped_output_flag,
              const bool &fasta_flag);

std::string random_string(size_t length);

}  // namespace spring

#endif  // SPRING_SPRING_H_
