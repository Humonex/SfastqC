/*
* Copyright 2018 University of Illinois Board of Trustees and Stanford
University. All Rights Reserved.
*We have used the reordering part of spring and improved them. Please refer
*to the paper to get the specific improvement methods.
* We promise that this algorithm is not used for any commercial purpose.
* And there is no exclusive behavior
*/

#ifndef SPRING_PREPROCESS_H_
#define SPRING_PREPROCESS_H_

#include <string>
#include "util.h"

namespace spring {

void preprocess(const std::string &infile_1, const std::string &infile_2,
                const std::string &temp_dir, compression_params &cp,
                const bool &gzip_flag, const bool &fasta_flag);

}  // namespace spring

#endif  // SPRING_PREPROCESS_H_
