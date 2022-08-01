/*
* Copyright 2018 University of Illinois Board of Trustees and Stanford
University. All Rights Reserved.
*We have used the reordering part of spring and improved them. Please refer
*to the paper to get the specific improvement methods.
* We promise that this algorithm is not used for any commercial purpose.
* And there is no exclusive behavior
*/

#ifndef SPRING_PE_ENCODE_H_
#define SPRING_PE_ENCODE_H_

#include <string>
#include "util.h"

namespace spring {

void pe_encode(const std::string &temp_dir, const compression_params &cp);

}  // namespace spring

#endif  // SPRING_PE_ENCODE_H_
