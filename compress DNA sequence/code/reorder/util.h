/*
* Copyright 2018 University of Illinois Board of Trustees and Stanford
University. All Rights Reserved.
*We have used the reordering part of spring and improved them. Please refer
*to the paper to get the specific improvement methods.
* We promise that this algorithm is not used for any commercial purpose.
* And there is no exclusive behavior
*/

#ifndef SPRING_UTIL_H_
#define SPRING_UTIL_H_

#include <fstream>
#include <string>

namespace spring {

static const char chartorevchar[128] = {
    0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0,
    0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0,
    0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 'T',
    0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, 0, 0, 0, 0, 'A', 0, 0, 0,
    0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0,
    0, 0,   0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0};
struct compression_params {
  bool paired_end;
  bool long_flag;//设置一个长度的标识符
  uint32_t num_reads;
  uint32_t num_reads_clean[2];
  uint32_t max_readlen;
  int num_reads_per_block;
  int num_thr;
};

uint32_t read_fastq_block(std::istream *fin, std::string *id_array,
                          std::string *read_array, std::string *quality_array,
                          const uint32_t &num_reads, const bool &fasta_flag);

void write_fastq_block(std::ofstream &fout, std::string *id_array,
                       std::string *read_array, std::string *quality_array,
                       const uint32_t &num_reads,
                       const bool preserve_quality, const int &num_thr,
                       const bool &gzip_flag);

void reverse_complement(char *s, char *s1, const int readlen);

std::string reverse_complement(const std::string &s, const int readlen);

void remove_CR_from_end(std::string &str);

bool is_permutation(uint32_t *order_array, const uint32_t &numreads);

void write_dna_in_bits(const std::string &read, std::ofstream &fout);

void read_dna_from_bits(std::string &read, std::ifstream &fin);

void write_dnaN_in_bits(const std::string &read, std::ofstream &fout);

void read_dnaN_from_bits(std::string &read, std::ifstream &fin);

size_t get_directory_size(const std::string &temp_dir);

}  // namespace spring

#endif  // SPRING_UTIL_H_
