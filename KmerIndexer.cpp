#include <KmerIndexer.hpp>
#include <kmer_util.hpp>
#include <iostream>

/*
    Constructor of KmerIndexer.
    input - string which should be indexed (reference)
    k_init - kmer length
    w_init - number of kmers from which the minimizer is generated
*/
KmerIndexer::KmerIndexer(std::string &input, int w_init, int k_init) {
    w = w_init;
    k = k_init;
    generate_minimizers(input);
}

/*
    Generates minimizer for given input string.
    input - input string
*/
void KmerIndexer::generate_minimizers(std::string &input) {
   reference_map =  KmerUtil::calculate_minimizers(input, k, w);
}

/*
    Generates indices where given kmer is found in the reference string (input).
    kmer - kmer for which indices should be found
    returns - vector of indices for given kmer
*/
std::vector<int> KmerIndexer::get_kmer_indices(std::string &kmer){
    return reference_map[kmer];
}
