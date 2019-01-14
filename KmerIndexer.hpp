#include <unordered_map>
#include <vector>
#include <string>

#ifndef KMER_INDEXER
#define KMER_INDEXER

/*
    Class used for indexing k-mer minimizers in reference string.
*/
class KmerIndexer {

    std::unordered_map<std::string, std::vector<int>> reference_map;
    int w;
    int k;

    public: 
        void generate_minimizers(std::string &input);
        KmerIndexer(std::string &input, int w_init, int k_init);
        std::vector<int> get_kmer_indices(std::string &kmer);        

};
#endif