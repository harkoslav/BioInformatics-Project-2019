#include <unordered_map>
#include <vector>
#include <string>

#ifndef KMER_INDEXER
#define KMER_INDEXER
class KmerIndexer {

    std::unordered_map<std::string, std::vector<int>> reference_map;
    int w;
    int k;
    int reference_length;

    private:
        void generate_minimizers(std::string &input);

    public: 
        KmerIndexer(std::string &input, int w_init, int k_init);
        std::vector<int> get_kmer_indices(std::string &kmer);
        int get_reference_length();
        

};

#endif