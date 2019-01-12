#include <KmerIndexer.hpp>
#include <kmer_util.hpp>
#include <iostream>

KmerIndexer::KmerIndexer(std::string &input, int w_init, int k_init) {
    w = w_init;
    k = k_init;
    reference_length = input.length();
    //std::cout << "Creating minimizers for reference string: " << std::endl << input << std::endl;
    generate_minimizers(input);
}



void KmerIndexer::generate_minimizers(std::string &input) {
   reference_map =  KmerUtil::calculate_minimizers(input, k, w);

   /*
    for(auto el: reference_map) {
        std::cout << el.first << ": ";

        for (auto i: el.second) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
    */
}

std::vector<int> KmerIndexer::get_kmer_indices(std::string &kmer){
    return reference_map[kmer];
}
