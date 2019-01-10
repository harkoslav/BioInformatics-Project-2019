#include <FastAReader.hpp>
#include <kmer_util.hpp>
#include <KmerIndexer.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <numeric>
#include <ctype.h>
#include <algorithm>
#include <tuple>


int main(int argc, char *argv[])
{

    //std::string reference_file ("./input/lambda.fasta");
    std::string sequence_file ("./input/lambda_simulated_reads.fasta");
    std::string reference_file ("./test_min.fasta");
    
    std::string reference = FastAReader::read_reference_file(reference_file);
    std::vector<std::string> sequence_list = FastAReader::read_sequence_file(sequence_file);

    int k = 3;
    int w = 3;



    KmerIndexer KmerIndex(reference, w, k);

    std::string test_minimizer("101");
    std::vector<int> indices =  KmerIndex.get_kmer_indices(test_minimizer);
    std::cout << "Indices for miminizer 101: " << std::endl;
    
    for (auto el: indices) {
        std::cout << el << std::endl;
    }

    //std::cout << "Sequence file [last]: " << sequence_list.back() << std::endl;
    //print seq readings
    /* for(auto it = sequence_list.begin(); it != sequence_list.end(); ++it){
        std::cout << *it << std::endl;
    } */

   // std::string testt ("2310343");
   // std::string ret = KmerUtil::minimizer_in_window(testt, 3).str;

    

    /*
    Test LogensIncreasingAlgo
    */
    int a[] = { 1, 9, 3, 8, 11, 4, 5, 6, 4, 19, 7, 1, 7 };
    std::vector<int> seq(std::begin(a), std::end(a));
    std::vector<int> result_subseq_indices;                              // result_subseq_indices : Vector containing indexes of longest subsequence 
    
    
    KmerUtil::LongestIncreasingSubsequence(seq, result_subseq_indices);

    std::cout << "1, 9, 3, 8, 11, 4, 5, 6, 4, 19, 7, 1, 7" << std::endl;
    std::cout << "Longes icreasing subseq: ";    
	for (size_t i = 0; i < result_subseq_indices.size(); i++)
		std::cout << seq[result_subseq_indices[i]];
 
    std::cout << std::endl;


    std::string s = "TCCAGAT";
    std::string t = "ACTCCTGAT";
    KmerUtil::globalAlignment(s,t, 1);


    return 0;

}