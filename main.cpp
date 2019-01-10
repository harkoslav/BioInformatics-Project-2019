#include <FastAReader.hpp>
#include <kmer_util.hpp>
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

    std::string reference_file ("./input/lambda.fasta");
    std::string sequence_file ("./input/lambda_simulated_reads.fasta");
    //std::string reference_file ("./input/test.fasta");
    
    std::string reference = FastAReader::read_reference_file(reference_file);
    std::vector<std::string> sequence_list = FastAReader::read_sequence_file(sequence_file);

    int k = 5;
    int w = 15;
    int window_length = w + k - 1;
    
    std::vector <std::string> kmer_minimizers = KmerUtil::calculate_minimizers(reference, k, w);


    //print all kmers
    /* for (auto it = kmers.begin() ; it != kmers.end(); ++it){
        std::cout << *it << std::endl;
    } */

    //print minimizers
    std::cout << "Minimizers: " << std::endl;
    for (auto it = kmer_minimizers.begin(); it != kmer_minimizers.end(); ++it){
        std::cout << *it << std::endl;
    }


    //std::cout << "Sequence file [last]: " << sequence_list.back() << std::endl;
    //print seq readings
    /* for(auto it = sequence_list.begin(); it != sequence_list.end(); ++it){
        std::cout << *it << std::endl;
    } */




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


    std::string x = "AATC";
    std::string y = "AAC";
    KmerUtil::globalAlignment(x,y);


    return 0;

}