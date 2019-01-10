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

    int k = 5;
    int w = 5;

    //KmerIndexer KmerIndex(reference, w, k);

    std::string test_ref("23103210123356764738274123314");
    std::string test_seq("5841233101");
    KmerIndexer KmerIndex(test_ref, w, k);


    std::pair<std::string, std::string> reference_sequence = KmerUtil::find_best_region(KmerIndex, test_seq,test_ref,  w, k);
    std::string ref_substr = reference_sequence.first;
    std::string seq_substr = reference_sequence.second;

    std::cout << "Found reference substr: " << std::endl << ref_substr << std::endl;
    std::cout << "Found sequence substr: " << std::endl << seq_substr << std::endl;
 

    std::string s = "TCCAGAT";
    std::string t = "ACTCCTGAT";
    KmerUtil::globalAlignment(s,t, 1);


    return 0;

}