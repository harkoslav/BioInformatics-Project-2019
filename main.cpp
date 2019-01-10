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

    int k = 4;
    int w = 4;



    //ovako bi trebao ici cijeli algoritam
    /* KmerIndexer KmerIndex(reference, w, k);
    for (auto seq_str: sequence_list) {
        std::tuple<std::string, std::string, int> reference_sequence_index = KmerUtil::find_best_region(KmerIndex, seq_str, reference,  w, k);
        
        std::string ref_substr = std::get<0>(reference_sequence_index);     // podstring pronadjen u glavnom reference genomu
        std::string seq_substr = std::get<1>(reference_sequence_index);     // podstring pronadjen u sequenc genomu
        int ref_index = std::get<2>(reference_sequence_index);   

        tuple<char, int, char> res = KmerUtil::globalAlignment(?,?, ref_index);

        //todo spremi res za evaluaciju
    } */

    std::string test_ref("23103210123356764738274123314");
    std::string test_seq("5841233101");
    KmerIndexer KmerIndex(test_ref, w, k);

    
    std::tuple<std::string, std::string, int> reference_sequence_index = KmerUtil::find_best_region(KmerIndex, test_seq, test_ref,  w, k);

    std::string ref_substr = std::get<0>(reference_sequence_index);     // podstring pronadjen u glavnom reference genomu
    std::string seq_substr = std::get<1>(reference_sequence_index);     // podstring pronadjen u sequenc genomu
    int ref_index = std::get<2>(reference_sequence_index);              // index od kojeg počinje ref_substr u glavnom reference genomu

    std::cout << "Found reference substr: " << std::endl << ref_substr << std::endl;
    std::cout << "Found sequence substr: " << std::endl << seq_substr << std::endl;
 

    std::string s = "TCCAGAT";
    std::string t = "ACTCCTGAT";
    KmerUtil::globalAlignment(s,t, 1);
    return 0;

}