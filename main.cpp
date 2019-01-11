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
    std::string reference_file ("./input/lambda.fasta");
    
    std::string reference = FastAReader::read_reference_file(reference_file);
    std::vector<std::string> sequence_list = FastAReader::read_sequence_file(sequence_file);

    int k = 15;
    int w = 15;
  
    std::ofstream csv_out;
    csv_out.open ("out.csv");

    KmerIndexer KmerIndex(reference, w, k);
    for (auto seq_str: sequence_list) {
        std::tuple<std::string, std::string, int> reference_sequence_index = KmerUtil::find_best_region(KmerIndex, seq_str, reference,  w, k);
        
        std::cout<< "Proslo " << std::endl;
        std::string ref_substr = std::get<0>(reference_sequence_index);     // podstring pronadjen u glavnom reference genomu
        std::string seq_substr = std::get<1>(reference_sequence_index);     // podstring pronadjen u sequenc genomu
        int ref_index = std::get<2>(reference_sequence_index);   


        // if there is no common minimizer, continue
        if (ref_index < 0) continue;

        std::vector < std::tuple<char, int, char> > res = KmerUtil::globalAlignment(ref_substr, seq_substr, ref_index);
        

        for (auto mutation : res) {

            //todo: remove in future
            if (std::get<2>(mutation) < 0) continue;
            
            csv_out << std::get<0>(mutation) << "," << std::get<1>(mutation) << "," << std::get<2>(mutation) << std::endl;
        }

    }

    csv_out.close();

    /* std::string test_ref("ACTGACGAGCAGTAGACGTACTGGGAA");
    std::string test_seq("GTACTACAG");
    KmerIndexer KmerIndex(test_ref, w, k);

    
    std::tuple<std::string, std::string, int> reference_sequence_index = KmerUtil::find_best_region(KmerIndex, test_seq, test_ref,  w, k);

    std::string ref_substr = std::get<0>(reference_sequence_index);     // podstring pronadjen u glavnom reference genomu
    std::string seq_substr = std::get<1>(reference_sequence_index);     // podstring pronadjen u sequenc genomu
    int ref_index = std::get<2>(reference_sequence_index);              // index od kojeg počinje ref_substr u glavnom reference genomu

    std::cout << "Found reference substr: " << std::endl << ref_substr << std::endl;
    std::cout << "Found sequence substr: " << std::endl << seq_substr << std::endl;
 

    std::string ref = "TCCAGAT";
    std::string seq = "CTGAT";
    std::vector < std::tuple<char, int, char> > outVector = KmerUtil::globalAlignment2(ref, seq, 1);


    std::ofstream csv_out;
    csv_out.open ("out.csv");

    for (auto mutation : outVector) {
        csv_out << std::get<0>(mutation) << "," << std::get<1>(mutation) << "," << std::get<2>(mutation) << std::endl;
    }
    csv_out.close(); */
    return 0;

}