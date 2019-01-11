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

    int k;
    int w;

    if (argc >= 2) {
        std::size_t pos;
        std::string arg = argv[1];
        std::string arg2 = argv[2];
        k = std::stoi(arg, &pos);
        w = std::stoi(arg2, &pos);
    } else {
        k = 15;
        w = 15;
    }

    std::cout << "k: " << k << " w: " << w << std::endl;
    //std::string reference_file ("./input/lambda.fasta");
    std::string sequence_file ("./input/lambda_simulated_reads.fasta");
    std::string reference_file ("./input/lambda.fasta");
    
    std::string reference = FastAReader::read_reference_file(reference_file);
    std::vector<std::string> sequence_list = FastAReader::read_sequence_file(sequence_file);

  


    


    /* std::ofstream csv_reg_test;    
    csv_reg_test.open ("reg.csv"); */


    std::unordered_map<int, std::unordered_map<char, std::unordered_map<char,int>>> ref_mutations;
    int seq_counter = 0;
    KmerIndexer KmerIndex(reference, w, k);
    for (auto seq_str: sequence_list) {
        seq_counter++;
    
        std::cout << "Seq num: " << seq_counter << std::endl;
        std::string seq_str_reversed = KmerUtil::to_reverse_complement(seq_str);
        std::tuple<std::string, std::string, int, int> ref_seq_refstart_refend = KmerUtil::find_best_region(KmerIndex, seq_str, reference, w, k);
        std::tuple<std::string, std::string, int, int> ref_seq_refstart_refend_comp = KmerUtil::find_best_region(KmerIndex, seq_str_reversed, reference, w, k);
        

        int begin_ref = std::get<2>(ref_seq_refstart_refend);
        int end_ref = std::get<3>(ref_seq_refstart_refend);
        int begin_ref_comp = std::get<2>(ref_seq_refstart_refend_comp);
        int end_ref_comp = std::get<3>(ref_seq_refstart_refend_comp);

        int region_size = end_ref - begin_ref;
        int region_size_comp = end_ref_comp - begin_ref_comp;

        std::string ref_substr;
        std::string seq_substr;


        int ref_index_start;
        int ref_index_end;
        bool was_comp = false;

        std::cout << "reg_size: " << region_size << "reg_size_comp: " << region_size_comp << std::endl;
        if (region_size > region_size_comp) {
            ref_substr = std::get<0>(ref_seq_refstart_refend);     // podstring pronadjen u glavnom reference genomu
            seq_substr = std::get<1>(ref_seq_refstart_refend);     // podstring pronadjen u sequenc genomu
            ref_index_start = begin_ref;   
            ref_index_end = end_ref;
        } else if (region_size < region_size_comp){
            ref_substr = std::get<0>(ref_seq_refstart_refend_comp);     // podstring pronadjen u glavnom reference genomu
            seq_substr = std::get<1>(ref_seq_refstart_refend_comp);     // podstring pronadjen u sequenc genomu
            ref_index_start = begin_ref_comp;   
            ref_index_end = end_ref_comp;
            std::cout << "using complementary strand " << std::endl;
            was_comp = true;
        } else {
            std::cout<< "No equal minimizer." << std::endl;
            continue;
        }

        // if there is no common minimizer, continue
        if (ref_index_start < 0) continue;
        //csv_reg_test << "Sequence_nr: " << seq_counter++  << ", ref_start: " << ref_index_start  << ", ref_end: " << ref_index_end << ", complement: " << was_comp << std::endl;



        std::vector < std::tuple<char, int, char> > res = KmerUtil::globalAlignment(ref_substr, seq_substr, ref_index_start);        
        

        //should be KmerUtil::add_mutation_to_result (adding to result)
        for (auto mutation : res) {
            char I_X_D = std::get<0>(mutation);
            int mutation_index_in_ref = std::get<1>(mutation);
            char base = std::get<2>(mutation);
            ref_mutations[mutation_index_in_ref][base][I_X_D]++;
        }

        //if (seq_counter >= 10) break;
    
        

    }

    //csv_reg_test.close(); 


    // should be KmerUtil::print_results_to_csv
    std::ofstream csv_out;
    csv_out.open ("outt_test.csv");

    for (int i = 0; i < reference.size(); i++) {
        auto it =  ref_mutations.find(i);
        if(it != ref_mutations.end()) {
            std::pair<char,char> action_mutation = KmerUtil::get_max_mutation_on_index(it->second);
            csv_out << action_mutation.first << "," << i << "," << action_mutation.second << std::endl;
        }
    }
    csv_out.close(); 
    
    
    //KmerUtil::print_results_to_csv(reference, ref_mutations); 











    /* std::string  indexes("012345678910");
    std::string test_ref("CTAGACGTAGCGACTGACTGACTACGATACTAGACGTAGCGACTGACTGACTACGATACTAGACGTAGCGACTGACTGACTACGATACTAGACGTAGCGACTGACTGACTACGATACTAGACGTAGCGACTGACTGACTACGATACTAGACGTAGCGACTGACTGACTACGATACTAGACGTAGCGACTGACTGACTACGATACTAGACGTAGCGACTGACTGACTACGATACTAGACGTAGCGACTGACTGACTACGATACTAGACGTAGCGACTGACTGACTACGATACTAGACGTAGCGACTGACTGACTACGATA");
    std::string test_seq("CTAGCCTGACGATCATGACGTACGCTAGCCTGACGATCATGACGTACGCTAGCCTGACGATCATGACGTACG");
    KmerIndexer KmerIndex(test_ref, w, k);
    
    std::cout << "ind: " << indexes << std::endl;
    std::cout << "ref: " << test_ref << std::endl;
    std::cout << "seq: " << test_seq << std::endl;
    std::tuple<std::string, std::string, int, int> reference_sequence_index = KmerUtil::find_best_region(KmerIndex, test_seq, test_ref,  w, k);

    std::string ref_substr = std::get<0>(reference_sequence_index);     // podstring pronadjen u glavnom reference genomu
    std::string seq_substr = std::get<1>(reference_sequence_index);     // podstring pronadjen u sequenc genomu
    int ref_index = std::get<2>(reference_sequence_index);              // index od kojeg počinje ref_substr u glavnom reference genomu
    int ref_size = std::get<3>(reference_sequence_index);
    
    std::cout << "Found reference substr: " << std::endl << ref_substr << std::endl;
    std::cout << "Found sequence substr: " << std::endl << seq_substr << std::endl;
    std::cout << "Ref string size: " << std::endl << ref_size << std::endl;

 

    

    std::cout << "ref_index: " << ref_index << std::endl;
    std::vector < std::tuple<char, int, char> > outVector = KmerUtil::globalAlignment(ref_substr, seq_substr, ref_index);


    std::unordered_map<int, std::unordered_map<char, std::unordered_map<char,int>>> ref_mutations_test;


    for (auto mutation : outVector) {
        
        char I_X_D = std::get<0>(mutation);
        int mutation_index_in_ref = std::get<1>(mutation);
        char base = std::get<2>(mutation);


        ref_mutations_test[mutation_index_in_ref][base][I_X_D]++;

        // csv_out << std::get<0>(mutation) << "," << std::get<1>(mutation) << "," << std::get<2>(mutation) << std::endl;
    }



    std::ofstream csv_out;
    csv_out.open ("out_test.csv");

    for (int i = 0; i < reference.size(); i++) {
        auto it =  ref_mutations_test.find(i);
        if(it != ref_mutations_test.end()) {
            std::pair<char,char> action_mutation = KmerUtil::get_max_mutation_on_index(it->second);
            csv_out << action_mutation.first << "," << i << "," << action_mutation.second << std::endl;
        }
    }
  

    csv_out.close();  */
     
    return 0;

}