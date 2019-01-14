#include <FastAReader.hpp>
#include <kmer_util.hpp>
#include <KmerIndexer.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <tuple>
#include <chrono>

using namespace std::chrono;


int main(int argc, char *argv[])
{
     /* --------------------------- INPUT --------------------------- */
    std::chrono::high_resolution_clock::time_point t1_start = high_resolution_clock::now();

    int k;
    int w;
    int c_treshold;
    int min_gap;
    int ecoli;
    if (argc >= 5) {
        std::size_t pos;
        std::string arg = argv[1];
        std::string arg2 = argv[2];
        std::string arg3 = argv[3];
        std::string arg4 = argv[4];
        std::string arg5 = argv[5];
        k = std::stoi(arg, &pos);
        w = std::stoi(arg2, &pos);
        c_treshold = std::stoi(arg3, &pos);
        min_gap = std::stoi(arg4, &pos);
        ecoli = std::stoi(arg5, &pos);
    } else {
        k = 15;
        w = 5;
        c_treshold = 12;
        min_gap = 15;
        ecoli = 0;
    }
    std::string file = (ecoli) ? "ecoli":"lambda";
    std::cout << "k: " << k << " w: " << w <<" c_tres: " << c_treshold << " gap: " << min_gap << " file: "<< file << std::endl;
    std::string sequence_file;
    std::string reference_file;

    if(ecoli == 1){
        sequence_file = "./input/ecoli_simulated_reads.fasta";
        reference_file = "./input/ecoli.fasta";
    } else {
       sequence_file = "./input/lambda_simulated_reads.fasta";
       reference_file = "./input/lambda.fasta"; 
    }
    
/*     std::string sequence_file ("./input/ecoli_simulated_reads.fasta");
    std::string reference_file ("./input/ecoli.fasta"); */
    
    std::string reference = FastAReader::read_reference_file(reference_file);
    std::vector<std::string> sequence_list = FastAReader::read_sequence_file(sequence_file);
    /* --------------------------- END OF INPUT --------------------------- */
    

    /* Time counters */
    int region_sum = 0;
    int align_sum = 0;
    int seq_counter = 0;

    /* Map used for storing number of mutation occurances for each index in reference string.*/
    std::unordered_map<int, std::unordered_map<char, std::unordered_map<char,int>>> ref_mutations;
    KmerIndexer KmerIndex(reference, w, k);

    /* ------------------------------------------------------ START OF MAIN LOOP ------------------------------------------------------ */
    for (auto seq_str: sequence_list) {
        seq_counter++;
    
        std::cout << "Seq num: " << seq_counter << std::endl;
        std::string seq_str_reversed = KmerUtil::to_reverse_complement(seq_str);
    
        /*Find best region in ref string using minimizers and get substrings of ref and seq representing similiar parts in which mutations should be found. */
        std::chrono::high_resolution_clock::time_point t1 = high_resolution_clock::now();
        std::tuple<std::string, std::string, int, int> ref_seq_refstart_refend = KmerUtil::find_best_region(KmerIndex, seq_str, reference, k, w, min_gap);
        std::tuple<std::string, std::string, int, int> ref_seq_refstart_refend_comp = KmerUtil::find_best_region(KmerIndex, seq_str_reversed, reference, k, w, min_gap);
        std::chrono::high_resolution_clock::time_point t2 = high_resolution_clock::now();

        auto duration = duration_cast<microseconds>( t2 - t1 ).count();
        region_sum += duration;

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

        /* Decide if complementary or regular sequenced substr had longer region in ref, take that for alignment.*/
        if (region_size > region_size_comp) {
            ref_substr = std::get<0>(ref_seq_refstart_refend);     // substr in ref genome 
            seq_substr = std::get<1>(ref_seq_refstart_refend);     // substr in seqc genome
            ref_index_start = begin_ref;   
            ref_index_end = end_ref;
        } else if (region_size < region_size_comp){
            ref_substr = std::get<0>(ref_seq_refstart_refend_comp);     
            seq_substr = std::get<1>(ref_seq_refstart_refend_comp);     
            ref_index_start = begin_ref_comp;   
            ref_index_end = end_ref_comp;
            std::cout << "using complementary strand " << std::endl;
            was_comp = true;
        } else {
            std::cout<< "No equal minimizer." << std::endl;
            continue;
        }

        /* Align two substring and calculate mutations.*/
        std::chrono::high_resolution_clock::time_point t1a = high_resolution_clock::now();
        std::vector < std::tuple<char, int, char> > res = KmerUtil::globalAlignment(ref_substr, seq_substr, ref_index_start);        
        std::chrono::high_resolution_clock::time_point t2a = high_resolution_clock::now();

        auto duration2 = duration_cast<microseconds>( t2a - t1a ).count();
        align_sum += duration2;

        KmerUtil::add_mutation_to_result_map(res, ref_mutations);
    }
    /* ------------------------------------------------------ END OF MAIN LOOP ------------------------------------------------------ */



    /* --------------------------- PRINTS MUTATIONS TO CSV FILE --------------------------- */
    std::ofstream csv_out;
    csv_out.open ("out.csv");

    for (int i = 0; i < reference.size(); i++) {
        auto it =  ref_mutations.find(i);
        if(it != ref_mutations.end()) {
            std::tuple<char,char, int> action_mutation = KmerUtil::get_max_mutation_on_index(it->second);
            int count_tres = std::get<2>(action_mutation);
            char action = std::get<0>(action_mutation);
            if( count_tres >= c_treshold){
                //csv_out << std::get<0>(action_mutation) << "," << i << "," << std::get<1>(action_mutation) << "," << std::get<2>(action_mutation) << std::endl;
                csv_out << std::get<0>(action_mutation) << "," << i << "," << std::get<1>(action_mutation) << std::endl;
            }
        }
    }
    csv_out.close(); 
    /* --------------------------- END OF PRINTS MUTATIONS TO CSV FILE --------------------------- */



    /* --------------------------- PRINT DURATION TIMES --------------------------- */
    std::chrono::high_resolution_clock::time_point t2_end = high_resolution_clock::now();
    auto duration3 = duration_cast<microseconds>( t2_end - t1_start ).count();
        
    std::cout << "Align duration avg: " << align_sum << std::endl;
    std::cout << "Region duration avg: " << region_sum << std::endl;
    std::cout << "Program duration: " << duration3 << std::endl;
    /* --------------------------- END OF PRINT DURATION TIMES --------------------------- */

    return 0;
}