#include <kmer_util.hpp>
#include <iostream>
#include <algorithm>    // std::min
#include <unordered_map>
#include <tuple>
#include <Kmer.hpp>
#include <iostream>
#include <fstream>

/*
    Helper map used for getting complement value of given nucleotide base.
*/
KmerUtil::ComplementMap KmerUtil::complement_map = { 
    { 'A', 'T' },
    { 'C', 'G' },
    { 'T', 'A' },
    { 'G', 'C' }
}; 

/*
    Used for calculating max number of occurences of particular mutation (I,X,D) on current index
    in the resulting map structure used for counting mutations.

    map - map which stores A,C,G,T values as keys & map of keys I,X,D with counters as values indicating
          how many times mutation occured 
    returns - tuple value (action which happened on the given mutation (I,X,D), nucleotide base, number of occurences) 
*/
std::tuple<char,char, int> KmerUtil::get_max_mutation_on_index(std::unordered_map <char, std::unordered_map<char, int>> map) {
    
    int max = -1;
    char base;
    char action;
    for (auto base_it = map.begin(); base_it != map.end(); ++base_it) {
        for(auto action_it = base_it->second.begin(); action_it != base_it->second.end(); ++action_it) {
            if(action_it->second > max) {
                max = action_it->second;
                action = action_it->first;
                base = base_it->first;
            }
        }
    }
    return std::make_tuple(action, base, max);
}


void KmerUtil::print_results_to_csv(std::string ref_string,  std::unordered_map<int, std::unordered_map<char, std::unordered_map<char,int>>> mutations_map, std::ofstream csv_out) {
    //std::ofstream csv_out;

    for (int i = 0; i < ref_string.size(); i++) {
        auto it =  mutations_map.find(i);
        if(it != mutations_map.end()) {
           // std::pair<char,char> action_mutation = get_max_mutation_on_index(it->second);
           // csv_out << action_mutation.first << "," << i << "," << action_mutation.second << ","  << std::endl;

        }
    }
    csv_out.close(); 
}

/*
    Adds mutation tuple to the result map.
    mutations_from_align - mutations found in align algorithm
    result_map - map in which mutation is added
*/
void KmerUtil::add_mutation_to_result_map(std::vector < std::tuple<char, int, char> > &mutations_from_align, 
    std::unordered_map<int, std::unordered_map<char, std::unordered_map<char,int>>> &result_map) {
    for (auto mutation : mutations_from_align) {
        char I_X_D = std::get<0>(mutation);
        int mutation_index_in_ref = std::get<1>(mutation);
        char base = std::get<2>(mutation);
        result_map[mutation_index_in_ref][base][I_X_D]++;
    }
}

/*
    Calculates reverse complement of input string.
    input - input string
*/
std::string KmerUtil::to_reverse_complement(std::string input) {
    std::string reverse_comp;
    for (int i = input.size() - 1; i >= 0; i--) {
        reverse_comp += KmerUtil::complement_map[input[i]];
    }
    return reverse_comp;
}

/*
    Calculates minimizers in input string of length k.
    input - input string
    k - kmer length
    returns - found kmer minimizer
*/
Kmer KmerUtil::minimizer_in_window(std::string &input, int k) {
    std::string max_str ("ZZZZZZZZZZZZZZZZ");
    Kmer minKmer(max_str, -1);
  //  std::cout << "trazim u: " << input << std::endl;
    for (int i=0; i < input.length() - k + 1; i++) {

        Kmer kmer(input.substr(i, k), i);

        if (kmer.str < minKmer.str) {
            minKmer.str = kmer.str;
            minKmer.index = i;
        }
    }
    return minKmer;
}

/*
    Calculates minimizers in given reference string using k and w.
    k - kmer length
    w - number of kmers from which the minimizer is generated
    returns - map of key minimizer, and value indicating vector of all indices found in reference string
*/
std::unordered_map<std::string, std::vector<int>> KmerUtil::calculate_minimizers(std::string &reference, int k, int w) {

    std::unordered_map<std::string, std::vector<int>> kmer_map;

    int win_size = w + k - 1;

    for (int i = 0; i <= reference.length() - win_size; i++) {
            std::string sub = reference.substr(i, win_size);
            Kmer minimizer = minimizer_in_window(sub, k);
            //std::cout << "  nasao: " << minimizer.str << " i: " << minimizer.index<<  std::endl; 
            auto iter = kmer_map.find(minimizer.str);
            //not found
            if(iter == kmer_map.end()) {
               // std::cout << "      nema: " << minimizer.str << std::endl;
                std::vector<int> kmer_indices = {i + minimizer.index};
                kmer_map.insert(std::pair<std::string, std::vector<int>>(minimizer.str, kmer_indices));
            } else {
                if (iter->second.back() != (i + minimizer.index)){
                    iter->second.push_back(i + minimizer.index);
                }
            }
    }

   return kmer_map;
}

/*
    Find longest incresing subsequence in vector of pairs based on the second value (index in ref string) of that pair.
    Stores the resulting indices in vector b.

    a - vector of pairs from which the LIS is calculated on second's pair value
    b - vector indicating longest increasing subsequence of ref indexes
    returns - true if LIS is found, false if not
*/
 bool KmerUtil::LongestIncreasingSubsequence(std::vector<std::pair<int, int>> &a, std::vector<int> &b) {
	std::vector<int> p(a.size());
	int u, v;
 
	if (a.empty()) return false;
 
	b.push_back(0);
 
	for (size_t i = 1; i < a.size(); i++){
		if (a[b.back()].second < a[i].second) {
			p[i] = b.back();
			b.push_back(i);
			continue;
		}
   
		for (u = 0, v = b.size()-1; u < v;) {
			int c = (u + v) / 2;
			if (a[b[c]].second < a[i].second) u=c+1; else v=c;
		}
 
		if (a[i].second < a[b[u]].second) {
			if (u > 0) p[i] = b[u-1];
			b[u] = i;
		}	
	}
 
	for (u = b.size(), v = b.back(); u--; v = p[v]) b[u] = v;
    return true;
} 

/*
    Helper function which gets index of given nucleotide char c
    used in scoring matrix W.
    c - nucleotide for which index is calculated
    returns - index of char c
*/
int mapChar(char c){
    if (c == 'A') return 0;
    else if (c == 'C') return 1;
    else if (c == 'G') return 2;
    else if (c == 'T') return 3;
}
/*
    Scoring matrix used for alignment.

    first - row in matrix
    second - col in matrix
    returns - score
*/
int W(char first, char second){
    int w[4][4] = {
        {4, -1, -1, -1},
        {-1, 4, -1, -1},
        {-1, -1, 4, -1},
        {-1, -1, -1, 4}
    };
    return w[mapChar(first)][mapChar(second)];
}

/*
    Global alignment function which aligns given strings s & t and returns number of found
    mutations.

    s - substring in ref string
    t - substring in sequenced string
    refIndex - index from which s substring starts in reference string
    returns - vector of found mutations in tuple format (action, index_in_ref, base)
*/
std::vector < std::tuple<char, int, char> > KmerUtil::globalAlignment(std::string &s, std::string &t, int refIndex){

    std::vector<std::vector<int>> V(s.length()+1, std::vector<int>(t.length()+1, 0));
    int d = -2;
    
    for (int i = 1; i <= s.length(); i++){
            V[i][0] = d*i;
    }
    for (int j = 1; j <= t.length(); j++){
        V[0][j] = d*j;
        //V[0][j] = 0;
    }

    for (int i = 1; i <= s.length();i++){
        #pragma omp parallel for
        for(int j = 1; j <= t.length(); j++){
                int MATCH = V[i-1][j-1] + W(s[i-1], t[j-1]);
                int INSERTION = V[i][j-1] + d;
                int DELETION = V[i-1][j] + d;
                V[i][j] = std::max(std::max(MATCH, INSERTION), DELETION);
                
        }
    }

    // print matrix
    /* for ( int i = 0; i <= s.length(); i ++){
        for( int j= 0; j <= t.length(); j++){
            std::cout << V[i][j] << " ";
        }
        std::cout << std::endl;
    } */
    
    std::string alignmentRef = "";
    std::string alignmentSeq = "";

    int i = s.length();
    int j = t.length();
        
    while (i>0 || j>0){
        if (i > 0 && j > 0 && V[i][j] == V[i-1][j-1] + W(s[i-1], t[j-1])){
            alignmentRef = s[i-1] + alignmentRef;
            alignmentSeq = t[j-1] + alignmentSeq;
            i -= 1;
            j -= 1;
        }
        else if (i > 0 && V[i][j] == V[i-1][j] + d){
            alignmentRef = s[i-1] + alignmentRef;
            alignmentSeq = "-" + alignmentSeq;
            i -=1;
        }
        else{
            alignmentRef = "-" + alignmentRef;
            alignmentSeq = t[j-1] + alignmentSeq;
            j -=1;
        } 
    }

    std::vector < std::tuple<char, int, char> > mutations;
    int realRefIndex = 0;

    for (int i = 0; i < alignmentRef.length(); i++){
            if (alignmentRef[i] == '-'){                    //insert
                mutations.push_back(std::make_tuple('I', refIndex+realRefIndex, alignmentSeq[i]));
            } 
            else if (alignmentSeq[i] == '-'){               //delete
                mutations.push_back(std::make_tuple('D', refIndex+realRefIndex, '-'));
                realRefIndex++;
            }
            else if (alignmentRef[i] != alignmentSeq[i]){   //substitucija
                mutations.push_back(std::make_tuple('X', refIndex+realRefIndex, alignmentSeq[i]));
                realRefIndex++;
            } 
            else {
                realRefIndex++;
            }  
    }

    //std::cout << alignmentRef << std::endl;
    //std::cout << alignmentSeq << std::endl; 
    //for (int i = 0; i < mutations.size(); i++){
    //    std::cout << std::get<0>(mutations[i]) << " " << std::get<1>(mutations[i]) << " "<< std::get<2>(mutations[i]) << std::endl;
    //}
 
    return mutations;
}

/*
    Calculates minimizers of sequence string while searching them in the reference string. 
    Returns longest region found in the tuple format (ref_substr, seq_substr, begin_ref, end_ref), so the
    sequence string can later be mapped to the reference in order to find mutations.

    refIndexer - used for indexing reference string kmer minimizers
    sequence - sequence string
    reference - reference string
    k - kmer length
    w - number of kmers from which the minimizer is generated
    min_gap - minimal allowed gap between minimizers in found region
    returns - ongest region found in the tuple format (ref_substr, seq_substr, begin_ref, end_ref)
*/
std::tuple<std::string, std::string, int, int> KmerUtil::find_best_region(KmerIndexer &refIndexer, std::string &sequence, std::string &reference, int k, int w, int min_gap) {

    std::unordered_map<std::string, std::vector<int>> kmer_map;
    std::vector <std::pair<int, int>> seqi_refi_pairs;
    int win_size = w + k - 1;

    for (int i = 0; i <= sequence.length() - win_size; i++) {
        std::string sub = sequence.substr(i, win_size);
        Kmer minimizer = minimizer_in_window(sub, k);
        int minimizer_index_in_seq = i + minimizer.index;
        std::vector<int> ref_indices =  refIndexer.get_kmer_indices(minimizer.str);
              
        if (seqi_refi_pairs.empty() && !ref_indices.empty()) {
            for (auto ref_index: ref_indices) {    
                seqi_refi_pairs.push_back(std::pair<int, int>(minimizer_index_in_seq, ref_index));
            }  
        } else if(!ref_indices.empty()){ 
            int last_inserted_index = seqi_refi_pairs.back().first;
            if (last_inserted_index != minimizer_index_in_seq) {                  
                  for (auto ref_index: ref_indices) {          
                    seqi_refi_pairs.push_back(std::pair<int, int>(minimizer_index_in_seq, ref_index));
                }     
            }
        }
    }

    /* std::cout << "All pairs: " << std::endl;
    for (auto el: seqi_refi_pairs) {
        std::cout << "(" << el.first << ", " << el.second << ") ";
    }
    std::cout << std::endl; */

    std::vector<int> longest_region_indices;
    bool success = LongestIncreasingSubsequence(seqi_refi_pairs, longest_region_indices);

    /* filter duplicated seq_index pairs */
    std::vector<std::pair<int,int>> lsi_filtered;
    int init = 1;
    int last_pushed;
    for (int i: longest_region_indices) {
        if (init || last_pushed != seqi_refi_pairs[i].first){
            init = 0;
            lsi_filtered.push_back(seqi_refi_pairs[i]);
            last_pushed = seqi_refi_pairs[i].first;
        }
    }
    seqi_refi_pairs = lsi_filtered;
    
    /* check if there is too big gap */
  
    int before_curr = -1;
    for (auto ind: longest_region_indices){
        if (before_curr < 0) {
            before_curr = seqi_refi_pairs[ind].second;
            continue;
        }
        int curr = seqi_refi_pairs[ind].second;

        if ( curr - before_curr > min_gap*(w+k-1) ){
            std::cout << "Too big gap, skipping...: " << curr - before_curr << std::endl;
            success = false;
            break;
        }
        before_curr = curr;
    }

    if (!success){
        return std::make_tuple("x","x",-1, -1);
    }

    /* std::cout << "LIS indices: " << std::endl;
    for (auto ind: longest_region_indices){
        std::cout << "( "  << seqi_refi_pairs[ind].first << ", " << seqi_refi_pairs[ind].second << " )"  << " ";
    }
    std::cout <<std::endl;
    */

    int begin_ref = seqi_refi_pairs[longest_region_indices.front()].second;
    int end_ref = seqi_refi_pairs[longest_region_indices.back()].second;
    int begin_seq = seqi_refi_pairs[longest_region_indices.front()].first;
    int end_seq = seqi_refi_pairs[longest_region_indices.back()].first;

    int end_str_ref_index = end_ref + k - 1;
    int end_str_seq_index = end_seq + k - 1;
    std::cout << "Found regions of best match:" << std::endl;
    std::cout << "Refrence str: ( begin_index: "  << begin_ref << ", end_index_start:" << end_str_ref_index << " )"  << std::endl;
    std::cout << "Sequence str: ( begin_index: "  << begin_seq << ", end_index_start:" << end_str_seq_index << " )"  << std::endl;
   
    std::string ref_substr = reference.substr(begin_ref, end_str_ref_index - begin_ref + 1);
    std::string seq_substr = sequence.substr(begin_seq, end_str_seq_index - begin_seq + 1);
     
    std::tuple<std::string, std::string, int, int> result (ref_substr, seq_substr, begin_ref, end_str_ref_index);

    return result; 
 
}