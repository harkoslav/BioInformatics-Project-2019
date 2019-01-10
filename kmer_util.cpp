#include <kmer_util.hpp>
#include <iostream>
#include <algorithm>    // std::min
#include <unordered_map>
#include <tuple>
#include <Kmer.hpp>
#include <iostream>


Kmer KmerUtil::minimizer_in_window(std::string &input, int k) {
    std::string max_str ("ZZZZ");
    Kmer minKmer(max_str, -1);
  //  std::cout << "trazim u: " << input << std::endl;
    for (int i=0,kmer_count=1; i < input.length() - k + 1; i++) {

        Kmer kmer(input.substr(i, k), i);

        if (kmer.str < minKmer.str) {
            minKmer.str = kmer.str;
            minKmer.index = i;
        }
    }
    return minKmer;
}

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


void KmerUtil::LongestIncreasingSubsequence(std::vector<int> &a, std::vector<int> &b)
{
	std::vector<int> p(a.size());
	int u, v;
 
	if (a.empty()) return;
 
	b.push_back(0);
 
	for (size_t i = 1; i < a.size(); i++){
		if (a[b.back()] < a[i]) {
			p[i] = b.back();
			b.push_back(i);
			continue;
		}
   
		for (u = 0, v = b.size()-1; u < v;) {
			int c = (u + v) / 2;
			if (a[b[c]] < a[i]) u=c+1; else v=c;
		}
 
		if (a[i] < a[b[u]]) {
			if (u > 0) p[i] = b[u-1];
			b[u] = i;
		}	
	}
 
	for (u = b.size(), v = b.back(); u--; v = p[v]) b[u] = v;
}


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

int mapChar(char c){
    if (c == 'A') return 0;
    else if (c == 'C') return 1;
    else if (c == 'G') return 2;
    else if (c == 'T') return 3;
}
int W(char first, char second){
    int w[4][4] = {
        {4, -1, -1, -1},
        {-1, 4, -1, -1},
        {-1, -1, 4, -1},
        {-1, -1, -1, 4}
    };
    return w[mapChar(first)][mapChar(second)];
}

std::vector < std::tuple<char, int, char> > KmerUtil::globalAlignment(std::string &s, std::string &t, int refIndex){
    int V[s.length()+1][t.length()+1];
    V[0][0]= 0;
    int d = -2;
    
    for (int i = 1; i <= s.length(); i++){
            V[i][0] = d*i;
    }
    for (int j = 1; j <= t.length(); j++){
        V[0][j] = d*j;
    }
    for ( int i = 1; i <= s.length();i++){
        for( int j = 1; j <= t.length(); j++){
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
            else{
                realRefIndex++;
            }  
    }
    /* std::cout << alignmentRef << std::endl;
    std::cout << alignmentSeq << std::endl;
    for (int i = 0; i < mutations.size(); i++){
        std::cout << std::get<0>(mutations[i]) << " " << std::get<1>(mutations[i]) << " "<< std::get<2>(mutations[i]) << std::endl;
    } */

    return mutations;
    //std::cout << alignmentSeq;
}

std::tuple<std::string, std::string, int> KmerUtil::find_best_region(KmerIndexer &refIndexer, std::string &sequence, std::string &reference, int k, int w) {
//std::pair<std::string, std::string> KmerUtil::find_best_region(KmerIndexer &refIndexer, std::string &sequence, std::string &reference, int k, int w) {

    std::unordered_map<std::string, std::vector<int>> kmer_map;
    std::vector <std::pair<int, int>> seqi_refi_pairs;
    int win_size = w + k - 1;

    for (int i = 0; i <= sequence.length() - win_size; i++) {
        std::string sub = sequence.substr(i, win_size);
        Kmer minimizer = minimizer_in_window(sub, k);
        //std::cout << "  nasao: " << minimizer.str << " i: " << i+minimizer.index<<  std::endl; 
        int minimizer_index_in_seq = i + minimizer.index;
        std::vector<int> ref_indices =  refIndexer.get_kmer_indices(minimizer.str);

        if (seqi_refi_pairs.empty()) {
            for (auto ref_index: ref_indices) {          
                    seqi_refi_pairs.push_back(std::pair<int, int>(minimizer_index_in_seq, ref_index));
                }  
        } else { 
            //std::cout << seqi_refi_pairs.back().first << ": " << minimizer_index_in_seq << std::endl;
            int last_inserted_index = seqi_refi_pairs.back().first;
        
            if (last_inserted_index != minimizer_index_in_seq) {
                for (auto ref_index: ref_indices) {          
                    seqi_refi_pairs.push_back(std::pair<int, int>(minimizer_index_in_seq, ref_index));
                }      
            }
        }
    }

    //std::cout << "Sequence string: " <<std::endl << sequence << std::endl;
    /* for (auto el: seqi_refi_pairs) {
        std::cout << "(" << el.first << ", " << el.second << ") ";
    } */
    //std::cout << std::endl;
    std::vector<int> longest_region_indices;
    
    bool success = LongestIncreasingSubsequence(seqi_refi_pairs, longest_region_indices);


    if (!success){
        std::cout << "no common minimizers!!!!" <<std::endl;
        return std::make_tuple("x","x",-1);
    }

    

    //std::cout << "LSI indices: " << std::endl;
   /*  for (auto ind: longest_region_indices){
        //std::cout << "( "  << seqi_refi_pairs[ind].first << ", " << seqi_refi_pairs[ind].second << " )"  << " ";
    } */
    //std::cout <<std::endl;


    
  

    int begin_ref = seqi_refi_pairs[longest_region_indices.front()].second;
    int end_ref = seqi_refi_pairs[longest_region_indices.back()].second;
    int begin_seq = seqi_refi_pairs[longest_region_indices.front()].first;
    int end_seq = seqi_refi_pairs[longest_region_indices.back()].first;

    std::cout << "Found regions of best match:" << std::endl;
    std::cout << "Refrence str: ( begin_index: "  << begin_ref << ", end_index:" << end_ref + k << " )"  << std::endl;
    std::cout << "Sequence str: ( begin_index: "  << begin_seq << ", end_index:" << end_seq + k << " )"  << std::endl;

    std::string ref_substr = reference.substr(begin_ref, end_ref - begin_ref + k + 1);
    std::string seq_substr = sequence.substr(begin_seq, end_seq - begin_seq + k + 1);

    std::tuple<std::string, std::string, int> result (ref_substr, seq_substr, begin_ref);
    return result; 
    //return std::make_pair(ref_substr, seq_substr);
 
}