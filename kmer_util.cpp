#include <kmer_util.hpp>
#include <iostream>
#include <algorithm>    // std::min
#include <unordered_map>
#include <tuple>
#include <Kmer.hpp>
#include <iostream>


Kmer KmerUtil::minimizer_in_window(std::string input, int k) {
    std::string max_str ("ZZZZ");
    Kmer minKmer(max_str, -1);
    std::cout << "trazim u: " << input << std::endl;
    for (int i=0,kmer_count=1; i < input.length() - k + 1; i++) {

        Kmer kmer(input.substr(i, k), i);

        if (kmer.str < minKmer.str) {
            minKmer.str = kmer.str;
            minKmer.index = i;
        }
    }
    return minKmer;
}

std::unordered_map<std::string, std::vector<int>> KmerUtil::calculate_minimizers(std::string reference, int k, int w) {

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

void  KmerUtil::globalAlignment(std::string &s, std::string &t, int refIndex){
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
    for ( int i = 0; i <= s.length(); i ++){
        for( int j= 0; j <= t.length(); j++){
            std::cout << V[i][j] << " ";
        }
        std::cout << std::endl;
    }
    
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

    std::vector < std::pair<char, int> > mutations;
    int realRefIndex = 0;

    for (int i = 0; i < alignmentRef.length(); i++){
            if (alignmentRef[i] == '-'){                    //insert
                mutations.push_back(std::make_pair('I', refIndex+realRefIndex));
            } 
            else if (alignmentSeq[i] == '-'){               //
                mutations.push_back(std::make_pair('D', refIndex+realRefIndex));
                realRefIndex++;
            }
            else if (alignmentRef[i] != alignmentSeq[i]){   //substitucija
                mutations.push_back(std::make_pair('X', refIndex+realRefIndex));
                realRefIndex++;
            } 
            else{
                realRefIndex++;
            }  
    }
    std::cout << alignmentRef << std::endl;
    std::cout << alignmentSeq << std::endl;
    for (int i = 0; i < mutations.size(); i++){
        std::cout << mutations[i].first << " " << mutations[i].second << std::endl;
    }
}