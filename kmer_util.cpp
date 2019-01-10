#include <kmer_util.hpp>
#include <iostream>
#include <algorithm>    // std::min

std::vector<std::string> KmerUtil::calculate_minimizers(std::string reference, int k, int w) {

    std::vector<std::string> kmers;
    std::vector<std::string> kmer_minimizers;
    std::string maxStr ("ZZZZZ");
    std::string min = maxStr;

    for (int i=0,kmer_count=1; i<reference.length() - k + 1; i++){
        std::string kmer = reference.substr(i, k);

        if (kmer < min) {            
            min = kmer;
        }

        //nadi min u novom windowu, pushaj trenutni min
        if (kmer_count == w) {
            kmer_count = 1;
            kmer_minimizers.push_back(min);
            min = maxStr;
        } else {
            kmer_count++;
        }
        
        kmers.push_back(kmer);
    }
    if (min != maxStr){
        kmer_minimizers.push_back(min);
    }

    return kmer_minimizers;
}

void KmerUtil::LongestIncreasingSubsequence(std::vector<int> &a, std::vector<int> &b)
{
	std::vector<int> p(a.size());
	int u, v;
 
	if (a.empty()) return;
 
	b.push_back(0);
 
	for (size_t i = 1; i < a.size(); i++) 
        {
		if (a[b.back()] < a[i]) 
                {
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

void  KmerUtil::globalAlignment(std::string &s, std::string &t){
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
    


}