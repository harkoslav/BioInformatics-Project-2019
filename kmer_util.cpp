#include <kmer_util.hpp>


std::vector<std::string> KmerUtil::calculate_minimizers(std::string reference, int k, int w) {

    std::vector<std::string> kmers;
    std::vector<std::string> kmer_minimizers;
    std::string maxStr ("ZZZZ");
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