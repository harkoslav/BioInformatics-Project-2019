#include <string>

#ifndef KMER
#define KMER

/*
Class representing k-mer substring and it's index
in the reference string.
*/
class Kmer {

    public:
        std::string str;
        int index;
        Kmer(std::string, int i);
    
};
#endif