#include <string>

#ifndef KMER
#define KMER
class Kmer {

    public:
        std::string str;
        int index;
        Kmer(std::string, int i);
    
};

#endif