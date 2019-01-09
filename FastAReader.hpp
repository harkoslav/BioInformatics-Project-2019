#include <string>
#include <vector>

class FastAReader {
    public:
        static std::string read_reference_file(std::string path);
        static std::vector<std::string> read_sequence_file(std::string path); 

};