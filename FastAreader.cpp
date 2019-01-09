#include <FastAReader.hpp>
#include <iostream>
#include <fstream>

std::string FastAReader::read_reference_file(std::string path) {

    std::string reference_string;
    std::ifstream file_stream(path);

    if(file_stream.is_open()){
        std::string line;
        for (int i = 0; getline(file_stream, line); i++) if (i != 0) reference_string += line;
        file_stream.close();
    }

    return reference_string;

}


std::vector<std::string> FastAReader::read_sequence_file(std::string path) {

    std::vector<std::string> sequence_list;
    std::ifstream file_stream(path);

    if (file_stream.is_open())
    {
        std::string line;
        for (int i=0; getline(file_stream, line); i++) {
                getline(file_stream, line);
                sequence_list.push_back(line);
        }
        file_stream.close();
    }
    return sequence_list;
}