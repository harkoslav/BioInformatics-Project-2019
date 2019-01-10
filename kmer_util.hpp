#include <string>
#include <vector>
#include <unordered_map>
#include <Kmer.hpp>

class KmerUtil {
public:
  static  void LongestIncreasingSubsequence(std::vector<int> &a, std::vector<int> &b);
  static void globalAlignment(std::string &s, std::string &t);
  static  std::unordered_map<std::string, std::vector<int>> calculate_minimizers(std::string reference, int k, int w);
  static  Kmer minimizer_in_window(std::string reference, int k);

};



