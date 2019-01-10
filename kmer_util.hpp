#include <string>
#include <vector>
#include <unordered_map>
#include <Kmer.hpp>
#include <KmerIndexer.hpp>

class KmerUtil {
public:
  static  void LongestIncreasingSubsequence(std::vector<int> &a, std::vector<int> &b);
  static  bool LongestIncreasingSubsequence(std::vector <std::pair<int, int>> &a, std::vector<int> &b);
  static  std::vector < std::tuple<char, int, char> > globalAlignment(std::string &s, std::string &t, int refIndex);
  static  std::unordered_map<std::string, std::vector<int>> calculate_minimizers(std::string &reference, int k, int w);
  static  Kmer minimizer_in_window(std::string &reference, int k);
  //static  std::pair<std::string, std::string> find_best_region(KmerIndexer &refIndexer, std::string &sequence, std::string &reference, int w, int k);
  static  std::tuple<std::string, std::string, int> find_best_region(KmerIndexer &refIndexer, std::string &sequence, std::string &reference, int w, int k);
};



