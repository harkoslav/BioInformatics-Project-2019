#include <string>
#include <vector>
#include <unordered_map>
#include <Kmer.hpp>
#include <KmerIndexer.hpp>

/*
 Util class containing all the relevant functions used in algorithm.
*/
class KmerUtil {

public:

  // map containing complement nucleotide values
  typedef std::unordered_map<char, char> ComplementMap;
  static ComplementMap complement_map;

  static  std::string to_reverse_complement(std::string input);
  static  bool LongestIncreasingSubsequence(std::vector <std::pair<int, int>> &a, std::vector<int> &b);
  static  Kmer minimizer_in_window(std::string &reference, int k);
  static  void print_results_to_csv(std::string ref_string,  std::unordered_map<int, std::unordered_map<char, std::unordered_map<char,int>>> mutations_map, std::ofstream csv_out);
  static  void add_mutation_to_result_map(std::vector < std::tuple<char, int, char> > &mutations_from_align, std::unordered_map<int, std::unordered_map<char, std::unordered_map<char,int>>> &result_map);
  static  std::vector < std::tuple<char, int, char> > globalAlignment(std::string &s, std::string &t, int refIndex);
  static  std::unordered_map<std::string, std::vector<int>> calculate_minimizers(std::string &reference, int k, int w);
  static  std::tuple<std::string, std::string, int, int> find_best_region(KmerIndexer &refIndexer, std::string &sequence, std::string &reference, int w, int k, int min_gap);
  static std::tuple <char, char, int> get_max_mutation_on_index(std::unordered_map <char, std::unordered_map<char, int>> map);
  
};



