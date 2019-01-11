#include <string>
#include <vector>
#include <unordered_map>
#include <Kmer.hpp>
#include <KmerIndexer.hpp>

class KmerUtil {
public:

  typedef std::unordered_map<char, char> ComplementMap;
  static ComplementMap complement_map;

  typedef std::unordered_map<char, int> BaseToTupleIndex;
  static BaseToTupleIndex base_to_index;

  typedef std::unordered_map<int, char> IndexToBase;
  static IndexToBase index_to_base;

  static  std::string to_reverse_complement(std::string input);
  static  void LongestIncreasingSubsequence(std::vector<int> &a, std::vector<int> &b);
  static  bool LongestIncreasingSubsequence(std::vector <std::pair<int, int>> &a, std::vector<int> &b);
  static  std::vector < std::tuple<char, int, char> > globalAlignment(std::string &s, std::string &t, int refIndex);
  static  std::vector < std::tuple<char, int, char> > globalAlignment2(std::string &s, std::string &t, int refIndex);
  static  std::unordered_map<std::string, std::vector<int>> calculate_minimizers(std::string &reference, int k, int w);
  static  Kmer minimizer_in_window(std::string &reference, int k);
  static  std::tuple<std::string, std::string, int, int> find_best_region(KmerIndexer &refIndexer, std::string &sequence, std::string &reference, int w, int k);

  static std::pair <char, char> get_max_mutation_on_index(std::unordered_map <char, std::unordered_map<char, int>> map);
  static void print_results_to_csv(std::string ref_string,  std::unordered_map<int, std::unordered_map<char, std::unordered_map<char,int>>> mutations_map, std::ofstream csv_out);
  static void add_mutation_to_result_map(std::vector < std::tuple<char, int, char> > mutations_from_align, std::unordered_map<int, std::unordered_map<char, std::unordered_map<char,int>>> result_map);
};



