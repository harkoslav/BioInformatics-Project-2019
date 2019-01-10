#include <string>
#include <vector>

class KmerUtil {
public:
  static  void LongestIncreasingSubsequence(std::vector<int> &a, std::vector<int> &b);
  static  std::vector<std::string> calculate_minimizers(std::string reference, int k, int w);
  static void globalAlignment(std::string &s, std::string &t);
};



