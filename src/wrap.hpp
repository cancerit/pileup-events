#include <vector>

std::vector<int> exec_pev (std::string aln_path,
                           std::string region_str,
                           bool no_overlaps = false,
                           int min_mapq = 25,
                           int min_baseq = 30,
                           int include_flag = 0,
                           int exclude_flag = 3844,
                           int max_depth = 1000000,
                           int clip_bound = 0);
