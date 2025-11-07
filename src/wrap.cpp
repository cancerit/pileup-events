#include "count.hpp"

std::vector<int> exec_pev (std::string aln_path,
                           std::string region_str,
                           bool no_overlaps = false,
                           int min_mapq = 25,
                           int min_baseq = 30,
                           int include_flag = 0,
                           int exclude_flag = 3844,
                           int max_depth = 1000000,
                           int clip_bound = 0) {
    htsFile *aln_in;
    bam_hdr_t *head;
    hts_region reg;
    count_params cp{min_baseq, min_mapq,     clip_bound,
                    max_depth, include_flag, exclude_flag};

    hts_idx_t *idx;
    int tid = -3;
    int64_t start, end;
    std::vector<int> result;
    try {
        aln_in = hts_open (aln_path.c_str(), "r");
        head = sam_hdr_read (aln_in);
        if (head == NULL) {
            throw std::runtime_error (
                "failed to get header from alignment file");
        }

        auto rp =
            sam_parse_region (head, region_str.c_str(), &tid, &start,
                              &end, HTS_PARSE_ONE_COORD);
        if (rp == NULL) {
            std::string msg;
            switch (tid) {
                case -2:
                    msg = "memory error";
                    break;
                case -1:
                    msg = "could not parse contig";
                    break;
                default:
                    msg = "specified range could not be parsed";
            }
            throw std::runtime_error (
                "parse failed for input region " + region_str +
                " - " + msg);
        }
        // converts to 0-indexed internal postions from 1-indexed
        // region str
        reg = hts_region::by_end (tid, start, end);

        idx = sam_index_load (aln_in, aln_path.c_str());
        if (!idx) {
            throw std::runtime_error ("failed to load index file");
        }

        safe_size_opts sso;
        sso.msg =
            "error in calculating cells needed for storing result";
        size_t n_cells = safe_size (
            static_cast<int64_t> (reg.rlen * N_FIELDS_PER_OBS));
        result.resize (n_cells, 0);

    } catch (std::exception &e) {
        throw std::runtime_error ("Error during setup: " +
                                  std::to_string (*e.what()));
    }

    try {
        AlleleEventCounter aev (cp, result, AEVSettings{no_overlaps});
        count (aln_in, idx, aev, reg, cp);
    } catch (std::exception &e) {
        throw std::runtime_error ("Error during calculation: " +
                                  std::to_string (*e.what()));
    }

    return result;
}
