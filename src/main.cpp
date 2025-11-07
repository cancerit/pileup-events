#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cxxopts.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "const.hpp"
#include "count.hpp"

int main (int argc,
          char *argv[]) {
    namespace fs = std::filesystem;

    fs::path aln_path;
    std::string region_str;
    htsFile *aln_in;
    bam_hdr_t *head;
    hts_region reg;
    count_params cp;

    // defaults
    bool print_head = false;
    bool print_row = false;
    bool no_overlaps = false;
    cp.min_mapq = 25;
    cp.min_baseq = 30;
    cp.include_flag = 0;
    cp.exclude_flag = 3844;
    cp.max_depth = 1000000;
    cp.clip_bound = 0;

    try {
        cxxopts::Options options (
            "pileup-events",
            "-----pileup-events:-------------------------------------"
            "|"
            "\n"
            "                                                        "
            "|"
            "\n"
            "  Count alleles and alignment events per position for   "
            "|"
            "\n"
            "  a specified genomic region.                           "
            "|"
            "\n"
            "                                                        "
            "|"
            "\n"
            "   -----------                                          "
            "|"
            "\n"
            "                                                        "
            "|"
            "\n"
            "  Where reference names contain colons, surround in     "
            "|"
            "\n"
            "  curly braces like {HLA-DRB1*12:17}:<start>{-<end>}.   "
            "|"
            "\n"
            "                                                        "
            "|"
            "\n"
            "  chr1:100 is treated as the single base pair region    "
            "|"
            "\n"
            "  chr1:100-100. chr1:-100 is shorthand for chr1:1-100   "
            "|"
            "\n"
            "  and chr1:100- is ch1:100-<end>. All co-ordinates are  "
            "|"
            "\n"
            "  1-based, end-inclusive; i.e. as reported in a VCF.    "
            "|"
            "\n"
            "                                                        "
            "|"
            "\n"
            "  A result matrix with end-start rows and 22 columns    "
            "|"
            "\n"
            "  of event counters (see --head) is printed to stdout   "
            "|"
            "\n"
            "  as a csv. The first 12 columns represent events on    "
            "|"
            "\n"
            "  the forward strand, the next 12 the reverse.          "
            "|"
            "\n"
            "                                                        "
            "|"
            "\n"
            "--------------------------------------------------------"
            "|"
            "\n");

        // clang-format off
        options.add_options()
            ("aln", "", cxxopts::value<fs::path>())  // positional
            ("region", "", cxxopts::value<std::string>())

            // parameters
            ("b,baseq",
             "Minimum base quality to treat base as unambiguous. (default 30)",
             cxxopts::value<int>())
            ("m,mapq",
             "Minimum mapping quality to include read (default 25)",
             cxxopts::value<int>())
            ("c,clip",
             "Treat bases within <clip> bases of read edges as ambiguous. (default 0)",
             cxxopts::value<int>())
            ("i,include",
             "Include only reads with all bits set in sam flag. Provide flag as integer. (default 0)",
             cxxopts::value<int>())
            ("e,exclude",
             "Exclude reads with any bits set in sam flag. Provide flag as integer. (default 3844)",
             cxxopts::value<int>())
            ("d,depth",
             "Maximum read depth (default 1000000)",
             cxxopts::value<int>())

            ("head", "Print header")
            ("row", "Print genomic position index for each row")
            ("discard-overlaps", "Avoid double counting of bases from the same template")
            ("h,help", "Print usage")
            ("version", "Print program version");  // ideally this would report the version of htslib compiled against
        // clang-format on

        options.parse_positional ({"aln", "region"});
        options.positional_help ("<.BAM/.CRAM> chr:start-end");
        auto parsed_args = options.parse (argc, argv);

        if (parsed_args.count ("help")) {
            std::cout << options.help() << std::endl;
            return 0; // nothing given nothing done
        }

        if (parsed_args.count ("version")) {
            std::cout << VERSION << std::endl;
            return 0;
        }

        if ((!parsed_args.count ("aln")) ||
            (!parsed_args.count ("region"))) {
            std::cout << "incorrect usage: all postional arguments "
                         "required. Try --help"
                      << std::endl;
            return 1;
        }

        aln_path = parsed_args["aln"].as<fs::path>();
        region_str = parsed_args["region"].as<std::string>();

        if (region_str.empty())
            throw std::runtime_error (
                "region string appears to be empty");

        if (parsed_args.count ("baseq")) {
            cp.min_baseq = parsed_args["baseq"].as<int>();
        }
        if (parsed_args.count ("mapq")) {
            cp.min_mapq = parsed_args["mapq"].as<int>();
        }
        if (parsed_args.count ("clip")) {
            cp.clip_bound = parsed_args["clip"].as<int>();
        }
        if (parsed_args.count ("include")) {
            cp.include_flag = parsed_args["include"].as<int>();
        }
        if (parsed_args.count ("exclude")) {
            cp.exclude_flag = parsed_args["exclude"].as<int>();
        }
        if (parsed_args.count ("depth")) {
            cp.max_depth = parsed_args["depth"].as<int>();
        }
        if (parsed_args.count ("head")) {
            print_head = true;
        }
        if (parsed_args.count ("row")) {
            print_row = true;
        }
        if (parsed_args.count ("discard-overlaps")) {
            no_overlaps = true;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error parsing CLI options: " << e.what()
                  << std::endl;
        return 1;
    }

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
        std::cerr << "Error during setup: " << e.what() << std::endl;
        return 1;
    }

    try {
        AlleleEventCounter aev (cp, result, AEVSettings{no_overlaps});
        count (aln_in, idx, aev, reg, cp);
    } catch (std::exception &e) {
        std::cerr << "Error during calculation: " << e.what()
                  << std::endl;
        return 1;
    }

    // NOTE: may also want to optionally include rid in output with
    // pos
    try {
        if (print_head) {
            if (print_row)
                std::cout << "pos,";
            std::cout << HEADER << "\n";
        }
        size_t i = 0;
        size_t row_counter = 1; // adds 1 for 1-indexed row to match
                                // input region string
        while (i < result.size()) {
            if (print_row)
                std::cout
                    << static_cast<uint64_t> (reg.start) + row_counter
                    << ",";
            size_t j = 0;
            while (j < (N_FIELDS_PER_OBS - 1)) {
                std::cout << result[i + j] << ",";
                ++j;
            }
            std::cout << result[i + j] << "\n";
            i += N_FIELDS_PER_OBS;
            ++row_counter;
        }
    } catch (std::exception &e) {
        std::cerr << "Error during write: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
