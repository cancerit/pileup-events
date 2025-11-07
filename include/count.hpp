#include "bounds.hpp"
#include "pileup.hpp"
#include "structs.hpp"

// nothing but C please
extern "C" {
struct pf_capture {
    htsFile *fh;
    hts_itr_t *it;
    const count_params *p;
};
inline int pileup_func (void *data,
                        bam1_t *b) {
    pf_capture *d = static_cast<pf_capture *> (data);
    int ret;
    // find the next good read
    while (1) {
        ret = sam_itr_next (d->fh, d->it, b);
        if (ret < 0) {
            break; // EOF/err
        }
        if (!(b->core.flag & d->p->exclude_flag) &&
            ((b->core.flag & d->p->include_flag) ==
             d->p->include_flag) &&
            b->core.qual >= d->p->min_mapq) {
            break; // found good read
        };
    }
    return ret;
};
}
// end nothing but C

// bam2R
// NOTE: does not at present include the max_mismatches functionality
// added to recent versions of deepsnv
inline void count (htsFile *aln_fh,
                   hts_idx_t *aln_idx,
                   AlleleEventCounter ctr,
                   const hts_region reg,
                   const count_params params) {
    bam_plp_t buf = NULL;
    bam1_t *b = NULL;
    bam_hdr_t *head = NULL;

    safe_size_opts sso_plp_pos;
    sso_plp_pos.msg = "error translating htslib pileup position into "
                      "appropriate index for results array";

    // fetch a read overlapping the query region;
    // then do a pileup per base for the total region
    // covered by the retrieved read;
    // then count events on those pileups which overlap
    // the original query region.
    hts_itr_t *iter =
        sam_itr_queryi (aln_idx, reg.rid, reg.start, reg.end);

    pf_capture pfc{aln_fh, iter, &params};
    buf = bam_plp_init (pileup_func,
                        &pfc); // initialize pileup
    bam_plp_set_maxcnt (buf, params.max_depth);

    int64_t plp_pos = -1;
    int plp_tid = -1, n_plp = -1;
    const bam_pileup1_t *pl;
    size_t pos_offset;
    while ((pl = bam_plp64_auto (buf, &plp_tid, &plp_pos, &n_plp)) !=
           0) {
        if (n_plp < 0 || plp_tid < 0 || plp_pos < 0) {
            throw std::runtime_error ("pileup failed");
        }
        if (!(plp_pos >= reg.start && plp_pos < reg.end)) {
            continue;
        }
        pos_offset = safe_size (plp_pos - reg.start, sso_plp_pos);
        ctr.count_pileup (pl, pos_offset, safe_size (n_plp));
    }

    sam_itr_destroy (iter);
    bam_destroy1 (b);
    bam_hdr_destroy (head);
    bam_plp_destroy (buf);
}
