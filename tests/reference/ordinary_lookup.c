/* Run original NCBI BLAST+ 2.17.0+ ordinary nucleotide lookup functions. */
#include <inttypes.h>
#include <strings.h>
#include "blast_nalookup.c"
#include "blast_lookup.c"
#include "blast_nascan.c"
#include "blast_util.c"

static void bytes(const Uint1 *p, int length)
{
    if (!length) putchar('-');
    for (int i = 0; i < length; ++i) printf("%02x", p[i]);
}

static QuerySetUpOptions *query_options(int mode, QuerySetUpOptions *opts, SBlastFilterOptions *filters)
{
    memset(opts, 0, sizeof(*opts)); memset(filters, 0, sizeof(*filters));
    if (mode == 0) return NULL;
    if (mode == 2 || mode == 3 || mode == 9 || mode == 10) opts->filtering_options = filters;
    filters->mask_at_hash = mode == 3 || mode == 9;
    switch (mode) {
    case 4: case 10: opts->filter_string = "m"; break;
    case 5: case 9: opts->filter_string = "F"; break;
    case 6: opts->filter_string = "F; m "; break;
    case 7: opts->filter_string = "am\0ignored"; break;
    case 8: opts->filter_string = "F\0m"; break;
    }
    return opts;
}

static void build_case(Uint1 *sequence, int length, int word, int width, int mode,
                      int *left, int *right, int count, int scan)
{
    BlastSeqLoc nodes[4] = {{0}};
    SSeqRange ranges[4];
    for (int i = 0; i < count; ++i) {
        ranges[i].left = left[i]; ranges[i].right = right[i];
        nodes[i].ssr = &ranges[i]; nodes[i].next = i+1 < count ? &nodes[i+1] : NULL;
    }
    QuerySetUpOptions opts; SBlastFilterOptions filters;
    QuerySetUpOptions *qopts = query_options(mode, &opts, &filters);
    LookupTableOptions lookup_options = {0}; lookup_options.word_size = word;
    if (length > 32768) {
        lookup_options.program_number = eBlastTypeBlastn;
        Int4 estimated = 0, maximum = 0, selected_width;
        for (int i = 0; i < count; ++i) {
            estimated += right[i]-left[i]; maximum = MAX(maximum,right[i]);
        }
        ASSERT(BlastChooseNaLookupTable(&lookup_options, estimated, maximum, &selected_width) == eNaLookupTable);
        ASSERT(selected_width == width);
    }
    BLAST_SequenceBlk query = {0}; query.sequence = sequence; query.length = length;
    BlastNaLookupTable *lookup;
    ASSERT(BlastNaLookupTableNew(&query, count ? nodes : NULL, &lookup, &lookup_options, qopts, width) == 0);
    LookupTableWrap wrap = {0}; wrap.lut_type = eNaLookupTable; wrap.lut = lookup;
    s_NaChooseScanSubject(&wrap);
    printf("na %d %d %d ", word, width, mode); bytes(sequence, length);
    printf(" %d", count);
    for (int i = 0; i < count; ++i) printf(" %d %d", left[i], right[i]);
    printf(" %d %d", lookup->longest_chain, lookup->overflow_size);
    int used = 0;
    for (int i = 0; i < lookup->backbone_size; ++i) if (lookup->thick_backbone[i].num_used) ++used;
    printf(" %d", used);
    for (int i = 0; i < lookup->backbone_size; ++i) if (lookup->thick_backbone[i].num_used) {
        NaLookupBackboneCell *cell = &lookup->thick_backbone[i];
        printf(" %d %d %d %d %d", i, cell->num_used,
            cell->payload.entries[0], cell->payload.entries[1], cell->payload.entries[2]);
        ASSERT(s_BlastLookupGetNumHits(lookup, i) == cell->num_used);
    }
    for (int i = 0; i < lookup->overflow_size; ++i) printf(" %d", lookup->overflow[i]);
    int pv_len = (lookup->backbone_size >> PV_ARRAY_BTS) + 1;
    used = 0;
    for (int i = 0; i < pv_len; ++i) if (lookup->pv[i]) ++used;
    printf(" %d %d", pv_len, used);
    for (int i = 0; i < pv_len; ++i) if (lookup->pv[i]) printf(" %d %" PRIu32, i, lookup->pv[i]);
    int masks = 0;
    for (BlastSeqLoc *loc = lookup->masked_locations; loc; loc = loc->next) ++masks;
    printf(" %d", masks);
    for (BlastSeqLoc *loc = lookup->masked_locations; loc; loc = loc->next)
        printf(" %d %d", loc->ssr->left, loc->ssr->right);
    printf(" %d\n", lookup->scansub_callback == (void *)s_BlastNaScanSubject_8_4);
    if (scan && length >= width) {
        Uint1 *packed = calloc((length+3)/4 + 2, 1);
        for (int i = 0; i < length; ++i) packed[i/4] |= (sequence[i]&3) << (6-2*(i%4));
        BLAST_SequenceBlk subject = {0}; subject.sequence = packed; subject.length = length;
        int capacities[] = {0, 1, MAX(1, lookup->longest_chain-1), MAX(1,lookup->longest_chain), lookup->longest_chain+5};
        for (int kind = 0; kind < 2; ++kind)
        for (int mask = 0; mask < 3; ++mask)
        for (int start = 0; start < 4; ++start)
        for (int c = 0; c < 5; ++c) {
            subject.mask_type = mask;
            int capacity = capacities[c];
            BlastOffsetPair *pairs = calloc(MAX(1,capacity), sizeof(*pairs));
            Int4 range[] = {start, length-width};
            TNaScanSubjectFunction scanner = kind == 0 ? (TNaScanSubjectFunction)lookup->scansub_callback : s_BlastNaScanSubject_Any;
            /* In the aligned no-masking path, original code rounds down
               unaligned start offsets. Keep those cases too. */
            Int4 hits = scanner(&wrap, &subject, pairs, capacity, range);
            printf("na_scan %d %d %d %d %d ", kind, mask, start, length-width, capacity);
            bytes(packed, (length+3)/4 + 2);
            printf(" %d %d", range[0], hits);
            for (int i = 0; i < hits; ++i) printf(" %d %d", pairs[i].qs_offsets.q_off, pairs[i].qs_offsets.s_off);
            putchar('\n'); free(pairs);
        }
        free(packed);
    }
    BlastNaLookupTableDestruct(lookup);
}

int main(void)
{
    for (int mode = 0; mode <= 10; ++mode) {
        QuerySetUpOptions opts; SBlastFilterOptions filters;
        QuerySetUpOptions *qopts = query_options(mode, &opts, &filters);
        printf("mask_at_hash %d %d\n", mode, s_HasMaskAtHashEnabled(qopts));
    }
    Uint1 sequence[128] = {0}; uint32_t state = 5831809u;
    int left[4] = {0}, right[4] = {0};
    /* Every loop branch and in-cell/overflow boundary. Widths 4..7 only
       exercise shared code; standard presets select ordinary width 8. */
    for (int width = 4; width <= 8; ++width)
    for (int stride = 1; stride <= 5; ++stride)
    for (int length = 8; length <= 40; length += 4)
    for (int pattern = 0; pattern < 3; ++pattern) {
        for (int i = 0; i < length; ++i) {
            state = state * 1664525u + 1013904223u;
            sequence[i] = pattern == 0 ? 0 : pattern == 1 ? i%4 : (state>>28);
        }
        right[0] = length-1;
        build_case(sequence,length,width+stride-1,width,3,left,right,1,1);
    }
    for (int mode = 0; mode <= 10; ++mode)
    for (int word = 11; word <= 28; word += 17)
    for (int intervals = 0; intervals < 5; ++intervals) {
        for (int i = 0; i < 128; ++i) sequence[i] = i%4;
        left[0] = 0; right[0] = 127; int count = 1;
        if (intervals == 0) count = 0;
        if (intervals == 2) { left[0] = 5; right[0] = 37; left[1] = 53; right[1] = 90; count = 2; }
        if (intervals == 3) { left[1] = 0; right[1] = 127; count = 2; }
        if (intervals == 4) { left[0] = 60; left[1] = 0; right[1] = 70; count = 2; }
        build_case(sequence,128,word,8,mode,left,right,count,intervals == 2 && mode == 3);
    }
    /* A legal ordinary-table selection: long query, few indexed words,
       offsets crossing 32768. Also retain offsets well beyond 65535. */
    Uint1 *large = malloc(100128); memset(large,14,100128);
    int positions[] = {32760,32768,65528,100000};
    for (int i = 0; i < 4; ++i) {
        left[i] = positions[i]; right[i] = positions[i]+31;
        for (int j = left[i]; j <= right[i]; ++j) large[j] = j%4;
    }
    for (int word = 11; word <= 28; word += 17) build_case(large,100128,word,8,3,left,right,4,0);
    free(large);
    return 0;
}
