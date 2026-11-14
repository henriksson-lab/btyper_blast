/* Differential driver for unchanged BLAST+ 2.17.0+ lookup kernels.
 * Unused toolkit functions are removed by the linker.
 */
#include <inttypes.h>
#include <strings.h>
#include "blast_nalookup.c"
#include "blast_lookup.c"
#include "lookup_util.c"
#include "blast_util.c"
#include "blast_nascan.c"

static void bytes(const Uint1 *p, int length)
{
    if (!length) putchar('-');
    for (int i = 0; i < length; ++i) printf("%02x", p[i]);
}

static void index_case(Uint1 *sequence, int length, int word, int width,
                       int *left, int *right, int count)
{
    BlastSeqLoc nodes[4] = {{0}};
    SSeqRange ranges[4];
    for (int i = 0; i < count; ++i) {
        ranges[i].left = left[i]; ranges[i].right = right[i];
        nodes[i].ssr = &ranges[i];
        nodes[i].next = i + 1 < count ? &nodes[i+1] : NULL;
    }
    BLAST_SequenceBlk query = {0};
    query.sequence = sequence; query.length = length;
    int size = 1 << (2 * width);
    Int4 **backbone = calloc(size, sizeof(Int4 *));
    BlastLookupIndexQueryExactMatches(backbone, word, 2, width, &query,
                                      count ? nodes : NULL);
    Int4 max_off;
    Int4 estimate = EstimateNumTableEntries(count ? nodes : NULL, &max_off);
    printf("index %d %d ", word, width); bytes(sequence, length);
    printf(" %d", count);
    for (int i = 0; i < count; ++i) printf(" %d %d", left[i], right[i]);
    printf(" %d %d", estimate, max_off);
    int occupied = 0;
    for (int i = 0; i < size; ++i) if (backbone[i]) ++occupied;
    printf(" %d", occupied);
    for (int i = 0; i < size; ++i) if (backbone[i]) {
        Int4 *chain = backbone[i];
        printf(" %d %d %d", i, chain[0], chain[1]);
        for (int j = 0; j < chain[1]; ++j) printf(" %d", chain[j+2]);
        free(chain);
    }
    free(backbone); putchar('\n');
    BlastCompressBlastnaSequence(&query);
    printf("rolling "); bytes(sequence, length); putchar(' ');
    bytes(query.compressed_nuc_seq_start, length + 3); putchar('\n');
    free(query.compressed_nuc_seq_start);
    if (width == 8) {
        query.compressed_nuc_seq_start = query.compressed_nuc_seq = NULL;
        BlastSmallNaLookupTable *lookup = NULL;
        LookupTableOptions options = {0}; options.word_size = word;
        Int4 status = BlastSmallNaLookupTableNew(&query, count ? nodes : NULL,
                                                &lookup, &options, NULL, width);
        printf("small %d %d ", word, width); bytes(sequence, length);
        printf(" %d", count);
        for (int i = 0; i < count; ++i) printf(" %d %d", left[i], right[i]);
        printf(" %d", status);
        if (!status) {
            printf(" %d %d", lookup->longest_chain, lookup->overflow_size);
            int used = 0;
            for (int i = 0; i < lookup->backbone_size; ++i) if (lookup->final_backbone[i] != -1) ++used;
            printf(" %d", used);
            for (int i = 0; i < lookup->backbone_size; ++i) if (lookup->final_backbone[i] != -1)
                printf(" %d %d", i, lookup->final_backbone[i]);
            for (int i = 2; i < lookup->overflow_size; ++i) printf(" %d", lookup->overflow[i]);
            int masks = 0;
            for (BlastSeqLoc *loc = lookup->masked_locations; loc; loc = loc->next) ++masks;
            printf(" %d", masks);
            for (BlastSeqLoc *loc = lookup->masked_locations; loc; loc = loc->next)
                printf(" %d %d", loc->ssr->left, loc->ssr->right);
            putchar('\n');
            if (length >= word && length <= 128 && count == 1) {
                Uint1 packed[260] = {0};
                for (int i = 0; i < length; ++i) packed[i/4] |= (sequence[i]&3) << (6-2*(i%4));
                BLAST_SequenceBlk subject = {0}; subject.sequence = packed; subject.length = length;
                LookupTableWrap wrap = {0}; wrap.lut_type = eSmallNaLookupTable; wrap.lut = lookup;
                int extras[] = {0,1,5,31};
                for (int start = 0; start < 4; ++start)
                for (int e = 0; e < 4; ++e) {
                    int first = word == 11 ? start*4 : start;
                    if (first > length-8) continue;
                    Int4 range[] = {first, length-8};
                    int capacity = MAX(1, lookup->longest_chain + extras[e]);
                    BlastOffsetPair *pairs = calloc(capacity, sizeof(*pairs));
                    do {
                        int previous = range[0];
                        Int4 hits = word == 11 ? s_BlastSmallNaScanSubject_8_4(&wrap, &subject, pairs, capacity, range)
                            : s_BlastSmallNaScanSubject_8_1Mod4(&wrap, &subject, pairs, capacity, range);
                        printf("scan %d %d %d ", previous, range[1], capacity);
                        bytes(packed, (length+3)/4);
                        printf(" %d %d", range[0], hits);
                        for (int i = 0; i < hits; ++i) printf(" %d %d", pairs[i].qs_offsets.q_off, pairs[i].qs_offsets.s_off);
                        putchar('\n');
                        ASSERT(range[0] > previous);
                    } while (range[0] <= range[1]);
                    free(pairs);
                }
            }
            free(query.compressed_nuc_seq_start);
            BlastSmallNaLookupTableDestruct(lookup);
        } else { putchar('\n'); }
    }
}

int main(void)
{
    int entries[] = {0,1,249,250,251,1249,1250,1251,8499,8500,8501,
        11999,12000,12001,17999,18000,18001,20999,21000,21001,
        32766,32767,32768,59999,60000,60001,179999,180000,180001,
        299999,300000,300001,899999,900000,900001};
    int offsets[] = {0,32766,32767,32768,32769,1000000};
    int words[] = {4,5,6,7,8,9,10,11,12,13,16,28,64};
    LookupTableOptions options = {0};
    for (int w = 0; w < (int)(sizeof(words)/sizeof(words[0])); ++w)
    for (int e = 0; e < (int)(sizeof(entries)/sizeof(entries[0])); ++e)
    for (int o = 0; o < (int)(sizeof(offsets)/sizeof(offsets[0])); ++o)
    for (int mode = 0; mode < 5; ++mode) {
        options.word_size = words[w];
        options.program_number = mode >= 2 ? eBlastTypeMapping : eBlastTypeBlastn;
        options.db_filter = mode == 1 || mode == 3;
        options.mb_template_length = mode == 4 ? 18 : 0;
        Int4 width;
        ELookupTableType type = BlastChooseNaLookupTable(&options, entries[e], offsets[o], &width);
        printf("choose %d %d %d %d %d %d %d %d\n", options.word_size,
            options.program_number, options.db_filter, options.mb_template_length,
            entries[e], offsets[o], type, width);
    }
    Uint1 sequence[1024] = {0};
    int left[4] = {0}, right[4] = {0};
    /* Exhaustive short words with ambiguity/fence values and exact final
       boundaries. Width 2 exercises the shared algorithm economically. */
    Uint1 alphabet[] = {0,1,2,3,4,14,15};
    for (int length = 0; length <= 5; ++length) {
        int combinations = 1;
        for (int i = 0; i < length; ++i) combinations *= 7;
        for (int code = 0; code < combinations; ++code) {
            int value = code;
            for (int i = 0; i < length; ++i) { sequence[i] = alphabet[value%7]; value /= 7; }
            right[0] = length - 1;
            index_case(sequence, length, 4, 2, left, right, length ? 1 : 0);
        }
    }
    /* Standard-task widths and chain growth, disjoint/overlapping ranges,
       unsorted traversal, mask-sized intervals, ambiguity at word edges. */
    uint32_t state = 183004u;
    for (int length = 1; length <= 1024; length = length < 40 ? length+1 : length*2) {
        for (int pattern = 0; pattern < 4; ++pattern) {
            for (int i = 0; i < length; ++i) {
                state = state * 1664525u + 1013904223u;
                sequence[i] = pattern == 0 ? 0 : pattern == 1 ? i%4 : (state>>28)% (pattern == 2 ? 4 : 16);
            }
            for (int word = 11; word <= 28; word += 17)
            for (int mode = 0; mode < 4; ++mode) {
                left[0] = 0; right[0] = length-1;
                int count = 1;
                if (mode == 1) { right[0] = length/2; left[1] = length/2; right[1] = length-1; count = 2; }
                if (mode == 2) { left[1] = 0; right[1] = length-1; count = 2; }
                if (mode == 3) { left[0] = length/2; left[1] = 0; right[1] = length/2; count = 2; }
                index_case(sequence, length, word, 8, left, right, count);
            }
        }
    }
    left[0] = 0; right[0] = 63;
    index_case(sequence, 64, 11, 8, left, right, 0);
    index_case(sequence, 64, 28, 8, left, right, 0);
    index_case(sequence, 64, 8, 8, left, right, 1);
    /* Exercise the constructor's propagated overflow failure as well. */
    Uint1 *large_query = calloc(32772, 1);
    right[0] = 32771;
    index_case(large_query, 32772, 11, 8, left, right, 1);
    free(large_query);
    /* Overflow cursors must fit in 15 bits, including chain terminators.
       Two moderate chains can exceed the limit despite each fitting alone. */
    int hit_counts[][2] = {{0,0},{1,0},{32763,0},{32764,0},{32765,0},
        {16381,16382},{16382,16382},{16382,16383}};
    for (int c = 0; c < 8; ++c) {
        Int4 *thin[4] = {0};
        for (int k = 0; k < 2; ++k) if (hit_counts[c][k]) {
            int hits = hit_counts[c][k];
            thin[k] = calloc(hits+2, sizeof(Int4));
            thin[k][0] = hits+2; thin[k][1] = hits;
            for (int j = 0; j < hits; ++j) thin[k][j+2] = j;
        }
        BLAST_SequenceBlk query = {0}; query.sequence = sequence; query.length = 40;
        BlastSmallNaLookupTable lookup = {0}; lookup.backbone_size = 4;
        int status = s_BlastSmallNaLookupFinalize(thin, &lookup, &query);
        printf("finalize %d %d %d %d %d", hit_counts[c][0], hit_counts[c][1],
               status, lookup.longest_chain, lookup.overflow_size);
        for (int i = 0; i < 4; ++i) { ASSERT(!thin[i]); printf(" %d", lookup.final_backbone ? lookup.final_backbone[i] : -99); }
        if (!status) for (int i = 2; i < lookup.overflow_size; ++i) printf(" %d", lookup.overflow[i]);
        printf(" %d\n", query.compressed_nuc_seq != NULL);
        free(query.compressed_nuc_seq_start); free(lookup.final_backbone); free(lookup.overflow);
    }
    return 0;
}
