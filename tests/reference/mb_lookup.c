/* Test unchanged original megablast constructors and scanners. */
#include <inttypes.h>
#include <strings.h>
#include "blast_nalookup.c"
#include "blast_lookup.c"
#include "blast_nascan.c"
#include "blast_util.c"
#include "lookup_util.c"

static void bytes(const Uint1 *p, int length)
{
    for (int i = 0; i < length; ++i) printf("%02x", p[i]);
}

int main(void)
{
    int configs[][2] = {{11,10},{11,11},{28,11},{28,12},{9,9},{13,10},{14,11},{15,12}};
    uint32_t state = 8927143u;
    for (int c = 0; c < 8; ++c)
    for (int compression = 0; compression < 6; ++compression)
    for (int pattern = 0; pattern < 4; ++pattern) {
        int length = 96 + pattern*9;
        Uint1 query_data[130] = {0};
        for (int i = 0; i < length; ++i) {
            state = state * 1664525u + 1013904223u;
            query_data[i+1] = pattern == 0 ? 0 : pattern == 1 ? i%4 : (state>>28)%4;
            if (pattern == 3 && i%17 == 0) query_data[i+1] = 14;
        }
        BLAST_SequenceBlk query = {0}; query.sequence_start = query_data;
        query.sequence = query_data+1; query.length = length;
        SSeqRange ranges[3] = {{0}}; BlastSeqLoc locations[3] = {{0}};
        ranges[0].left = 0; ranges[0].right = length-1;
        int location_count = 1;
        if (pattern == 3) {
            ranges[0].left = 5; ranges[0].right = 42;
            ranges[1].left = 55; ranges[1].right = 83;
            ranges[2].left = 98; ranges[2].right = length-1; location_count = 3;
        }
        for (int i = 0; i < location_count; ++i) {
            locations[i].ssr = &ranges[i]; locations[i].next = i+1 < location_count ? &locations[i+1] : NULL;
        }
        LookupTableOptions options = {0}; options.word_size = configs[c][0]; options.program_number = eBlastTypeBlastn;
        int estimates[] = {10000,15001,15000,799999,800000,800001};
        int approx = estimates[compression];
        QuerySetUpOptions qopts = {0}; SBlastFilterOptions filters = {0};
        qopts.filtering_options = &filters; filters.mask_at_hash = pattern % 2;
        BlastMBLookupTable *built = NULL;
        ASSERT(BlastMBLookupTableNew(&query,locations,&built,&options,&qopts,approx,configs[c][1],NULL) == 0);
        BlastMBLookupTable lookup = *built; free(built);
        int pv_size = lookup.hashsize >> lookup.pv_array_bts;
        LookupTableWrap wrap = {0}; wrap.lut_type = eMBLookupTable; wrap.lut = &lookup;
        s_MBChooseScanSubject(&wrap);
        int selected;
        if (lookup.scansub_callback == (void*)s_MBScanSubject_10_2) selected = 1;
        else if (lookup.scansub_callback == (void*)s_MBScanSubject_11_1Mod4) selected = 2;
        else if (lookup.scansub_callback == (void*)s_MBScanSubject_11_2Mod4) selected = 3;
        else { ASSERT(lookup.scansub_callback == (void*)s_MBScanSubject_Any); selected = 0; }
        printf("mb %d %d ",approx,filters.mask_at_hash); bytes(query.sequence,length);
        printf(" %d",location_count);
        for (int i = 0; i < location_count; ++i) printf(" %d %d",ranges[i].left,ranges[i].right);
        putchar(' ');
        printf("%d %d %d %d %d %d %d", lookup.word_length,lookup.lut_word_length,lookup.scan_step,
            length,lookup.longest_chain,lookup.pv_array_bts,selected);
        int used = 0;
        for (Int8 i = 0; i < lookup.hashsize; ++i) if (lookup.hashtable[i]) ++used;
        printf(" %d",used);
        for (Int8 i = 0; i < lookup.hashsize; ++i) if (lookup.hashtable[i]) printf(" %" PRId64 " %d",i,lookup.hashtable[i]);
        for (int i = 0; i <= length; ++i) printf(" %d",lookup.next_pos[i]);
        used = 0;
        for (int i = 0; i < pv_size; ++i) if (lookup.pv_array[i]) ++used;
        printf(" %d %d",pv_size,used);
        for (int i = 0; i < pv_size; ++i) if (lookup.pv_array[i]) printf(" %d %" PRIu32,i,lookup.pv_array[i]);
        int masks = 0;
        for (BlastSeqLoc *loc = lookup.masked_locations; loc; loc = loc->next) ++masks;
        printf(" %d",masks);
        for (BlastSeqLoc *loc = lookup.masked_locations; loc; loc = loc->next) printf(" %d %d",loc->ssr->left,loc->ssr->right);
        putchar('\n');
        BlastOffsetPair probes[130];
        for (Int8 i = 0; i < lookup.hashsize; ++i) if (lookup.hashtable[i]) {
            Int8 indices[] = {i, i^32};
            for (int k = 0; k < 2; ++k) {
                Int4 hits = s_BlastMBLookupRetrieve(&lookup,indices[k],probes,79);
                printf("mb_probe %" PRId64 " %d %d",indices[k],s_BlastMBLookupHasHits(&lookup,indices[k]),hits);
                for (int j = 0; j < hits; ++j) printf(" %d %d",probes[j].qs_offsets.q_off,probes[j].qs_offsets.s_off);
                putchar('\n');
            }
        }
        Uint1 packed[36] = {0};
        for (int i = 0; i < length; ++i) packed[i/4] |= (query_data[i+1]&3) << (6-2*(i%4));
        BLAST_SequenceBlk subject = {0}; subject.sequence = packed; subject.length = length;
        int capacities[] = {0,1,lookup.longest_chain,lookup.longest_chain+1,lookup.longest_chain+7,(length+1)*(length+1)};
        for (int generic = 0; generic < 2; ++generic)
        for (int mask = 0; mask < 3; ++mask)
        for (int first = 0; first < 4; ++first)
        for (int cap = 0; cap < 6; ++cap) {
            if (compression >= 2) continue; /* Additional estimates test construction/PV boundaries. */
            if (cap == 5 && (mask != 0 || first != 0)) continue;
            subject.mask_type = mask;
            int capacity = capacities[cap];
            BlastOffsetPair *pairs = calloc(MAX(1,capacity),sizeof(*pairs));
            Int4 scan_range[] = {first,length-lookup.lut_word_length};
            TNaScanSubjectFunction scanner = generic ? s_MBScanSubject_Any : (TNaScanSubjectFunction)lookup.scansub_callback;
            int previous, calls = 0;
            do {
                previous = scan_range[0];
                Int4 hits = scanner(&wrap,&subject,pairs,capacity,scan_range);
                printf("mb_scan %d %d %d %d %d ",generic,mask,previous,scan_range[1],capacity);
                bytes(packed,(length+3)/4+2);
                printf(" %d %d",scan_range[0],hits);
                for (int i = 0; i < hits; ++i) printf(" %d %d",pairs[i].qs_offsets.q_off,pairs[i].qs_offsets.s_off);
                putchar('\n');
                calls++;
            } while (calls < 4 && scan_range[0] > previous && scan_range[0] <= scan_range[1]);
            free(pairs);
        }
        free(lookup.hashtable); free(lookup.next_pos); free(lookup.pv_array); BlastSeqLocFree(lookup.masked_locations);
    }
    for (int bit = 0; bit < 63; ++bit) {
        Int8 value = (Int8)1 << bit;
        printf("ilog2 %" PRId64 " %d\n",value,ilog2(value));
        if (value > 1) printf("ilog2 %" PRId64 " %d\n",value-1,ilog2(value-1));
        printf("ilog2 %" PRId64 " %d\n",value+1,ilog2(value+1));
    }
    printf("ilog2 0 %d\n",ilog2(0));
    printf("ilog2 %" PRId64 " %d\n",INT64_MAX,ilog2(INT64_MAX));
    BLAST_SequenceBlk empty_query = {0}; LookupTableOptions empty_options = {0};
    empty_options.word_size = 28; empty_options.program_number = eBlastTypeBlastn;
    BlastMBLookupTable *empty = NULL;
    printf("mb_empty %d\n",BlastMBLookupTableNew(&empty_query,NULL,&empty,&empty_options,NULL,0,11,NULL));
    return 0;
}
