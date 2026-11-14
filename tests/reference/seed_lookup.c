/* Standard-preset nucleotide lookup orchestration. All indexing, construction,
 * selection and scan algorithms execute their unchanged BLAST+ 2.17.0+ bodies.
 * This driver supplies the nucleotide branch of lookup_wrap.c orchestration;
 * it does not claim to execute the full multi-program wrapper or search engine.
 */
#include <inttypes.h>
#include <strings.h>
#include "blast_nalookup.c"
#include "blast_lookup.c"
#include "blast_nascan.c"
#include "blast_util.c"
#include "lookup_util.c"

static void bytes(const Uint1 *p, int length)
{
    for (int i = 0; i < length; ++i) printf("%02x",p[i]);
}

static void query_case(int length, int word, int mode)
{
    Uint1 *data = calloc(length+2,1);
    uint32_t state = length+79u;
    for (int i = 0; i < length; ++i) {
        state = state*1664525u+1013904223u;
        data[i+1] = mode == 1 ? 0 : (state>>28)%4;
    }
    BLAST_SequenceBlk query = {0}; query.sequence_start = data; query.sequence = data+1; query.length = length;
    BlastSeqLoc locations[3] = {{0}}; SSeqRange ranges[3] = {{0}};
    int count = 1;
    ranges[0].left = 0; ranges[0].right = length-1;
    if (mode == 2) { ranges[0].left = 40000; ranges[0].right = 40063; }
    if (mode == 3) { ranges[0].left = length; ranges[0].right = length-1; }
    if (mode == 4) {
        ranges[0].left = 0; ranges[0].right = 7;
        ranges[1].left = 15; ranges[1].right = 25;
        ranges[2].left = 40; ranges[2].right = length-1; count = 3;
    }
    if (mode == 5) count = 0;
    for (int i = 0; i < count; ++i) { locations[i].ssr = &ranges[i]; locations[i].next = i+1 < count ? &locations[i+1] : NULL; }
    LookupTableOptions *options = NULL;
    ASSERT(LookupTableOptionsNew(eBlastTypeBlastn,&options) == 0); options->word_size = word;
    QuerySetUpOptions qopts = {0}; SBlastFilterOptions filters = {0};
    filters.mask_at_hash = TRUE; qopts.filtering_options = &filters;
    Int4 max_offset, width;
    int entries = EstimateNumTableEntries(count ? locations : NULL,&max_offset);
    LookupTableWrap wrap = {0};
    wrap.lut_type = BlastChooseNaLookupTable(options,entries,max_offset,&width);
    int status;
    if (wrap.lut_type == eSmallNaLookupTable) {
        status = BlastSmallNaLookupTableNew(&query,count ? locations : NULL,
            (BlastSmallNaLookupTable**)&wrap.lut,options,&qopts,width);
        if (status != 0) {
            wrap.lut_type = eNaLookupTable;
            status = BlastNaLookupTableNew(&query,count ? locations : NULL,
                (BlastNaLookupTable**)&wrap.lut,options,&qopts,width);
        }
    } else if (wrap.lut_type == eNaLookupTable) {
        status = BlastNaLookupTableNew(&query,count ? locations : NULL,
            (BlastNaLookupTable**)&wrap.lut,options,&qopts,width);
    } else {
        ASSERT(wrap.lut_type == eMBLookupTable);
        status = BlastMBLookupTableNew(&query,count ? locations : NULL,
            (BlastMBLookupTable**)&wrap.lut,options,&qopts,entries,width,NULL);
    }
    ASSERT(status == 0);
    BlastChooseNucleotideScanSubject(&wrap);
    int longest, step; TNaScanSubjectFunction scanner;
    if (wrap.lut_type == eSmallNaLookupTable) {
        BlastSmallNaLookupTable *table = wrap.lut;
        longest = table->longest_chain; step = table->scan_step; scanner = (TNaScanSubjectFunction)table->scansub_callback;
    } else if (wrap.lut_type == eNaLookupTable) {
        BlastNaLookupTable *table = wrap.lut;
        longest = table->longest_chain; step = table->scan_step; scanner = (TNaScanSubjectFunction)table->scansub_callback;
    } else {
        BlastMBLookupTable *table = wrap.lut;
        longest = table->longest_chain; step = table->scan_step; scanner = (TNaScanSubjectFunction)table->scansub_callback;
    }
    printf("seed_query %d ",word); bytes(query.sequence,length); printf(" %d",count);
    for (int i = 0; i < count; ++i) printf(" %d %d",ranges[i].left,ranges[i].right);
    printf(" %d %d %d %d\n",wrap.lut_type,width,step,longest);
    int lengths[] = {1,2,3,7,8,9,10,11,12,15,16,27,28,29,31,32,33,63,64,65,127,128,129};
    for (int k = 0; k < 23; ++k) {
        int subject_len = lengths[k]; Uint1 raw[129]; Uint1 packed[37] = {0};
        for (int i = 0; i < subject_len; ++i) {
            raw[i] = query.sequence[(mode == 2 ? 40000+i : i)%length]&3;
            packed[i/4] |= raw[i] << (6-2*(i%4));
        }
        packed[subject_len/4] |= subject_len%4;
        BLAST_SequenceBlk subject = {0}; subject.sequence = packed; subject.length = subject_len;
        int minimum = longest + (wrap.lut_type == eMBLookupTable);
        int capacities[] = {MAX(1,minimum),minimum+7,4096+longest};
        for (int cap = 0; cap < 3; ++cap) {
            int capacity = capacities[cap];
            BlastOffsetPair *pairs = calloc(capacity,sizeof(*pairs));
            printf("seed_subject "); bytes(raw,subject_len); printf(" %d\n",capacity);
            Int4 scan_range[] = {0,subject_len-width};
            while (scan_range[0] <= scan_range[1]) {
                int previous = scan_range[0];
                Int4 hits = scanner(&wrap,&subject,pairs,capacity,scan_range);
                printf("seed_batch %d %d %d",previous,scan_range[0],hits);
                for (int i = 0; i < hits; ++i) printf(" %d %d",pairs[i].qs_offsets.q_off,pairs[i].qs_offsets.s_off);
                putchar('\n'); ASSERT(scan_range[0] > previous);
            }
            printf("seed_done %d\n",scan_range[0]);
            free(pairs);
        }
    }
    if (wrap.lut_type == eSmallNaLookupTable) BlastSmallNaLookupTableDestruct(wrap.lut);
    else if (wrap.lut_type == eNaLookupTable) BlastNaLookupTableDestruct(wrap.lut);
    else BlastMBLookupTableDestruct(wrap.lut);
    free(query.compressed_nuc_seq_start); free(options); free(data);
}

int main(void)
{
    for (int word = 11; word <= 28; word += 17) {
        for (int mode = 0; mode < 6; ++mode) query_case(mode == 2 || mode == 3 ? 100128 : 64,word,mode);
    }
    int lengths[] = {8500,8501,12000,12001,15001,15002,180000,180001,300000,300001,800001};
    for (int i = 0; i < 11; ++i) query_case(lengths[i],i < 2 || i > 7 ? 28 : 11,0);
    return 0;
}
