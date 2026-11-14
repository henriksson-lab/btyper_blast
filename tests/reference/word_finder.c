/* Original file-local nucleotide kernels are included unchanged. The test
 * driver supplies arrays/options, then unused toolkit entry points are stripped.
 */
#include <inttypes.h>
#include <strings.h>
#include "na_ungapped.c"
#include "blast_nalookup.c"
#include "blast_lookup.c"
#include "blast_util.c"
#include "blast_extend.c"
#include "blast_nascan.c"

static uint32_t state=1605347u;
static unsigned next(void) {state=state*1664525u+1013904223u;return state>>28;}
static void bytes(const Uint1* p,int length)
{
    if (!length) putchar('-');
    for (int i=0;i<length;++i) printf("%02x",p[i]);
}
static void result(const BlastUngappedData* data)
{
    printf(" %d %d %d %d\n",data->q_start,data->s_start,data->length,data->score);
}
static void matrix_case(int reward,int penalty,Int4 cells[16][16],Int4* rows[16],int table[256])
{
    for (int i=0;i<16;++i) {rows[i]=cells[i];for (int j=0;j<16;++j) cells[i][j]=12345;}
    SBlastScoreMatrix matrix={0};matrix.data=rows;matrix.nrows=matrix.ncols=16;
    BlastScoreBlk sbp={0};sbp.reward=reward;sbp.penalty=penalty;sbp.alphabet_size=16;sbp.matrix=&matrix;
    ASSERT(BlastScoreBlkNuclMatrixCreate(&sbp)==0);
    printf("matrix %d %d",reward,penalty);
    for (int i=0;i<16;++i) for (int j=0;j<16;++j) printf(" %d",cells[i][j]);
    putchar('\n');
    /* Exact block from BlastInitialWordParametersNew, supplied as an input
       score table here; no Rust parameter-constructor parity is claimed. */
    for (int i=0;i<256;++i) {
        int score=0;
        if (i&3) score+=penalty;else score+=reward;
        if ((i>>2)&3) score+=penalty;else score+=reward;
        if ((i>>4)&3) score+=penalty;else score+=reward;
        if (i>>6) score+=penalty;else score+=reward;
        table[i]=score;
    }
    printf("score_table");for (int i=0;i<256;++i) printf(" %d",table[i]);putchar('\n');
}
typedef Int4 (*MiniExtendFunction)(const BlastOffsetPair*, Int4, const BlastInitialWordParameters*, LookupTableWrap*, BLAST_SequenceBlk*, BLAST_SequenceBlk*, Int4**, BlastQueryInfo*, Blast_ExtendWord*, BlastInitHitList*, Uint4);
static void finder_cases(void)
{
    Int4 cells[16][16],*cells_ptr[16];int table[256];matrix_case(1,-2,cells,cells_ptr,table);
    int kinds[][3]={{eSmallNaLookupTable,11,8},{eSmallNaLookupTable,28,8},
        {eNaLookupTable,11,8},{eNaLookupTable,28,8},{eMBLookupTable,11,10},
        {eMBLookupTable,11,11},{eMBLookupTable,28,11},{eMBLookupTable,28,12},{eSmallNaLookupTable,8,8},{eSmallNaLookupTable,9,8}};
    for (int k=0;k<8;++k) for (int pattern=0;pattern<4;++pattern) {
        Uint1 qb[259]={0},packed[69]={0};qb[0]=qb[258]=15;
        for (int i=0;i<257;++i) {
            qb[i+1]=pattern==0?0:pattern==1?i%4:next()%4;
            if (pattern==3 && i%23==5) qb[i+1]=14;
        }
        qb[129]=15;
        for (int i=0;i<257;++i) packed[i/4]|=(qb[i+1]&3) << (6-2*(i%4));
        packed[64]|=1;
        BLAST_SequenceBlk query={0},subject={0};query.sequence=qb+1;query.sequence_start=qb;query.length=257;
        subject.sequence=packed;subject.length=257;
        BlastContextInfo contexts[2]={{0}};contexts[0].query_length=128;contexts[1].query_offset=129;contexts[1].query_length=128;
        BlastQueryInfo info={0};info.contexts=contexts;info.last_context=1;info.num_queries=1;info.min_length=info.max_length=128;
        BlastSeqLoc locations[4]={{0}};SSeqRange ranges[4]={{0,127},{129,256}};int count=2;
        if (pattern&1) {count=4;ranges[0]=(SSeqRange){4,48};ranges[1]=(SSeqRange){56,127};ranges[2]=(SSeqRange){133,199};ranges[3]=(SSeqRange){215,256};}
        for (int i=0;i<count;++i) {locations[i].ssr=&ranges[i];locations[i].next=i+1<count?&locations[i+1]:NULL;}
        LookupTableOptions options={0};options.word_size=kinds[k][1];options.program_number=eBlastTypeBlastn;
        QuerySetUpOptions qopts={0};SBlastFilterOptions filters={0};filters.mask_at_hash=TRUE;qopts.filtering_options=&filters;
        LookupTableWrap wrap={0};wrap.lut_type=kinds[k][0];int width=kinds[k][2];BlastSeqLoc* masked=NULL;
        if (wrap.lut_type==eSmallNaLookupTable) {
            ASSERT(BlastSmallNaLookupTableNew(&query,locations,(BlastSmallNaLookupTable**)&wrap.lut,&options,&qopts,width)==0);
            wrap.lookup_callback=(void*)s_SmallNaLookup;masked=((BlastSmallNaLookupTable*)wrap.lut)->masked_locations;
        } else if (wrap.lut_type==eNaLookupTable) {
            ASSERT(BlastNaLookupTableNew(&query,locations,(BlastNaLookupTable**)&wrap.lut,&options,&qopts,width)==0);
            wrap.lookup_callback=(void*)s_NaLookup;masked=((BlastNaLookupTable*)wrap.lut)->masked_locations;
        } else {
            ASSERT(BlastMBLookupTableNew(&query,locations,(BlastMBLookupTable**)&wrap.lut,&options,&qopts,20000,width,NULL)==0);
            wrap.lookup_callback=(void*)s_MBLookup;masked=((BlastMBLookupTable*)wrap.lut)->masked_locations;
        }
        printf("lookup %d %d %d ",wrap.lut_type,options.word_size,width);bytes(query.sequence,257);putchar(' ');bytes(packed,69);printf(" %d",count);
        for (int i=0;i<count;++i) printf(" %d %d",ranges[i].left,ranges[i].right);printf(" %d\n",masked!=NULL);



        BlastChooseNucleotideScanSubject(&wrap);BlastChooseNaExtend(&wrap);
        int longest=wrap.lut_type==eSmallNaLookupTable?((BlastSmallNaLookupTable*)wrap.lut)->longest_chain:wrap.lut_type==eNaLookupTable?((BlastNaLookupTable*)wrap.lut)->longest_chain:((BlastMBLookupTable*)wrap.lut)->longest_chain;
        int profiles[][6]={{0,0,1,50,25,0},{40,3,1,300,200,0},{0,0,0,0,0,0}};
        for(int hash=0;hash<2;hash++)for(int profile=0;profile<3;profile++)for(int cap=0;cap<2;cap++) {
            int *v=profiles[profile];
            BlastInitialWordOptions opt={0};opt.window_size=v[0];opt.scan_range=v[1];opt.program_number=eBlastTypeBlastn;
            BlastUngappedCutoffs cuts[2]={{0}};
            for(int c=0;c<2;c++){cuts[c].x_dropoff_init=cuts[c].x_dropoff=c?7:20;cuts[c].cutoff_score=v[3]+c;cuts[c].reduced_nucl_cutoff_score=v[4]+c;}
            BlastInitialWordParameters params={0};params.options=&opt;params.cutoffs=cuts;params.ungapped_extension=v[2];params.container_type=hash?eDiagHash:eDiagArray;
            memcpy(params.nucl_score_table,table,sizeof(table));
            Blast_ExtendWord *ewp=NULL;BlastExtendWordNew(257,&params,&ewp);
            BlastInitHitList *hits=BLAST_InitHitListNew();BlastUngappedStats stats={0};
            int capacity=cap?4096+longest:MAX(1,longest+(wrap.lut_type==eMBLookupTable));
            BlastOffsetPair* pairs=malloc(capacity*sizeof(*pairs));
            printf("initial_new %d %d %d %d %d %d %d\n",hash,v[0],v[1],v[2],v[3],v[4],v[5]);
            puts("finder_stats_new");
            int lengths[]={0,1,7,8,9,11,27,28,29,63,64,65,127,128,129,257};
            for(int mode=0;mode<3;mode++)for(int n=0;n<16;n++) {
                int length=lengths[n];Uint1 data[69]={0};
                for(int i=0;i<length;i++){Uint1 base=query.sequence[i]&3;if(mode==1)base=(base+1)%4;if(mode==2&&i%11==5)base=(base+1)%4;data[i/4]|=base<<(6-2*(i%4));}
                data[length/4]|=length%4;
                BLAST_SequenceBlk subj={0};SSeqRange range={0,length};subj.sequence=data;subj.length=length;subj.seq_ranges=&range;subj.num_seq_ranges=1;subj.mask_type=eNoSubjMasking;
                BlastInitHitListReset(hits);
                if(mode==1&&n==5){if(hash)ewp->hash_table->offset=INT4_MAX/4;else ewp->diag_table->offset=INT4_MAX/4;puts("finder_rollover");}
                int with_stats=(mode+n)%2;
                int status=BlastNaWordFinder(&subj,&query,&info,&wrap,cells_ptr,&params,ewp,pairs,capacity,hits,with_stats?&stats:NULL);
                printf("finder_subject %d ",length);bytes(data,69);printf(" %d %d %d %d",capacity,with_stats,status,hits->total);
                for(int i=0;i<hits->total;i++){BlastInitHSP* h=hits->init_hsp_array+i;printf(" %d %d %d",h->offsets.qs_offsets.q_off,h->offsets.qs_offsets.s_off,h->ungapped_data!=NULL);if(h->ungapped_data)printf(" %d %d %d %d",h->ungapped_data->q_start,h->ungapped_data->s_start,h->ungapped_data->length,h->ungapped_data->score);}
                printf(" %" PRId64 " %d %d %d %d\n",stats.lookup_hits,stats.num_seqs_lookup_hits,stats.init_extends,stats.good_init_extends,stats.num_seqs_passed);
                printf("finder_offset %d\n",hash?ewp->hash_table->offset:ewp->diag_table->offset);
            if(n==5 || n==15) {
            if(hash) {
                BLAST_DiagHash* h=ewp->hash_table;printf("initial_hash %u %u",h->capacity,h->occupancy);
                for(int i=0;i<512;i++)printf(" %u",h->backbone[i]);
                for(unsigned i=0;i<h->occupancy;i++)printf(" %d %d %d %u %u",h->chain[i].diag,h->chain[i].level,h->chain[i].hit_len,h->chain[i].hit_saved,h->chain[i].next);
            } else {
                BLAST_DiagTable* d=ewp->diag_table;printf("initial_array %d",d->diag_array_length);
                for(int i=0;i<d->diag_array_length;i++)printf(" %d %u %d",d->hit_level_array[i].last_hit,d->hit_level_array[i].flag,d->hit_len_array?d->hit_len_array[i]:-1);
            }
            putchar('\n');
            }

            }
            free(pairs);BLAST_InitHitListFree(hits);BlastExtendWordFree(ewp);
        }
        if (wrap.lut_type==eSmallNaLookupTable) BlastSmallNaLookupTableDestruct(wrap.lut);
        else if (wrap.lut_type==eNaLookupTable) BlastNaLookupTableDestruct(wrap.lut);
        else BlastMBLookupTableDestruct(wrap.lut);
        free(query.compressed_nuc_seq_start);
    }
}

static void ranges(void)
{
    SSeqRange ranges[]={{0,5},{7,9},{30,60}};BLAST_SequenceBlk subject={0};subject.seq_ranges=ranges;subject.num_seq_ranges=3;
    for(int word=8;word<=28;word+=10)for(int start=0;start<70;start+=3)for(int end=0;end<70;end+=5){
        Int4 r[]={0,start,end};int ok=s_DetermineScanningOffsets(&subject,word,8,r);
        printf("finder_ranges %d 8 %d %d %d %d %d %d\n",word,start,end,ok,r[0],r[1],r[2]);
    }
}

static void stats_cases(void)
{
    for(int present=0;present<2;present++)for(int total=0;total<3;total++)for(int saved=-1;saved<=1;saved++) {
        BlastUngappedStats stats={INT64_MAX-1,INT32_MAX,INT32_MAX-1,INT32_MAX,INT32_MAX};
        Blast_UngappedStatsUpdate(present?&stats:NULL,total,3,saved);
        printf("finder_stats_case %d %d 3 %d %" PRId64 " %d %d %d %d\n",present,total,saved,stats.lookup_hits,stats.num_seqs_lookup_hits,stats.init_extends,stats.good_init_extends,stats.num_seqs_passed);
    }
}
int main(void) {stats_cases();ranges();finder_cases();return 0;}

