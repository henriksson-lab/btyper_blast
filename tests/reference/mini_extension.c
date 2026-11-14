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
static void mini_cases(void)
{
    Int4 cells[16][16],*cells_ptr[16];int table[256];matrix_case(1,-2,cells,cells_ptr,table);
    int kinds[][3]={{eSmallNaLookupTable,11,8},{eSmallNaLookupTable,28,8},
        {eNaLookupTable,11,8},{eNaLookupTable,28,8},{eMBLookupTable,11,10},
        {eMBLookupTable,11,11},{eMBLookupTable,28,11},{eMBLookupTable,28,12},{eSmallNaLookupTable,8,8},{eSmallNaLookupTable,9,8}};
    for (int k=0;k<10;++k) for (int pattern=0;pattern<4;++pattern) {
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
        TNaScanSubjectFunction scan;
        MiniExtendFunction selected;
        int longest;
        if(wrap.lut_type==eSmallNaLookupTable){BlastSmallNaLookupTable* t=wrap.lut;scan=(TNaScanSubjectFunction)t->scansub_callback;selected=(MiniExtendFunction)t->extend_callback;longest=t->longest_chain;}
        else if(wrap.lut_type==eNaLookupTable){BlastNaLookupTable* t=wrap.lut;scan=(TNaScanSubjectFunction)t->scansub_callback;selected=(MiniExtendFunction)t->extend_callback;longest=t->longest_chain;}
        else{BlastMBLookupTable* t=wrap.lut;scan=(TNaScanSubjectFunction)t->scansub_callback;selected=(MiniExtendFunction)t->extend_callback;longest=t->longest_chain;}
        int capacity=4096+longest,total=0,allocated=1024;BlastOffsetPair *pairs=malloc(capacity*sizeof(*pairs)),*seeds=malloc(allocated*sizeof(*seeds));
        int offsets[2]={0,257-width};
        while(offsets[0]<=offsets[1]){
            int before=offsets[0];int count=scan(&wrap,&subject,pairs,capacity,offsets);ASSERT(offsets[0]>before);
            if(total+count>allocated){allocated=2*(total+count);seeds=realloc(seeds,allocated*sizeof(*seeds));}
            memcpy(seeds+total,pairs,count*sizeof(*pairs));total+=count;
        }
        printf("mini_seeds %d",total);for(int i=0;i<total;i++)printf(" %d %d",seeds[i].qs_offsets.q_off,seeds[i].qs_offsets.s_off);putchar('\n');
        MiniExtendFunction callbacks[]={s_BlastNaExtendDirect,s_BlastNaExtend,s_BlastNaExtendAligned,s_BlastSmallNaExtendAlignedOneByte,s_BlastSmallNaExtend};
        int selected_id=-1;for(int i=0;i<5;i++)if(selected==callbacks[i])selected_id=i;ASSERT(selected_id>=0);
        int profiles[][6]={{0,0,1,50,25,0},{40,3,1,300,200,0},{0,0,0,0,0,0}};
        for(int mode=0;mode<2;mode++)for (int hash=0;hash<2;hash++) for(int profile=0;profile<3;profile++) {
            int *v=profiles[profile];
            BlastInitialWordOptions opt={0};opt.window_size=v[0];opt.scan_range=v[1];opt.program_number=eBlastTypeBlastn;
            BlastUngappedCutoffs cuts[2]={{0}};
            for(int c=0;c<2;c++){cuts[c].x_dropoff_init=cuts[c].x_dropoff=c?7:20;cuts[c].cutoff_score=v[3]+c;cuts[c].reduced_nucl_cutoff_score=v[4]+c;}
            BlastInitialWordParameters params={0};params.options=&opt;params.cutoffs=cuts;params.ungapped_extension=v[2];params.container_type=hash?eDiagHash:eDiagArray;
            memcpy(params.nucl_score_table,table,sizeof(table));
            Blast_ExtendWord *ewp=NULL;BlastExtendWordNew(257,&params,&ewp);
            BlastInitHitList *hits=BLAST_InitHitListNew();
            printf("initial_new %d %d %d %d %d %d %d\n",hash,v[0],v[1],v[2],v[3],v[4],v[5]);
            int id=mode==0?selected_id:wrap.lut_type==eSmallNaLookupTable?4:1;
            printf("mini_mode %d %d\n",id,selected_id);
            for(int start=0;start<total;){
                int count=MIN(start%3==0?1:start%3==1?7:127,total-start);int previous=hits->total;
                int ready=callbacks[id](seeds+start,count,&params,&wrap,&query,&subject,cells_ptr,&info,ewp,hits,257);
                printf("mini_batch %d %d %d %d",start,count,ready,hits->total);
                for(int i=previous;i<hits->total;i++){BlastInitHSP* h=hits->init_hsp_array+i;printf(" %d %d %d",h->offsets.qs_offsets.q_off,h->offsets.qs_offsets.s_off,h->ungapped_data!=NULL);if(h->ungapped_data)printf(" %d %d %d %d",h->ungapped_data->q_start,h->ungapped_data->s_start,h->ungapped_data->length,h->ungapped_data->score);}
                putchar('\n');start+=count;
            }
            if(hash) {
                BLAST_DiagHash* h=ewp->hash_table;printf("initial_hash %u %u",h->capacity,h->occupancy);
                for(int i=0;i<512;i++)printf(" %u",h->backbone[i]);
                for(unsigned i=0;i<h->occupancy;i++)printf(" %d %d %d %u %u",h->chain[i].diag,h->chain[i].level,h->chain[i].hit_len,h->chain[i].hit_saved,h->chain[i].next);
            } else {
                BLAST_DiagTable* d=ewp->diag_table;printf("initial_array %d",d->diag_array_length);
                for(int i=0;i<d->diag_array_length;i++)printf(" %d %u %d",d->hit_level_array[i].last_hit,d->hit_level_array[i].flag,d->hit_len_array?d->hit_len_array[i]:-1);
            }
            putchar('\n');
            BLAST_InitHitListFree(hits);BlastExtendWordFree(ewp);
        }
        free(pairs);free(seeds);
        if (wrap.lut_type==eSmallNaLookupTable) BlastSmallNaLookupTableDestruct(wrap.lut);
        else if (wrap.lut_type==eNaLookupTable) BlastNaLookupTableDestruct(wrap.lut);
        else BlastMBLookupTableDestruct(wrap.lut);
        free(query.compressed_nuc_seq_start);
    }
}
int main(void) {printf("mini_tables");for(int i=0;i<256;i++)printf(" %u %u",s_ExactMatchExtendLeft[i],s_ExactMatchExtendRight[i]);putchar('\n');mini_cases();return 0;}
