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
static void initial_cases(void)
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

        int profiles[][6]={{0,0,1,0,0,0},{0,0,1,50,25,0},{0,0,1,300,301,0},{40,0,1,100,80,0},{40,3,1,300,200,0},{40,3,0,0,0,0},{0,0,1,30,60,1},{0,0,0,0,0,0},{40,3,1,0,10,1}};
        for (int hash=0;hash<2;hash++) for(int profile=0;profile<9;profile++) {
            int *v=profiles[profile];
            BlastInitialWordOptions opt={0};opt.window_size=v[0];opt.scan_range=v[1];opt.program_number=eBlastTypeBlastn;
            BlastUngappedCutoffs cuts[2]={{0}};
            for(int c=0;c<2;c++){cuts[c].x_dropoff_init=cuts[c].x_dropoff=c?7:20;cuts[c].cutoff_score=v[3]+c;cuts[c].reduced_nucl_cutoff_score=v[4]+c;}
            BlastInitialWordParameters params={0};params.options=&opt;params.cutoffs=cuts;params.ungapped_extension=v[2];params.matrix_only_scoring=v[5];params.container_type=hash?eDiagHash:eDiagArray;
            memcpy(params.nucl_score_table,table,sizeof(table));
            Blast_ExtendWord *ewp=NULL;BlastExtendWordNew(257,&params,&ewp);
            BlastInitHitList *hits=BLAST_InitHitListNew();
            printf("initial_new %d %d %d %d %d %d %d\n",hash,v[0],v[1],v[2],v[3],v[4],v[5]);
            if(v[0]==40 && v[1]==3) for(int side=-1;side<=1;side+=2) {
                BlastExtendWordFree(ewp);BlastExtendWordNew(257,&params,&ewp);BlastInitHitListReset(hits);
                int qpos=60,soff=60,masks=0,range=soff+options.word_size;
                int level=soff+options.word_size+1;
                if(hash) s_BlastDiagHashInsert(ewp->hash_table,side,level,options.word_size,0,0,100000);
                else {int slot=(ewp->diag_table->diag_array_length+side)&ewp->diag_table->diag_mask;ewp->diag_table->hit_level_array[slot].last_hit=level;ewp->diag_table->hit_len_array[slot]=options.word_size;}
                printf("initial_prime %d %d %d\n",hash,side,level);
                int ready=hash?s_BlastnDiagHashExtendInitialHit(&query,&subject,qpos,soff,masked,&info,range,options.word_size,width,&wrap,&params,cells_ptr,ewp->hash_table,hits,masks)
                              :s_BlastnDiagTableExtendInitialHit(&query,&subject,qpos,soff,masked,&info,range,options.word_size,width,&wrap,&params,cells_ptr,ewp->diag_table,hits,masks);
                printf("initial_hit %d %d %d %d %d %d",qpos,soff,masks,range,ready,hits->total);
                if(hits->total){BlastInitHSP* h=hits->init_hsp_array+hits->total-1;printf(" %d %d %d",h->offsets.qs_offsets.q_off,h->offsets.qs_offsets.s_off,h->ungapped_data!=NULL);if(h->ungapped_data)printf(" %d %d %d %d",h->ungapped_data->q_start,h->ungapped_data->s_start,h->ungapped_data->length,h->ungapped_data->score);}
                putchar('\n');
            }
            BlastExtendWordFree(ewp);BlastExtendWordNew(257,&params,&ewp);BlastInitHitListReset(hits);
            puts("initial_restart");
            for(int pass=0;pass<2;pass++) for(int qpos=1;qpos<245;qpos+=7) for(int delta=-2;delta<=2;delta++) {
                int qend=qpos<129?128:257;int soff=qpos+delta;
                if(qpos+options.word_size>qend || soff<0 || soff+options.word_size>257)continue;
                int masks=(qpos+delta+pass)%2;int range=257;
                int ready=hash?s_BlastnDiagHashExtendInitialHit(&query,&subject,qpos,soff,masked,&info,range,options.word_size,width,&wrap,&params,cells_ptr,ewp->hash_table,hits,masks)
                              :s_BlastnDiagTableExtendInitialHit(&query,&subject,qpos,soff,masked,&info,range,options.word_size,width,&wrap,&params,cells_ptr,ewp->diag_table,hits,masks);
                printf("initial_hit %d %d %d %d %d %d",qpos,soff,masks,range,ready,hits->total);
                if(hits->total){BlastInitHSP* h=hits->init_hsp_array+hits->total-1;printf(" %d %d %d",h->offsets.qs_offsets.q_off,h->offsets.qs_offsets.s_off,h->ungapped_data!=NULL);if(h->ungapped_data)printf(" %d %d %d %d",h->ungapped_data->q_start,h->ungapped_data->s_start,h->ungapped_data->length,h->ungapped_data->score);}
                putchar('\n');
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
        if (wrap.lut_type==eSmallNaLookupTable) BlastSmallNaLookupTableDestruct(wrap.lut);
        else if (wrap.lut_type==eNaLookupTable) BlastNaLookupTableDestruct(wrap.lut);
        else BlastMBLookupTableDestruct(wrap.lut);
        free(query.compressed_nuc_seq_start);
    }
}
int main(void) {initial_cases();return 0;}
