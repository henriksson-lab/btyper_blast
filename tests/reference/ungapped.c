/* Original file-local nucleotide kernels are included unchanged. The test
 * driver supplies arrays/options, then unused toolkit entry points are stripped.
 */
#include <inttypes.h>
#include <strings.h>
#include "na_ungapped.c"
#include "blast_nalookup.c"
#include "blast_lookup.c"
#include "blast_util.c"

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
static void extension_cases(void)
{
    Int4 cells[16][16],*matrix[16];int table[256];
    int rewards[][2]={{1,-2},{2,-3},{2,-1},{5,-1},{0,-3}};
    int lengths[]={1,2,3,4,7,8,9,15,16,17,31,32,33,63,64,65,127,128,129,256};
    for (int profile=0;profile<5;++profile) {
        matrix_case(rewards[profile][0],rewards[profile][1],cells,matrix,table);
        for (int k=0;k<20;++k) for (int pattern=0;pattern<4;++pattern) {
            int qlen=lengths[k],slen=MAX(1,qlen+(pattern==1?-3:pattern==2?3:0));
            Uint1 qb[260]={0},packed[70]={0};qb[0]=15;qb[qlen+1]=15;
            for (int i=0;i<qlen;++i) {
                qb[i+1]=pattern==0?0:pattern==1?i%4:next()%4;
                if (pattern==3 && i%13==5) qb[i+1]=4+next()%11;
                if (pattern==3 && qlen>32 && i==qlen/2) qb[i+1]=15;
            }
            for (int i=0;i<slen;++i) {
                Uint1 base=qb[1+i%qlen]&3;
                if (pattern==2 && i%11==3) base=(base+1)%4;
                if (pattern==3 && i%7==2) base=next()%4;
                packed[i/4] |= base << (6-2*(i%4));
            }
            packed[slen/4]|=slen%4;
            BLAST_SequenceBlk query={0},subject={0};query.sequence=qb+1;query.length=qlen;
            subject.sequence=packed;subject.length=slen;
            printf("sequences ");bytes(query.sequence,qlen);printf(" %d ",slen);bytes(packed,slen/4+4);putchar('\n');
            int starts[]={0,1,2,3,qlen/2,qlen-1};int drops[]={0,-1,-7,-20};
            for (int a=0;a<6;++a) for (int delta=-1;delta<=1;++delta) {
                int qoff=starts[a],soff=qoff+delta;
                if (qoff>=qlen || soff<0 || soff>=slen) continue;
                for (int d=0;d<4;++d) {
                    int x=drops[d];BlastUngappedData data={-997,-998,-999,-1000};
                    s_NuclUngappedExtendExact(&query,&subject,matrix,qoff,soff,x,&data);
                    printf("exact %d %d %d",qoff,soff,x);result(&data);
                    if (qoff+4>qlen || soff+4>slen) continue;
                    int match_end=MIN(slen,soff+11);
                    BlastUngappedData baseline;
                    s_NuclUngappedExtend(&query,&subject,matrix,qoff,match_end,soff,x,&baseline,table,1000000);
                    int cutoffs[]={0,1,20,baseline.score,baseline.score+1};
                    for (int t=0;t<5;++t) {
                        data=(BlastUngappedData){-997,-998,-999,-1000};
                        s_NuclUngappedExtend(&query,&subject,matrix,qoff,match_end,soff,x,&data,table,cutoffs[t]);
                        printf("approx %d %d %d %d %d",qoff,match_end,soff,x,cutoffs[t]);result(&data);
                    }
                }
            }
        }
    }
}
static void word_cases(void)
{
    int kinds[][3]={{eSmallNaLookupTable,11,8},{eSmallNaLookupTable,28,8},
        {eNaLookupTable,11,8},{eNaLookupTable,28,8},{eMBLookupTable,11,10},
        {eMBLookupTable,11,11},{eMBLookupTable,28,11},{eMBLookupTable,28,12}};
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
        for (int qpos=0;qpos+width<=257;++qpos) {
            int index=ComputeTableIndex(width,2,query.sequence+qpos);
            int positions[]={-1,0,qpos,qpos+1,264};int indices[]={index,index|(int)0x80000000u,index^(1<<28)};
            for (int a=0;a<5;++a) for (int b=0;b<3;++b)
                printf("member %d %d %d\n",indices[b],positions[a],((T_Lookup_Callback)wrap.lookup_callback)(&wrap,indices[b],positions[a]));
            printf("masked %d %d %d\n",qpos,qpos,s_IsSeedMasked(&wrap,&subject,qpos,width,qpos));
        }
        int starts[]={0,1,2,3,4,5,10,24,48,80,100,129,130,131,132,140,160,180,229};
        for (int a=0;a<19;++a) for (int delta=-1;delta<=1;++delta) {
            int qoff=starts[a],soff=qoff+delta;
            int qend=qoff<129?128:257;
            if (qoff+options.word_size>qend || soff<0 || soff+options.word_size>257) continue;
            for (int loc=0;loc<2;++loc) for (int dbl=0;dbl<2;++dbl) {
                int q=qoff,s=soff,extended=-999;
                int type=s_TypeOfWord(&query,&subject,&q,&s,loc?masked:NULL,&info,257,options.word_size,width,&wrap,dbl,&extended);
                printf("word %d %d %d %d %d %d %d %d %d\n",qoff,soff,loc,dbl,type,q,s,extended,masked!=NULL);
            }
        }
        if (wrap.lut_type==eSmallNaLookupTable) BlastSmallNaLookupTableDestruct(wrap.lut);
        else if (wrap.lut_type==eNaLookupTable) BlastNaLookupTableDestruct(wrap.lut);
        else BlastMBLookupTableDestruct(wrap.lut);
        free(query.compressed_nuc_seq_start);
    }
}
static void hash_snapshot(const BLAST_DiagHash* table)
{
    unsigned used=0;for (unsigned i=0;i<DIAGHASH_NUM_BUCKETS;++i) if (table->backbone[i]) ++used;
    printf("hash_snapshot %u %u %u",table->capacity,table->occupancy,used);
    for (unsigned i=0;i<DIAGHASH_NUM_BUCKETS;++i) if (table->backbone[i]) printf(" %u %u",i,table->backbone[i]);
    for (unsigned i=0;i<table->occupancy;++i) {
        const DiagHashCell* c=&table->chain[i];
        printf(" %d %d %d %d %u",c->diag,(int)c->level,c->hit_len,(int)c->hit_saved,c->next);
    }
    putchar('\n');
}
static void hash_probe(BLAST_DiagHash* table,int diag)
{
    Int4 level=-707,len=-708,saved=-709;
    int found=s_BlastDiagHashRetrieve(table,diag,&level,&len,&saved);
    printf("hash_get %d %d %d %d %d\n",diag,found,level,len,saved);
}
static void hash_insert(BLAST_DiagHash* table,int diag,int level,int length,int flag,int off,int window)
{
    int ok=s_BlastDiagHashInsert(table,diag,level,length,flag,off,window);
    printf("hash_put %d %d %d %d %d %d %d %u %u\n",diag,level,length,flag,off,window,ok,table->capacity,table->occupancy);
    hash_probe(table,diag);hash_probe(table,INT32_MAX);
}
static void hash_cases(void)
{
    BLAST_DiagHash table={0};table.num_buckets=DIAGHASH_NUM_BUCKETS;table.occupancy=table.capacity=1;
    table.backbone=calloc(DIAGHASH_NUM_BUCKETS,sizeof(Uint4));table.chain=calloc(1,sizeof(DiagHashCell));
    printf("hash_new %u %u %u\n",table.num_buckets,table.occupancy,table.capacity);
    hash_probe(&table,0);hash_snapshot(&table);
    for (int i=0;i<600;++i) {
        hash_insert(&table,i*512,1000+i,11+i%28,i%4,0,100000);
        if (i%31==0 || i==599) hash_snapshot(&table);
    }
    for (int i=0;i<200;++i) {
        hash_insert(&table,i*512,5000+i,28,i%2,0,100000);
        if (i%31==0) hash_snapshot(&table);
    }
    for (int i=0;i<200;++i) {
        hash_insert(&table,(1000+i)*512,100000+i,17,i%2,100000+i,40);
        if (i%31==0) hash_snapshot(&table);
    }
    for (int i=0;i<400;++i) {
        int diag=(int)(next()*257+next())-2000;
        hash_insert(&table,diag,200000+i,11+i%45,i%3,200000+i,i%2?40:0);
        if (i%31==0) hash_snapshot(&table);
    }
    int levels[]={INT32_MIN,-1073741825,-1073741824,-1,0,1073741823,1073741824,INT32_MAX};
    for (int i=0;i<8;++i) {
        hash_insert(&table,-(i+1)*512,levels[i],i,-1+i,0,100000);
        hash_snapshot(&table);
    }
    free(table.backbone);free(table.chain);
}
int main(void) {extension_cases();word_cases();hash_cases();return 0;}
