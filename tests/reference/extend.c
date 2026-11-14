/* Executes the original extension-state and hit-list definitions unchanged. */
#include "blast_extend.c"
static void state(Blast_ExtendWord* e)
{
    if(e->diag_table) {
        BLAST_DiagTable* d=e->diag_table;
        printf("state %d %d %d %d %d",d->offset,d->diag_array_length,d->diag_mask,d->multiple_hits,d->actual_window);
        for(int i=0;i<d->diag_array_length;i++) printf(" %d %u %d",d->hit_level_array[i].last_hit,d->hit_level_array[i].flag,d->hit_len_array?d->hit_len_array[i]:-1);
    } else {
        BLAST_DiagHash* h=e->hash_table;
        printf("hash %d %u %u %u",h->offset,h->num_buckets,h->occupancy,h->capacity);
        for(unsigned i=0;i<h->num_buckets;i++) printf(" %u",h->backbone[i]);
        for(int i=0;i<3;i++) printf(" %d %d %u %d %u",h->chain[i].diag,h->chain[i].level,h->chain[i].hit_saved,h->chain[i].hit_len,h->chain[i].next);
    }
    putchar('\n');
}
static void list(BlastInitHitList* l)
{
    printf("list %d %d %d %d",l->total,l->allocated,l->do_not_reallocate,Blast_InitHitListIsSortedByScore(l));
    for(int i=0;i<l->total;i++) {
        BlastInitHSP* h=l->init_hsp_array+i;
        printf(" %d %d %d",h->offsets.qs_offsets.q_off,h->offsets.qs_offsets.s_off,h->ungapped_data!=NULL);
        if(h->ungapped_data) printf(" %d %d %d %d",h->ungapped_data->q_start,h->ungapped_data->s_start,h->ungapped_data->length,h->ungapped_data->score);
    }
    putchar('\n');
}
int main(void)
{
    int lengths[]={0,1,15,16,17,127};int windows[]={0,1,40};
    for(int hash=0;hash<2;hash++) for(int k=0;k<6;k++) for(int w=0;w<3;w++) {
        BlastInitialWordOptions opt={0};opt.window_size=windows[w];
        BlastInitialWordParameters p={0};p.options=&opt;p.container_type=hash?eDiagHash:eDiagArray;
        Blast_ExtendWord* e=NULL;ASSERT(BlastExtendWordNew(lengths[k],&p,&e)==0);
        printf("new %d %d %d\n",lengths[k],hash,windows[w]);state(e);
        for(int n=0;n<5;n++) {
            int offset=n==0?0:n==1?INT4_MAX/4-1:n==2?INT4_MAX/4:n==3?INT4_MAX/4+1:123;
            printf("dirty %d\n",offset);
            if(e->diag_table) {
                BLAST_DiagTable* d=e->diag_table;d->offset=offset;
                for(int i=0;i<d->diag_array_length;i++) {d->hit_level_array[i].last_hit=i-17;d->hit_level_array[i].flag=i%2;if(d->hit_len_array)d->hit_len_array[i]=(Uint1)(i+129);}
            } else {
                BLAST_DiagHash* h=e->hash_table;h->offset=offset;h->occupancy=3;
                for(int i=0;i<512;i++) h->backbone[i]=i%3;
                for(int i=0;i<3;i++) {h->chain[i].diag=-i;h->chain[i].level=i-13;h->chain[i].hit_saved=i%2;h->chain[i].hit_len=i+23;h->chain[i].next=i%3;}
            }
            printf("exit 77 %d\n",Blast_ExtendWordExit(e,77));state(e);
        }
        BlastExtendWordFree(e);
    }
    printf("null %d %d\n",Blast_ExtendWordExit(NULL,2),s_BlastDiagClear(NULL));
    BlastInitHitList *l=BLAST_InitHitListNew(),*dst=BLAST_InitHitListNew();puts("list_new");list(l);
    for(int i=0;i<501;i++) {
        if(i==100||i==200){l->do_not_reallocate=TRUE;puts("lock");}
        if(i==101||i==201){l->do_not_reallocate=FALSE;puts("unlock");}
        BlastUngappedData* d=NULL;
        if(i%7){d=calloc(1,sizeof(*d));d->q_start=i%17;d->s_start=i%13;d->length=1+i%31;d->score=i%19-9;}
        Boolean ok=BLAST_SaveInitialHit(l,i,-i,d);
        printf("save %d %d %d %d\n",i,-i,i%7!=0,ok);if(!ok) free(d);
        if(i%100==0)list(l);
    }
    puts("sort");Blast_InitHitListSortByScore(l);list(l);
    puts("move");BlastInitHitListMove(dst,l);list(l);list(dst);
    puts("reset");BlastInitHitListReset(dst);list(dst);
    puts("hsp");BlastSaveInitHsp(dst,3,5,7,9,11,13);list(dst);
    BLAST_InitHitListFree(l);BLAST_InitHitListFree(dst);
    return 0;
}
