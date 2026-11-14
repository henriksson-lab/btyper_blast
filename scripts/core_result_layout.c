/* NCBI BLAST+ 2.17.0+ ABI reporter; no BLAST implementation is linked/modified. */
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/blast_extend.h>
#define TYPE(T) printf("\"" #T "\":{\"size\":%zu,\"fields\":{",sizeof(T))
#define FIELD(T,F) printf("\"" #F "\":{\"offset\":%zu,\"size\":%zu}",offsetof(T,F),sizeof(((T*)0)->F))
#define NEXT(T,F) putchar(',');FIELD(T,F)
#define END() printf("}}")
int main(void){
 uint16_t byte_order=1;
 printf("{\"pointer_size\":%zu,\"double_size\":%zu,\"byte_order\":\"%s\",\"op_size\":%zu,\"types\":{",sizeof(void*),sizeof(double),*(unsigned char*)&byte_order?"little":"big",sizeof(EGapAlignOpType));
 TYPE(BlastHSPResults);FIELD(BlastHSPResults,num_queries);NEXT(BlastHSPResults,hitlist_array);END();putchar(',');
 TYPE(BlastHitList);FIELD(BlastHitList,hsplist_count);NEXT(BlastHitList,hsplist_max);NEXT(BlastHitList,worst_evalue);NEXT(BlastHitList,low_score);NEXT(BlastHitList,heapified);NEXT(BlastHitList,hsplist_array);NEXT(BlastHitList,hsplist_current);NEXT(BlastHitList,num_hits);END();putchar(',');
 TYPE(BlastHSPList);FIELD(BlastHSPList,oid);NEXT(BlastHSPList,query_index);NEXT(BlastHSPList,hsp_array);NEXT(BlastHSPList,hspcnt);NEXT(BlastHSPList,allocated);NEXT(BlastHSPList,hsp_max);NEXT(BlastHSPList,do_not_reallocate);NEXT(BlastHSPList,best_evalue);END();putchar(',');
 TYPE(BlastHSP);FIELD(BlastHSP,score);NEXT(BlastHSP,num_ident);NEXT(BlastHSP,num_positives);NEXT(BlastHSP,bit_score);NEXT(BlastHSP,evalue);NEXT(BlastHSP,query);NEXT(BlastHSP,subject);NEXT(BlastHSP,context);NEXT(BlastHSP,gap_info);NEXT(BlastHSP,num);NEXT(BlastHSP,comp_adjustment_method);NEXT(BlastHSP,pat_info);NEXT(BlastHSP,map_info);END();putchar(',');
 TYPE(BlastSeg);FIELD(BlastSeg,frame);NEXT(BlastSeg,offset);NEXT(BlastSeg,end);NEXT(BlastSeg,gapped_start);END();putchar(',');
 TYPE(GapEditScript);FIELD(GapEditScript,op_type);NEXT(GapEditScript,num);NEXT(GapEditScript,size);END();
 putchar(',');TYPE(BLAST_SequenceBlk);FIELD(BLAST_SequenceBlk,oid);NEXT(BLAST_SequenceBlk,length);NEXT(BLAST_SequenceBlk,frame);END();
 putchar(',');TYPE(BlastInitHitList);FIELD(BlastInitHitList,total);NEXT(BlastInitHitList,allocated);NEXT(BlastInitHitList,init_hsp_array);END();
 putchar(',');TYPE(BlastInitHSP);FIELD(BlastInitHSP,offsets);NEXT(BlastInitHSP,ungapped_data);END();
 putchar(',');TYPE(BlastOffsetPair);FIELD(BlastOffsetPair,qs_offsets.q_off);NEXT(BlastOffsetPair,qs_offsets.s_off);END();
 putchar(',');TYPE(BlastUngappedData);FIELD(BlastUngappedData,q_start);NEXT(BlastUngappedData,s_start);NEXT(BlastUngappedData,length);NEXT(BlastUngappedData,score);END();
 putchar(',');TYPE(BlastInitialWordParameters);FIELD(BlastInitialWordParameters,options);NEXT(BlastInitialWordParameters,x_dropoff_max);NEXT(BlastInitialWordParameters,cutoff_score_min);NEXT(BlastInitialWordParameters,cutoffs);NEXT(BlastInitialWordParameters,container_type);END();
 putchar(',');TYPE(BlastInitialWordOptions);FIELD(BlastInitialWordOptions,gap_trigger);NEXT(BlastInitialWordOptions,window_size);NEXT(BlastInitialWordOptions,scan_range);NEXT(BlastInitialWordOptions,x_dropoff);END();
 putchar(',');TYPE(BlastUngappedCutoffs);FIELD(BlastUngappedCutoffs,x_dropoff_init);NEXT(BlastUngappedCutoffs,x_dropoff);NEXT(BlastUngappedCutoffs,cutoff_score);NEXT(BlastUngappedCutoffs,reduced_nucl_cutoff_score);END();
 putchar(',');TYPE(BlastQueryInfo);FIELD(BlastQueryInfo,num_queries);NEXT(BlastQueryInfo,first_context);NEXT(BlastQueryInfo,last_context);NEXT(BlastQueryInfo,contexts);END();
 putchar(',');TYPE(BlastContextInfo);FIELD(BlastContextInfo,query_offset);NEXT(BlastContextInfo,query_length);NEXT(BlastContextInfo,frame);NEXT(BlastContextInfo,is_valid);END();
 putchar(',');TYPE(BlastExtensionParameters);FIELD(BlastExtensionParameters,gap_x_dropoff);NEXT(BlastExtensionParameters,gap_x_dropoff_final);END();
 printf("}}\n");return 0;
}
