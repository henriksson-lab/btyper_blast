#include <inttypes.h>
#include "blast_stat.c"
#include "blast_parameters.c"
int main(void){
 for(int profile=0;profile<2;profile++){
  Int4 m[16][16],*rows[16];for(int i=0;i<16;i++)rows[i]=m[i];SBlastScoreMatrix matrix={0};matrix.data=rows;matrix.nrows=matrix.ncols=16;
  BlastScoreBlk sb={0};sb.matrix=&matrix;sb.alphabet_code=BLASTNA_SEQ_CODE;sb.alphabet_size=16;sb.reward=profile?2:1;sb.penalty=profile?-3:-2;sb.scale_factor=1;
  Uint1 ambiguous=14;sb.ambiguous_res=&ambiguous;sb.ambig_occupy=1;BlastScoreBlkNuclMatrixCreate(&sb);BlastScoreBlkMaxScoreSet(&sb);
  BlastContextInfo contexts[2]={{0}};contexts[0].query_length=32;contexts[1].query_offset=33;contexts[1].query_length=32;contexts[0].is_valid=contexts[1].is_valid=TRUE;
  BlastQueryInfo info={0};info.contexts=contexts;info.last_context=1;info.num_queries=1;
  Uint1 query[65];for(int i=0;i<32;i++)query[i]=i%4;query[32]=15;memset(query+33,14,32);
  sb.sfp=calloc(2,sizeof(*sb.sfp));sb.kbp_std=calloc(2,sizeof(*sb.kbp_std));sb.kbp_psi=calloc(2,sizeof(*sb.kbp_psi));
  int status=Blast_ScoreBlkKbpUngappedCalc(eBlastTypeBlastn,&sb,query,&info,NULL);
  Blast_KarlinBlk gap={0};Boolean round=FALSE;int gs=Blast_KarlinBlkNuclGappedCalc(&gap,profile?5:0,profile?2:0,sb.reward,sb.penalty,sb.kbp_std[0],&round,NULL);
  Blast_KarlinBlk* gaps[]={&gap,NULL};sb.kbp_gap=gaps;
  BlastInitialWordOptions *opts=NULL;BlastInitialWordOptionsNew(eBlastTypeBlastn,&opts);
  BlastGappedCutoffs cuts[2]={{0}};cuts[0].cutoff_score_max=23;cuts[1].cutoff_score_max=INT4_MAX;BlastHitSavingParameters hits={0};hits.cutoffs=cuts;
  BlastInitialWordParameters* p=NULL;int ps=BlastInitialWordParametersNew(eBlastTypeBlastn,opts,&hits,NULL,&sb,&info,1000,&p);
  printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",profile,status,contexts[0].is_valid,contexts[1].is_valid,gs,ps,round,p->cutoff_score_min,p->x_dropoff_max,p->cutoffs[0].cutoff_score,p->cutoffs[0].x_dropoff_init,p->cutoffs[0].reduced_nucl_cutoff_score,p->cutoffs[1].cutoff_score,p->container_type);
 }
}
