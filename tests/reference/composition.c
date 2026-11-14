/* Original composition and score-frequency definitions execute unchanged. */
#include <inttypes.h>
#include "blast_stat.c"
static void value(double x){uint64_t b;memcpy(&b,&x,8);printf(" %016" PRIx64,b);}
static void bytes(const Uint1* p,int n){if(!n)putchar('-');for(int i=0;i<n;i++)printf("%02x",p[i]);}
static void freq(const Blast_ResFreq* p,int n){for(int i=0;i<n;i++)value(p->prob0[i]);putchar('\n');}
static void block(const Blast_KarlinBlk* k){value(k->Lambda);value(k->K);value(k->logK);value(k->H);value(k->paramC);putchar('\n');}
static void sf(const Blast_ScoreFreq* s){printf(" %d %d",s->obs_min,s->obs_max);value(s->score_avg);for(int i=s->score_min;i<=s->score_max;i++)value(s->sprob[i]);putchar('\n');}
static void alphabet(int protein,int reward,int penalty,int ambigmode)
{
    int size=protein?28:16;Int4 cells[28][28],*rows[28];
    for(int i=0;i<28;i++){rows[i]=cells[i];for(int j=0;j<28;j++)cells[i][j]=i==j?5:-4;}
    BlastScoreBlk sb={0};SBlastScoreMatrix matrix={0};matrix.data=rows;matrix.nrows=matrix.ncols=size;
    sb.alphabet_size=size;sb.alphabet_start=0;sb.alphabet_code=protein?BLASTAA_SEQ_CODE:BLASTNA_SEQ_CODE;sb.protein_alphabet=protein;sb.matrix=&matrix;
    if(!protein){sb.reward=reward;sb.penalty=penalty;BlastScoreBlkNuclMatrixCreate(&sb);}
    Uint1 ambigs[4]={protein?21:14,0,protein?25:15,protein?3:7};sb.ambiguous_res=ambigs;sb.ambig_occupy=ambigmode;
    printf("alphabet %d %d %d",sb.alphabet_code,size,ambigmode);for(int i=0;i<ambigmode;i++)printf(" %u",ambigs[i]);for(int i=0;i<size;i++)for(int j=0;j<size;j++)printf(" %d",cells[i][j]);putchar('\n');
    BlastScoreBlkMaxScoreSet(&sb);printf("bounds %d %d\n",sb.loscore,sb.hiscore);
    Blast_ResFreq *p=Blast_ResFreqNew(&sb),*standard=Blast_ResFreqNew(&sb);Blast_ResComp *comp=BlastResCompNew(&sb);
    printf("new_frequency %u",p->alphabet_code);freq(p,size);printf("new_comp %u",comp->alphabet_code);for(int i=0;i<size;i++)printf(" %d",comp->comp0[i]);putchar('\n');
    Blast_ResFreqStdComp(&sb,standard);printf("standard");freq(standard,size);
    for(int dirty=0;dirty<2;dirty++){for(int i=0;i<size;i++)p->prob[i]=dirty?.1*(i+1):0;int status=Blast_ResFreqStdComp(&sb,p);printf("standard_dirty %d %d",dirty,status);freq(p,size);}
    Blast_ScoreBlkKbpIdealCalc(&sb);printf("ideal");block(sb.kbp_ideal);sb.kbp_ideal=Blast_KarlinBlkFree(sb.kbp_ideal);
    {int lo=sb.loscore,hi=sb.hiscore;sb.loscore=sb.hiscore=0;int status=Blast_ScoreBlkKbpIdealCalc(&sb);printf("ideal_invalid %d",status);block(sb.kbp_ideal);sb.kbp_ideal=Blast_KarlinBlkFree(sb.kbp_ideal);sb.loscore=lo;sb.hiscore=hi;}
    int lengths[]={0,1,2,3,15,16,17,27,28,29,127,128,129,1024};
    for(int k=0;k<14;k++)for(int mode=0;mode<6;mode++){
        int length=lengths[k];Uint1 sequence[1024];for(int i=0;i<length;i++)sequence[i]=mode==0?0:mode==1?size-1:mode==2?i%size:mode==3?(i*7+11)%size:mode==4?(protein?21:14):(protein?i%size:(i*13+17)%256);
        for(int i=0;i<size;i++)comp->comp[i]=991;
        int status=BlastResCompStr(&sb,comp,(char*)sequence,length);
        printf("sequence ");bytes(sequence,length);printf(" %d",status);for(int i=0;i<size;i++)printf(" %d",comp->comp[i]);putchar('\n');
        for(int i=0;i<size;i++)p->prob[i]=991;
        status=Blast_ResFreqResComp(&sb,p,comp);printf("composition_freq %d",status);freq(p,size);
        for(int i=0;i<size;i++)p->prob[i]=991;
        status=Blast_ResFreqString(&sb,p,(char*)sequence,length);printf("string_freq %d",status);freq(p,size);
        Blast_ScoreFreq* s=Blast_ScoreFreqNew(sb.loscore,sb.hiscore);status=BlastScoreFreqCalc(&sb,s,p,standard);printf("score_freq %d",status);sf(s);
        Blast_KarlinBlk out={91,92,93,94,95};status=Blast_KarlinBlkUngappedCalc(&out,s);printf("karlin %d",status);block(&out);
        Blast_ScoreFreqFree(s);
    }
    double scales[]={0,0.001,0.01,0.1,1.0};
    for(int mode=0;mode<5;mode++){
        for(int i=0;i<size;i++)p->prob[i]=standard->prob[i]*scales[mode];
        Blast_ScoreFreq* s=Blast_ScoreFreqNew(sb.loscore,sb.hiscore);int status=BlastScoreFreqCalc(&sb,s,p,p);
        printf("scaled_score");value(scales[mode]);printf(" %d",status);sf(s);Blast_ScoreFreqFree(s);
    }
    for(int missing=0;missing<4;missing++){
        Blast_ScoreFreq* s=Blast_ScoreFreqNew(sb.loscore+1,sb.hiscore);s->obs_min=-111;s->obs_max=222;s->score_avg=333;for(int i=s->score_min;i<=s->score_max;i++)s->sprob[i]=444;
        int status=BlastScoreFreqCalc(missing&1?NULL:&sb,missing&2?NULL:s,p,standard);
        printf("score_error %d %d",missing,status);sf(s);Blast_ScoreFreqFree(s);
    }
    double norms[]={0,1,-1,2,NAN,INFINITY};
    for(int mode=0;mode<5;mode++)for(int n=0;n<6;n++) {
        for(int i=0;i<size;i++)p->prob[i]=mode==0?0:mode==1?(i+1):mode==2?(i==3?-1:i):mode==3?(i==3?NAN:.1):.1;
        int status=Blast_ResFreqNormalize(&sb,p,norms[n]);printf("normalize %d",mode);value(norms[n]);printf(" %d",status);freq(p,size);
    }
    for(int missing=0;missing<8;missing++) {
        for(int i=0;i<size;i++)comp->comp[i]=991;
        Uint1 seq[]={0};int status=BlastResCompStr(missing&1?NULL:&sb,missing&2?NULL:comp,missing&4?NULL:(char*)seq,1);
        printf("comp_null %d %d",missing,status);for(int i=0;i<size;i++)printf(" %d",comp->comp[i]);putchar('\n');
    }
    comp->alphabet_code=0;for(int i=0;i<size;i++)comp->comp[i]=991;Uint1 seq[]={0};int status=BlastResCompStr(&sb,comp,(char*)seq,1);printf("comp_mismatch %d",status);for(int i=0;i<size;i++)printf(" %d",comp->comp[i]);putchar('\n');
    for(int i=0;i<size;i++)p->prob[i]=991;status=Blast_ResFreqResComp(&sb,p,comp);printf("freq_mismatch %d",status);freq(p,size);
    for(int missing=1;missing<4;missing++){status=Blast_ResFreqResComp(&sb,missing&1?NULL:p,missing&2?NULL:comp);printf("freq_null %d %d",missing,status);freq(p,size);}
    BlastResCompDestruct(comp);Blast_ResFreqFree(p);Blast_ResFreqFree(standard);
}
int main(void)
{
    printf("null_frequency %d\n",Blast_ResFreqNew(NULL)!=NULL);printf("ideal_null %d\n",Blast_ScoreBlkKbpIdealCalc(NULL));
    for(int code=0;code<3;code++)for(int len=0;len<=21;len++){Uint1 out[21];memset(out,173,21);int alpha=code==0?BLASTAA_SEQ_CODE:code==1?BLASTNA_SEQ_CODE:0;int status=Blast_GetStdAlphabet(alpha,out,len);printf("std_alphabet %d %d %d ",alpha,len,status);bytes(out,21);putchar('\n');}
    for(int i=0;i<20;i++){printf("background %u",Robinson_prob[i].ch);value(Robinson_prob[i].p);putchar('\n');}
    for(int i=0;i<4;i++){printf("nt_background %u",nt_prob[i].ch);value(nt_prob[i].p);putchar('\n');}
    for(int protein=0;protein<2;protein++)for(int amb=0;amb<3;amb++)alphabet(protein,1,-2,amb);
    alphabet(0,2,-3,1);
    Int4 bounds[][4]={{BLAST_SCORE_MIN,BLAST_SCORE_MAX,INT32_MIN,INT32_MAX},{-1,1,BLAST_SCORE_MIN,BLAST_SCORE_MAX},{BLAST_SCORE_MIN+1,BLAST_SCORE_MAX-1,BLAST_SCORE_MIN,BLAST_SCORE_MAX}};
    for(int k=0;k<3;k++){Int4* rows[]={bounds[k],bounds[k]+2};SBlastScoreMatrix m={0};m.data=rows;BlastScoreBlk sb={0};sb.matrix=&m;sb.alphabet_size=2;BlastScoreBlkMaxScoreSet(&sb);printf("matrix_bounds");for(int i=0;i<4;i++)printf(" %d",bounds[k][i]);printf(" %d %d\n",sb.loscore,sb.hiscore);}
    return 0;
}
