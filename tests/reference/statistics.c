/* Selected file-local statistical routines execute unchanged. */
#include <inttypes.h>
#include "blast_stat.c"
static void value(double v) {uint64_t b;memcpy(&b,&v,sizeof(b));printf(" %016" PRIx64,b);}
static void block(const Blast_KarlinBlk* k){value(k->Lambda);value(k->K);value(k->logK);value(k->H);value(k->paramC);putchar('\n');}
static void frequency(int low,int high,int mode)
{
    Blast_ScoreFreq* s=Blast_ScoreFreqNew(low,high);ASSERT(s);s->obs_min=low;s->obs_max=high;
    if(mode<4) {
        double p=mode==0?.25:mode==1?.1:mode==2?.7:.5;
        s->sprob[low]=1-p;s->sprob[high]=p;
    } else {
        double sum=0;
        for(int i=low;i<=high;i++) {double p=(double)((i-low+mode)%7+1);s->sprob[i]=p;sum+=p;}
        for(int i=low;i<=high;i++)s->sprob[i]/=sum;
    }
    if(mode==7)s->sprob[0]=NAN;
    for(int i=low;i<=high;i++)s->score_avg+=i*s->sprob[i];
    printf("frequency %d %d",low,high);value(s->score_avg);for(int i=low;i<=high;i++)value(s->sprob[i]);putchar('\n');
    double guesses[]={0.5,0.0,-1.0,1.0,100.0,INFINITY,NAN};
    for(int i=0;i<7;i++){printf("lambda");value(guesses[i]);value(Blast_KarlinLambdaNR(s,guesses[i]));putchar('\n');}
    double lambda=Blast_KarlinLambdaNR(s,.5),h=BlastKarlinLtoH(s,lambda);
    double inputs[]={lambda,-1.0,0.0,0.1,0.5,1.0,1000.0,NAN};
    for(int i=0;i<8;i++){printf("h");value(inputs[i]);value(BlastKarlinLtoH(s,inputs[i]));putchar('\n');}
    double hs[]={h,-1.0,0.0,0.1,NAN};
    for(int i=0;i<8;i++)for(int j=0;j<5;j++) {
        printf("k");value(inputs[i]);value(hs[j]);value(BlastKarlinLHtoK(s,inputs[i],hs[j]));putchar('\n');
    }
    int d=-low;for(int i=1;i<=high-low&&d>1;i++)if(s->sprob[low+i]!=0)d=BLAST_Gcd(d,i);
    int itmaxes[]={0,1,5,20};int newtons[]={0,2,37};
    for(int i=0;i<4;i++)for(int j=0;j<3;j++) {
        int itn=-101;double result=NlmKarlinLambdaNR(s->sprob,d,low,high,.5,1e-5,itmaxes[i],newtons[j],&itn);
        printf("nr %d %d %d %d",d,itmaxes[i],newtons[j],itn);value(result);putchar('\n');
    }
    Blast_KarlinBlk k={91,92,93,94,95};int status=Blast_KarlinBlkUngappedCalc(&k,s);printf("block %d",status);block(&k);
    printf("null %d %d %d\n",Blast_KarlinBlkUngappedCalc(NULL,s),Blast_KarlinBlkUngappedCalc(&k,NULL),Blast_KarlinBlkUngappedCalc(NULL,NULL));
    Blast_ScoreFreqFree(s);
}
int main(void)
{
    EBlastProgramType programs[]={eBlastTypeBlastn,eBlastTypeTblastn};
    for(int i=0;i<2;i++){BlastInitialWordOptions* opt=NULL;ASSERT(BlastInitialWordOptionsNew(programs[i],&opt)==0);printf("options %d %d %d",programs[i],opt->window_size,opt->scan_range);value(opt->gap_trigger);value(opt->x_dropoff);printf(" %d\n",opt->program_number);BlastInitialWordOptionsFree(opt);}
    int checks[][2]={{-1,1},{0,1},{-1,0},{-32768,32767},{-32769,1},{-1,32768},{INT32_MIN,INT32_MAX},{-3,2}};
    for(int i=0;i<8;i++) {
        int lo=checks[i][0],hi=checks[i][1];Blast_ScoreFreq* s=Blast_ScoreFreqNew(lo,hi);
        printf("new %d %d %d %d",lo,hi,BlastScoreChk(lo,hi),s!=NULL);
        if(s){printf(" %d %d %d %d",s->score_min,s->score_max,s->obs_min,s->obs_max);value(s->score_avg);for(int j=lo;j<=hi;j++)value(s->sprob[j]);Blast_ScoreFreqFree(s);}putchar('\n');
    }
    Blast_KarlinBlk *zero=Blast_KarlinBlkNew();printf("zero");block(zero);Blast_KarlinBlkFree(zero);
    Blast_KarlinBlk from={1,2,3,4,5},to={9,9,9,9,9};
    for(int i=0;i<4;i++){int status=Blast_KarlinBlkCopy(i&1?&to:NULL,i&2?&from:NULL);printf("copy %d %d",i,status);block(&to);}
    int scores[][2]={{-2,1},{-3,2},{-1,1},{-1,3},{-3,1},{-6,4},{-4,2},{-5,3},{-3,4},{-8,5},{-10,2}};
    for(int i=0;i<11;i++)for(int mode=0;mode<8;mode++)frequency(scores[i][0],scores[i][1],mode);
    return 0;
}
