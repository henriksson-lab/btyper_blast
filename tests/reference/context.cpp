// Unchanged setup functions plus adapters for query data and object ownership.
#include "reference_context.hpp"
#include <iomanip>
#include <climits>

static void bytes(const Uint1* data,size_t size)
{
    if (!size) cout << '-';
    for (size_t i=0;i<size;++i) cout << hex << setw(2) << setfill('0') << unsigned(data[i]);
    cout << dec;
}
static void intervals(const BlastSeqLoc* loc)
{
    int count=0; for (auto p=loc;p;p=p->next) ++count;
    cout << count;
    for (;loc;loc=loc->next) cout << ' ' << loc->ssr->left << ' ' << loc->ssr->right;
}
static void masks(const BlastMaskLoc* mask)
{
    for (int i=0;i<mask->total_size;++i) { cout << ' '; intervals(mask->seqloc_array[i]); }
}
static void contexts(const BlastQueryInfo* info)
{
    cout << ' ' << info->first_context << ' ' << info->last_context << ' ' << info->num_queries
         << ' ' << info->max_length << ' ' << info->min_length;
    for (int i=0;i<=info->last_context;++i) {
        const auto& c=info->contexts[i];
        cout << ' ' << c.query_offset << ' ' << c.query_length << ' ' << c.eff_searchsp
             << ' ' << c.length_adjustment << ' ' << c.query_index << ' ' << int(c.frame)
             << ' ' << int(c.is_valid) << ' ' << c.segment_flags;
    }
}
static void complement_case(BlastQueryInfo* info,EBlastProgramType program,int variant)
{
    BlastMaskLoc* mask=BlastMaskLocNew(info->last_context+1);
    for (int i=0;i<=info->last_context;++i) {
        int length=info->contexts[i].query_length;
        if (length<1 || variant==0) continue;
        BlastSeqLoc** head=&mask->seqloc_array[i];
        if (variant==1) BlastSeqLocNew(head,0,length-1);
        if (variant==2) BlastSeqLocNew(head,0,0);
        if (variant==3) BlastSeqLocNew(head,length-1,length-1);
        if (variant>=4) {
            // Ordered, disjoint masks; retain short leading/trailing gaps.
            for (int start=variant%3;start<length;start+=5) BlastSeqLocNew(head,start,MIN(start+1,length-1));
        }
    }
    bool present=variant!=6, output=variant!=7;
    cout << "complement " << int(program) << ' ' << present << ' ' << output;
    masks(mask);
    BlastSeqLoc* result=NULL;
    int status=BLAST_ComplementMaskLocations(program,info,present?mask:NULL,output?&result:NULL);
    cout << ' ' << status << ' ';intervals(result);masks(mask);cout << '\n';
    BlastSeqLocFree(result);BlastMaskLocFree(mask);
}
static IBlastSeqVector strands(const vector<Uint1>& plus)
{
    IBlastSeqVector sv;sv.plus=plus;
    Uint1* minus=NULL;ASSERT(GetReverseNuclSequence(plus.data(),plus.size(),&minus)==0);
    sv.minus.assign(minus+1,minus+plus.size()+1);free(minus);
    return sv;
}
static void strand_case(const vector<Uint1>& plus)
{
    auto sv=strands(plus);
    for (auto enc:{eBlastEncodingNucleotide,eBlastEncodingNcbi4na})
    for (auto sent:{eSentinels,eNoSentinels})
    for (auto strand:{eNa_strand_plus,eNa_strand_minus,eNa_strand_both}) {
        auto value=strand==eNa_strand_both?GetSequenceNucleotideBothStrands(sv,enc,sent):GetSequenceSingleNucleotideStrand(sv,enc,strand,sent);
        cout << "strand " << int(enc) << ' ' << int(sent) << ' ' << int(strand) << ' ';
        bytes(sv.plus.data(),sv.size());cout << ' ';bytes(sv.minus.data(),sv.size());cout << ' ';
        bytes(value.data.get(),value.length);cout << '\n';
    }
}
static void batch_case(const vector<string>& input)
{
    auto info=BlastQueryInfoNew(eBlastTypeBlastn,input.size());
    vector<IBlastSeqVector> data;
    auto mask=BlastMaskLocNew(info->last_context+1);
    info->min_length=INT_MAX;
    cout << "batch " << input.size();
    for (size_t i=0;i<input.size();++i) {
        cout << ' ';bytes(reinterpret_cast<const Uint1*>(input[i].data()),input[i].size());
        vector<Uint1> plus;
        for (unsigned char letter:input[i]) plus.push_back(IUPACNA_TO_NCBI4NA[letter]);
        data.push_back(strands(plus));
        s_QueryInfo_SetContext(info,2*i,plus.size());s_QueryInfo_SetContext(info,2*i+1,plus.size());
        info->max_length=MAX(info->max_length,plus.size());info->min_length=MIN(info->min_length,plus.size());
        CSymDustMasker masker;
        auto dust=masker(objects::CSeqVector(input[i]));
        for (const auto& range:*dust) {
            BlastSeqLocNew(&mask->seqloc_array[2*i],range.first,range.second);
            BlastSeqLocNew(&mask->seqloc_array[2*i+1],range.first,range.second);
        }
    }
    Uint4 buflen=QueryInfo_GetSeqBufLen(info);
    vector<Uint1> buffer(buflen+1);
    for (size_t i=0;i<input.size();++i) {
        auto sequence=GetSequenceNucleotideBothStrands(data[i],eBlastEncodingNucleotide,eSentinels);
        memcpy(buffer.data()+info->contexts[2*i].query_offset,sequence.data.get(),sequence.length);
    }
    contexts(info);cout << ' ';bytes(buffer.data(),buflen);
    BlastSeqLoc* included=NULL;
    ASSERT(BLAST_ComplementMaskLocations(eBlastTypeBlastn,info,mask,&included)==0);
    cout << ' ';intervals(included);cout << '\n';
    for (int word:{11,28}) {
        LookupTableOptions options={0};options.word_size=word;options.program_number=eBlastTypeBlastn;
        Int4 maxoff,width;
        int estimate=EstimateNumTableEntries(included,&maxoff);
        int kind=BlastChooseNaLookupTable(&options,estimate,maxoff,&width);
        cout << "batch_lookup " << word << ' ' << kind << ' ' << width << '\n';
    }
    BlastSeqLocFree(included);BlastMaskLocFree(mask);BlastQueryInfoFree(info);
}
int main()
{
    EBlastProgramType programs[]={eBlastTypeBlastn,eBlastTypeTblastn,eBlastTypeBlastx,eBlastTypeTblastx,eBlastTypeMapping};
    for (auto program:programs) {
        for (int n:{0,-1}) cout << "new_none " << int(program) << ' ' << n << ' ' << (BlastQueryInfoNew(program,n)==NULL) << '\n';
        for (int count:{1,2,3}) for (int variant=0;variant<6;++variant) {
            auto info=BlastQueryInfoNew(program,count);
            cout << "new " << int(program) << ' ' << count;contexts(info);cout << '\n';
            vector<Uint4> lengths;
            for (int i=0;i<=info->last_context;++i) {
                unsigned length=(i*7+variant*3)%19;
                if (variant==0) length=0;
                lengths.push_back(length);s_QueryInfo_SetContext(info,i,length);
            }
            int reassign=variant==2?info->last_context:-1;
            if (reassign>=0) { s_QueryInfo_SetContext(info,reassign,0);s_QueryInfo_SetContext(info,reassign,3); }
            info->max_length=0;info->min_length=INT_MAX;
            for (int i=0;i<=info->last_context;++i) {
                info->max_length=MAX(info->max_length,info->contexts[i].query_length);
                info->min_length=MIN(info->min_length,info->contexts[i].query_length);
            }
            if (variant==4) info->first_context=MIN(1,info->last_context);
            cout << "info " << int(program) << ' ' << count << ' ' << lengths.size();
            for (auto length:lengths) cout << ' ' << length;
            cout << ' ' << reassign << ' ' << info->first_context;contexts(info);
            cout << ' ' << QueryInfo_GetSeqBufLen(info);
            cout << '\n';
            for (int offset=0;offset<int(QueryInfo_GetSeqBufLen(info))+3;++offset)
                cout << "offset " << offset << ' ' << BSearchContextInfo(offset,info) << ' ' << Blast_GetQueryIndexFromQueryOffset(offset,program,info) << '\n';
            for (int v=0;v<8;++v) complement_case(info,program,v);
            BlastQueryInfoFree(info);
        }
    }
    for (Uint4 length:{0u,1u,2u,3u,31u,1024u,0x7fffffffu,0xffffffffu}) {
        for (auto encoding:{eBlastEncodingNucleotide,eBlastEncodingNcbi4na})
        for (auto strand:{eNa_strand_unknown,eNa_strand_plus,eNa_strand_minus,eNa_strand_both,eNa_strand_both_rev,eNa_strand_other})
        for (auto sentinel:{eSentinels,eNoSentinels})
            cout << "buffer " << length << ' ' << int(encoding) << ' ' << int(strand) << ' ' << int(sentinel) << ' ' << CalculateSeqBufferLength(length,encoding,strand,sentinel) << '\n';
        cout << "buffer " << length << " 0 0 0 " << CalculateSeqBufferLength(length,eBlastEncodingProtein,eNa_strand_unknown,eSentinels) << '\n';
        cout << "buffer " << length << " 3 1 1 " << CalculateSeqBufferLength(length,eBlastEncodingNcbi2na,eNa_strand_plus,eNoSentinels) << '\n';
    }
    for (auto encoding:{eBlastEncodingProtein,eBlastEncodingNucleotide,eBlastEncodingNcbi4na})
        cout << "sentinel " << int(encoding) << ' ' << int(GetSentinelByte(encoding)) << '\n';
    for (int length=1;length<=2;++length) {
        int combinations=1<<(4*length);
        for (int word=0;word<combinations;++word) {
            vector<Uint1> plus(length);int value=word;
            for (auto& base:plus) { base=value&15;value>>=4; }
            strand_case(plus);
        }
    }
    uint32_t state=3120971u;
    for (int length:{3,7,8,9,31,32,33,127,128,129,1024}) {
        vector<Uint1> plus(length);
        for (auto& base:plus) { state=state*1664525u+1013904223u;base=(state>>28); }
        strand_case(plus);
    }
    for (int length:{1,2,7,8,9,15,16,31,32,63,64,65,127,128,129,512,4096}) {
        string random,repeated,ambiguous;
        for (int i=0;i<length;++i) {
            state=state*1664525u+1013904223u;random+="ACGT"[(state>>28)%4];
            repeated+="ACG"[i%3];ambiguous+="ACGTNRYSWKMBVDH"[(state>>28)%15];
        }
        batch_case({random});batch_case({repeated});batch_case({ambiguous});
        batch_case({random,repeated,ambiguous});
        if (length==4096) batch_case({random,random});
    }
    batch_case({string(17000,'N'),"ACGTACGTACGTACGTACGTACGT"});
    batch_case({"ACGT"+string(80,'A')+"CGATC",string(90,'N'),"GGGTCGTGCACTGATGCTACGATCGTAGCATCGT"});
    return 0;
}
