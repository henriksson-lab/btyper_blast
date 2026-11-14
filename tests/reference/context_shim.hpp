// Test-only input/ownership adapters for the original setup helper functions.
// Sequence conversion, stitching, context setup and masking bodies are original.
#include <algo/blast/core/blast_encoding.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_query_info.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_nalookup.h>
#include <algo/blast/core/lookup_util.h>
#include <cstring>
#include <memory>
#include <string>
#include <vector>
#include "reference_symdust.hpp"

#define THROWS(x)
using TSeqPos = Uint4;
#include "reference_na_strand.hpp"
namespace objects {
struct CSeq_data { enum { e_Ncbi4na=4 }; };
}
using namespace objects;
namespace NStr { static std::string IntToString(unsigned value) { return std::to_string(value); } }
struct FreeUint1 { void operator()(Uint1* p) const { free(p); } };
using TAutoUint1Ptr = std::unique_ptr<Uint1,FreeUint1>;
struct SBlastSequence {
    TAutoUint1Ptr data;
    TSeqPos length;
    SBlastSequence(Uint1* p,TSeqPos n):data(p),length(n) {}
};
class IBlastSeqVector {
public:
    std::vector<Uint1> plus,minus;
    TSeqPos size() const { return plus.size(); }
    void SetCoding(int coding) { assert(coding==CSeq_data::e_Ncbi4na); }
    Uint1 operator[](TSeqPos pos) const { return plus[pos]; }
    void GetStrandData(ENa_strand strand,Uint1* output) const {
        const auto& input = strand == eNa_strand_plus ? plus : minus;
        memcpy(output,input.data(),input.size());
    }
};
