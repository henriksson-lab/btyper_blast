// Test-only sequence-vector cursor and toolkit type declarations. The masker,
// all its nested algorithm classes, and method bodies are original source.
#include <algorithm>
#include <deque>
#include <iterator>
#include <list>
#include <memory>
#include <utility>
using Uint1 = uint8_t;
#define NCBI_XALGODUSTMASK_EXPORT
template<class T> class CConstRef;
template<class T> class CRef;

namespace objects {
class CSeq_id;
class CSeq_loc;
class CPacked_seqint;
class CSeqVector {
public:
    using size_type = uint32_t; // Original CSeqVectorTypes::size_type = TSeqPos.
    explicit CSeqVector(std::string input) : data_(std::move(input)) {}
    bool empty() const { return data_.empty(); }
    size_type size() const { return data_.size(); }
    class const_iterator {
    public:
        const_iterator(const CSeqVector& sequence, size_type position)
            : sequence_(sequence), position_(position) {}
        unsigned char operator*() const { return sequence_.data_.at(position_); }
        const_iterator& operator++() { ++position_; return *this; }
        size_type GetPos() const { return position_; }
        void SetPos(size_type position) { position_ = position; }
    private:
        const CSeqVector& sequence_;
        size_type position_;
    };
private:
    std::string data_;
};
}
