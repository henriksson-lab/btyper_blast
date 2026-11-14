// Original CSymDustMasker and nested algorithms, with a byte-vector input cursor.
#include "reference_symdust.hpp"
#include <iomanip>

static void bytes(const string& sequence)
{
    if (sequence.empty()) cout << '-';
    for (unsigned char value : sequence) cout << hex << setw(2) << setfill('0') << unsigned(value);
    cout << dec;
}

static void result(const CSymDustMasker::TMaskList& masks)
{
    if (masks.empty()) cout << '-';
    bool first = true;
    for (auto interval : masks) {
        if (!first) cout << ',';
        cout << interval.first << ':' << interval.second;
        first = false;
    }
}

static void record(const string& sequence, uint32_t level = 20, uint32_t window = 64, uint32_t linker = 1)
{
    CSymDustMasker masker(level, window, linker);
    objects::CSeqVector input(sequence);
    cout << "mask " << level << ' ' << window << ' ' << linker << ' ';
    bytes(sequence); cout << ' '; result(*masker(input)); cout << '\n';
}

int main()
{
    CRandom random;
    for (int i = 0; i < 1000; ++i) cout << "default_rng " << i << ' ' << random.GetRand() << '\n';
    random.Reset();
    for (int i = 0; i < 1000; ++i) cout << "reset_rng " << i << ' ' << random.GetRand() << '\n';
    random.SetSeed(4123);
    random.Reset();
    cout << "after_seed_reset " << random.GetRand() << '\n';

    // Every five-letter (ACGTN) word up to length six.
    for (int length = 0; length <= 6; ++length) {
        uint32_t combinations = 1;
        for (int i = 0; i < length; ++i) combinations *= 5;
        for (uint32_t word = 0; word < combinations; ++word) {
            string sequence(length, 'A');
            uint32_t value = word;
            for (int i = 0; i < length; ++i) { sequence[i] = "ACGTN"[value % 5]; value /= 5; }
            record(sequence);
        }
    }
    for (const string motif : {"A", "C", "G", "T", "N", "AT", "ACG", "ACGT", "RYSWKMBDHVN"}) {
        for (int length : {7, 8, 9, 15, 16, 31, 32, 62, 63, 64, 65, 66, 127, 128, 129, 512, 4096}) {
            string sequence;
            for (int i = 0; i < length; ++i) sequence += motif[i % motif.size()];
            record(sequence);
            // Enter and leave single-triplet windows and drain perfect intervals.
            record("ACGTGCACTG" + sequence + "TGCACTGAGC");
        }
    }
    uint32_t state = 0x15d42bc9;
    for (int variant = 0; variant < 200; ++variant) {
        string sequence;
        for (int i = 0; i < 64 + variant * 3; ++i) {
            state = state * 1664525u + 1013904223u;
            sequence += "ACGTNRYSWKMBVDH"[(state >> 16) % 15];
        }
        record(sequence);
        sequence.insert(sequence.size() / 2, string(64, "ACGTN"[variant % 5]));
        record(sequence);
    }
    const string special = "NNNACGT" + string(70, 'A') + "NNNN" + string(70, 'C') + "RYSWKMBDHVNNNN";
    for (uint32_t level : {0u, 1u, 2u, 20u, 64u, 65u})
        for (uint32_t window : {0u, 7u, 8u, 32u, 64u, 65u})
            for (uint32_t linker : {0u, 1u, 4u, 32u, 33u}) record(special, level, window, linker);

    CSymDustMasker reused;
    for (int iteration = 0; iteration < 20; ++iteration) {
        string sequence = string(100, 'N') + string(iteration * 3, 'A') + "ACGTNNNNNN";
        cout << "reuse " << iteration << ' '; bytes(sequence); cout << ' ';
        result(*reused(objects::CSeqVector(sequence))); cout << '\n';
    }
    for (uint32_t start : {0u, 1u, 3u, 63u, 64u, 128u, 1000u})
        for (uint32_t stop : {0u, 3u, 64u, 128u, 1000u}) {
            CSymDustMasker masker;
            cout << "range " << start << ' ' << stop << ' '; bytes(special); cout << ' ';
            result(*masker(objects::CSeqVector(special), start, stop)); cout << '\n';
        }
    return 0;
}
