// Exercises selected unchanged NCBI C++ methods compiled via reference_database.hpp.
#include "reference_database.hpp"
#include <iomanip>

static void bytes(const string& value)
{
    if (value.empty()) cout << '-';
    for (unsigned char byte : value) {
        cout << hex << setw(2) << setfill('0') << unsigned(byte);
    }
    cout << dec;
}

static void record(const vector<unsigned char>& sequence)
{
    string ncbi4na((sequence.size() + 1) / 2, '\0'), seq, amb;
    for (size_t i = 0; i < sequence.size(); ++i)
        ncbi4na[i / 2] |= sequence[i] << (i % 2 ? 0 : 4);
    WriteDB_Ncbi4naToBinary(ncbi4na.data(), ncbi4na.size(), sequence.size(), seq, amb);
    cout << "record ";
    bytes(string(sequence.begin(), sequence.end())); cout << ' ';
    bytes(seq); cout << ' '; bytes(amb); cout << '\n';
}

int main()
{
    for (uint32_t seed : {0u, 1u, 2u, 3u, 127u, 4096u, 0xffffffffu}) {
        CRandom generator(seed);
        for (int i = 0; i < 1000; ++i)
            cout << "rng " << seed << ' ' << i << ' ' << generator.GetRand() << '\n';
        generator.SetSeed(seed);
        cout << "reseed " << seed << ' ' << generator.GetRand() << '\n';
        for (int i = 0; i < 100; ++i)
            cout << "rng32 " << seed << ' ' << i << ' ' << generator.x_GetRand32Bits() << '\n';
    }
    record({});
    // Every combination of ambiguity nibbles, including the raw gap code zero.
    for (int length = 1; length <= 4; ++length) {
        uint32_t combinations = 1u << (length * 4);
        for (uint32_t word = 0; word < combinations; ++word) {
            vector<unsigned char> sequence(length);
            for (int i = 0; i < length; ++i) sequence[i] = (word >> (i * 4)) & 15;
            record(sequence);
        }
    }
    for (int length : {5, 7, 8, 11, 14, 15, 16, 17, 31, 32, 33, 127, 128, 129, 4094, 4095, 4096, 4097}) {
        for (int code : {3, 5, 7, 10, 11, 12, 13, 14, 15}) {
            record(vector<unsigned char>(length, code));
        }
        vector<unsigned char> sequence(length);
        uint32_t state = 0x3f075a19;
        for (int i = 0; i < length; ++i) {
            state = state * 1664525u + 1013904223u;
            sequence[i] = state >> 28;
        }
        record(sequence);
    }
    // Exercise the large-sequence format decision without allocating a huge sequence.
    for (int size : {0xffffff, 0x1000000, 0x1000001}) {
        CAmbigDataBuilder builder(size);
        for (int offset : {0, 17, 0xfffffe}) builder.Check(15, offset);
        string amb;
        builder.GetAmbig(amb);
        cout << "large_ambiguity " << size << ' '; bytes(amb); cout << '\n';
    }
    return 0;
}
