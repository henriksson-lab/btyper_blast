// Test-only declarations and exception adapters. All selected algorithm bodies
// are appended verbatim from the pinned source by the reference generator.
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using Uint4 = uint32_t;
#define _ASSERT assert
#define NCBI_THROW(type, code, message) throw std::runtime_error(message)

struct ReferenceRandomSupplier {
    bool IsInitialized() const { return false; }
};
static ReferenceRandomSupplier reference_random_supplier;
static ReferenceRandomSupplier* s_RandomSupplier = &reference_random_supplier;

class CRandom {
public:
    using TValue = Uint4;
    enum EGetRandMethod { eGetRand_LFG, eGetRand_Sys };
    explicit CRandom(EGetRandMethod method = eGetRand_LFG);
    explicit CRandom(TValue seed);
    void Reset(void);
    void SetSeed(TValue seed);
    TValue GetRand(void);
    TValue x_GetRand32Bits(void);
private:
    enum { kStateSize = 33 };
    EGetRandMethod m_RandMethod;
    TValue m_State[kStateSize];
    int m_RJ;
    int m_RK;
    TValue m_Seed;
    TValue x_GetSysRand32Bits(void) const {
        throw std::logic_error("System RNG is outside the selected reference scope");
    }
};
