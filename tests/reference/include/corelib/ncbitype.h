/* Test-only substitute for generated toolkit configuration.
 * Only fixed-width typedefs; no algorithm implementation or BLAST macros.
 * The real BLAST core headers and C sources are used unchanged. */
#ifndef BTYPER_REFERENCE_NCBITYPE_H
#define BTYPER_REFERENCE_NCBITYPE_H
#include <stdint.h>
typedef int8_t Int1;
typedef uint8_t Uint1;
typedef int16_t Int2;
typedef uint16_t Uint2;
typedef int32_t Int4;
typedef uint32_t Uint4;
typedef int64_t Int8;
typedef uint64_t Uint8;
#ifdef __cplusplus
static_assert(sizeof(Int4) == 4 && sizeof(Int8) == 8, "NCBI integer widths");
#else
_Static_assert(sizeof(Int4) == 4 && sizeof(Int8) == 8, "NCBI integer widths");
#endif
#endif
