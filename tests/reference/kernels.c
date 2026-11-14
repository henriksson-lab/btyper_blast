/* Test driver for unmodified NCBI BLAST+ 2.17.0+ kernels.
 * Including blast_util.c permits testing its file-local s_CodonToAA as well.
 * Unused toolkit entry points are discarded by the linker, not reimplemented.
 */
#include <inttypes.h>
#include <strings.h>
#include "blast_util.c"
#include <algo/blast/core/gapinfo.h>
#include <algo/blast/core/ncbi_math.h>

static void bytes(const void *data, size_t length)
{
    const unsigned char *p = data;
    if (!length) putchar('-');
    for (size_t i = 0; i < length; ++i) printf("%02x", p[i]);
}

static uint64_t bits(double value)
{
    uint64_t result;
    memcpy(&result, &value, sizeof(result));
    return result;
}

#define TABLE(name, length) do { \
    printf("table %s ", #name); bytes(name, length); putchar('\n'); \
} while (0)

int main(void)
{
    /* gc.prt code 1, in T,C,A,G codon order, encoded with upstream table. */
    const char *code1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    Uint1 codes[64];
    for (int i = 0; i < 64; ++i) codes[i] = AMINOACID_TO_NCBISTDAA[(int)code1[i]];
    TABLE(NCBI4NA_TO_BLASTNA, 16);
    TABLE(BLASTNA_TO_NCBI4NA, 16);
    TABLE(BLASTNA_TO_IUPACNA, 16);
    TABLE(NCBI4NA_TO_IUPACNA, 16);
    TABLE(IUPACNA_TO_BLASTNA, 128);
    TABLE(IUPACNA_TO_NCBI4NA, 128);
    TABLE(AMINOACID_TO_NCBISTDAA, 128);
    TABLE(NCBISTDAA_TO_AMINOACID, 28);
    TABLE(codes, 64);
    printf("sentinels %u %u\n", kProtSentinel, kNuclSentinel);
    Uint1 *translation_tables[2];
    for (int reverse = 0; reverse <= 1; ++reverse) {
        translation_tables[reverse] = s_BlastGetTranslationTable(codes, reverse);
        printf("packed_table %d ", reverse); bytes(translation_tables[reverse], 64); putchar('\n');
    }
    printf("packed_table_null %d\n", s_BlastGetTranslationTable(NULL, FALSE) == NULL);

    EBlastProgramType programs[] = {
        eBlastTypeUndefined, eBlastTypeBlastp, eBlastTypeBlastn, eBlastTypeBlastx,
        eBlastTypeTblastn, eBlastTypeTblastx, eBlastTypePsiBlast, eBlastTypePsiTblastn,
        eBlastTypeRpsBlast, eBlastTypeRpsTblastn, eBlastTypePhiBlastp,
        eBlastTypePhiBlastn, eBlastTypeMapping
    };
    for (unsigned i = 0; i < sizeof(programs)/sizeof(programs[0]); ++i) {
        EBlastProgramType p = programs[i];
        unsigned flags = Blast_QueryIsProtein(p) | (Blast_QueryIsNucleotide(p) << 1)
            | (Blast_SubjectIsProtein(p) << 2) | (Blast_SubjectIsNucleotide(p) << 3)
            | (Blast_QueryIsTranslated(p) << 4) | (Blast_SubjectIsTranslated(p) << 5)
            | (Blast_ProgramIsNucleotide(p) << 6) | (Blast_ProgramIsValid(p) << 7);
        printf("program %u %u %u", p, flags, BLAST_GetNumberOfContexts(p));
        for (unsigned j = 0; j < 24; ++j) printf(" %d", BLAST_ContextToFrame(p, j));
        putchar('\n');
        if (Blast_QueryIsTranslated(p) || Blast_SubjectIsTranslated(p)) {
            for (Int2 frame = -3; frame <= 3; ++frame)
                if (frame) printf("frame %u %d %d\n", p, frame, BLAST_FrameToContext(frame, p));
        } else if (Blast_QueryIsNucleotide(p) || Blast_SubjectIsNucleotide(p)) {
            for (Int2 frame = -1; frame <= 1; frame += 2)
                printf("frame %u %d %d\n", p, frame, BLAST_FrameToContext(frame, p));
        } else {
            printf("frame %u 0 %d\n", p, BLAST_FrameToContext(0, p));
        }
    }
    Uint1 alphabet[17];
    for (int i = 0; i < 16; ++i) alphabet[i] = i;
    alphabet[16] = FENCE_SENTRY;
    for (int a = 0; a < 17; ++a)
        for (int b = 0; b < 17; ++b)
            for (int c = 0; c < 17; ++c) {
                Uint1 codon[3] = {alphabet[a], alphabet[b], alphabet[c]};
                printf("codon %u %u %u %u\n", codon[0], codon[1], codon[2], s_CodonToAA(codon, codes));
            }
    Uint1 sequence[4096], protein[4096];
    uint32_t state = 0x12345678;
    for (unsigned variant = 0; variant < 3; ++variant) {
        for (int i = 0; i < 4096; ++i) {
            state = state * UINT32_C(1664525) + UINT32_C(1013904223);
            sequence[i] = variant == 0 ? (1 << (state >> 30)) : (state >> 28);
            if (variant == 2 && i % 7 == 0) sequence[i] = FENCE_SENTRY;
        }
        for (int n = 0; n <= 129; ++n) {
            int length = n < 127 ? n : (n == 127 ? 511 : (n == 128 ? 512 : 4096));
            if (variant != 2) {
                for (int encoding = 1; encoding <= 2; ++encoding) {
                    Uint1 *packed = NULL;
                    if (BLAST_PackDNA(sequence, length, encoding, &packed)) return 1;
                    printf("pack %d ", encoding); bytes(sequence, length); putchar(' ');
                    bytes(packed, length / 4 + 1); putchar('\n');
                    free(packed);
                }
            }
            Uint1 *reverse = NULL;
            if (GetReverseNuclSequence(sequence, length, &reverse)) return 1;
            printf("reverse "); bytes(sequence, length); putchar(' ');
            bytes(reverse, length + 2); putchar('\n');
            for (Int2 frame = -3; frame <= 3; ++frame) {
                if (!frame) continue;
                Int4 prot_length = BLAST_GetTranslation(sequence, reverse, length, frame, protein, codes);
                printf("translation %d ", frame); bytes(sequence, length);
                printf(" %d ", prot_length); bytes(protein, prot_length + 2); putchar('\n');
            }
            free(reverse);
            if (length <= 5 || length == 11 || length == 64 || length == 126 || length >= 511) {
                for (int use_mixed = 0; use_mixed <= 1; ++use_mixed) {
                    if (use_mixed && length < 2) continue;
                    for (Int2 frame = -3; frame <= 3; ++frame) {
                        if (!frame) continue;
                        Uint1 *partial = NULL, *mixed = NULL;
                        Int4 prot_length = -1;
                        if (Blast_GetPartialTranslation(sequence, length, frame, codes, &partial,
                            &prot_length, use_mixed ? &mixed : NULL)) return 1;
                        printf("partial %d %d ", frame, use_mixed); bytes(sequence, length);
                        printf(" %d ", prot_length);
                        bytes(partial, use_mixed ? length + 2 : prot_length + 2);
                        putchar(' '); bytes(mixed, use_mixed ? length + 1 : 0); putchar('\n');
                        free(partial); free(mixed);
                    }
                }
            }
            if (length >= 4) {
                for (int encoding = eBlastEncodingNcbi4na; encoding <= eBlastEncodingNcbi2na; ++encoding) {
                    Uint1 *packed = NULL;
                    const Uint1 *input = sequence;
                    int input_length = length;
                    if (encoding == eBlastEncodingNcbi2na) {
                        if (variant == 2) continue; /* Fences are not valid NCBI4na table indexes. */
                        if (BLAST_PackDNA(sequence, length, eBlastEncodingNcbi4na, &packed)) return 1;
                        input = packed;
                        input_length = length / 4 + 1;
                        for (Int2 frame = -3; frame <= 3; ++frame) {
                            if (!frame) continue;
                            memset(protein, 0xee, sizeof(protein));
                            Int4 prot_length = BLAST_TranslateCompressedSequence(
                                translation_tables[frame < 0], length, input, frame, protein);
                            printf("compressed %d %d ", frame, length); bytes(input, input_length);
                            printf(" %d ", prot_length); bytes(protein, prot_length + 2); putchar('\n');
                        }
                    }
                    Uint1 *all = NULL, *mixed = NULL;
                    Uint4 *offsets = NULL;
                    if (BLAST_GetAllTranslations(input, encoding, length, codes, &all, &offsets, &mixed)) return 1;
                    printf("alltranslations %d %d ", encoding, length); bytes(input, input_length);
                    for (int context = 0; context <= NUM_FRAMES; ++context) printf(" %u", offsets[context]);
                    putchar(' '); bytes(all, offsets[NUM_FRAMES]);
                    putchar(' '); bytes(mixed, 2 * length + 3); putchar('\n');
                    free(all); free(offsets); free(mixed); free(packed);
                }
            }
        }
    }
    for (int length = 0; length <= 129; ++length) {
        printf("protein_length %d", length);
        for (unsigned context = 0; context < 18; ++context)
            printf(" %zu", BLAST_GetTranslatedProteinLength(length, context));
        putchar('\n');
    }
    for (int length = 0; length <= 2; ++length)
        for (Int2 frame = -3; frame <= 3; ++frame) {
            if (!frame) continue;
            memset(protein, 0xee, sizeof(protein));
            Int4 prot_length = BLAST_TranslateCompressedSequence(translation_tables[frame < 0], length, sequence, frame, protein);
            printf("compressed %d %d ", frame, length); bytes(sequence, length / 4 + 1);
            printf(" %d ", prot_length); bytes(protein, prot_length + 2); putchar('\n');
        }
    memset(protein, 0xee, sizeof(protein));
    printf("compressed_null %d %d %u\n",
        BLAST_TranslateCompressedSequence(translation_tables[0], 12, NULL, 1, protein),
        BLAST_TranslateCompressedSequence(translation_tables[0], 12, sequence, 1, NULL), protein[0]);
    free(translation_tables[0]); free(translation_tables[1]);
    for (int size = -1; size <= 8; ++size) {
        GapEditScript *script = GapEditScriptNew(size);
        printf("edit %d %d", size, script != NULL);
        if (script) {
            printf(" %d", script->size);
            for (int i = 0; i < script->size; ++i)
                printf(" %d:%d", script->op_type[i], script->num[i]);
        }
        printf(" %d\n", GapEditScriptDelete(script) == NULL);
    }
    GapPrelimEditBlock *block = GapPrelimEditBlockNew();
    for (int i = 0; i < 1300; ++i) {
        /* Adjacent equal ops, zero counts, several allocation thresholds and reuse. */
        int op = (i / 3) % 3 == 0 ? eGapAlignDel : ((i / 3) % 3 == 1 ? eGapAlignSub : eGapAlignIns);
        int count = i % 5 == 0 ? 0 : i % 11;
        if (i == 900) GapPrelimEditBlockReset(block);
        else GapPrelimEditBlockAdd(block, op, count);
        /* Snapshot selected steps, including boundaries; inputs are recorded. */
        if (i < 12 || i % 50 == 0 || (i >= 890 && i <= 910) || i == 1299) {
            printf("prelim %d %d %d %d", i, block->num_ops_allocated, block->num_ops, block->last_op);
            for (int j = 0; j < block->num_ops; ++j)
                printf(" %d:%d", block->edit_ops[j].op_type, block->edit_ops[j].num);
            putchar('\n');
        }
    }
    printf("free %d\n", GapPrelimEditBlockFree(block) == NULL);

    const double inputs[] = {
        -100.0, -10.0, -1.0, -0.999999, -0.5, -0.33000000000000007, -0.33,
        -0.2, -0.19999999999999998, -0.01, -1.e-15, -1.e-16, -1.e-17,
        -0.0, 0.0, 1.e-17, 1.e-16, 1.e-15, 0.01, 0.19999999999999998, 0.2,
        0.33, 0.33000000000000007, 0.5, 1.0, 10.0, 100.0
    };
    for (unsigned i = 0; i < sizeof(inputs)/sizeof(inputs[0]); ++i) {
        double x = inputs[i];
        printf("expm1 %016" PRIx64 " %016" PRIx64 "\n", bits(x), bits(BLAST_Expm1(x)));
        printf("log1p %016" PRIx64 " %016" PRIx64 "\n", bits(x), bits(BLAST_Log1p(x)));
        printf("nint %016" PRIx64 " %ld\n", bits(x), BLAST_Nint(x));
        for (int n = -10; n <= 10; ++n)
            printf("powi %016" PRIx64 " %d %016" PRIx64 "\n", bits(x), n, bits(BLAST_Powi(x, n)));
    }
    for (int a = -15; a <= 15; ++a)
        for (int b = -15; b <= 15; ++b) {
            printf("gcd %d %d %d\n", a, b, BLAST_Gcd(a, b));
            for (int c = -3; c <= 3; ++c) {
                int aa = a, bb = b, cc = c;
                int g = BLAST_Gdb3(&aa, &bb, &cc);
                printf("gdb3 %d %d %d %d %d %d %d\n", a, b, c, g, aa, bb, cc);
            }
        }
    return 0;
}
