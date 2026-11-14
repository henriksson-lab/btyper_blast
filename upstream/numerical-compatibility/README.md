# Known extended-exponential compatibility limit

The pinned BLAST+ **2.17.0+** `boost_erf.c` evaluates four expressions of the form `result *= expl(-z * z) / z;`. Its coefficients, polynomial intermediates, and the argument `-z*z` are double precision; `expl`, division, and multiplication then use the platform's `long double`. Our portable high/low adapter reproduces 64-bit-significand intermediate rounding with exponent scaling. That choice follows the measured GNU x86-64 reference. The C source itself does not prescribe the same long-double format or libm implementation on every platform.

One known isolated `ErfC` input, **6.78671432911633232**, differs from the original on this Ubuntu/glibc 2.35 host by one final double ULP:

| Stage | Original GNU C | Native adapter |
| --- | --- | --- |
| `expl` | `0xb.b6dd156e5e32552p-70` | `0xb.b6dd156e5e32551p-70` |
| Final double | `0x1.ed5acef190cf0p-71` | `0x1.ed5acef190cefp-71` |

Independent 160-digit Decimal arithmetic and 120-digit mpmath arithmetic agree that the native exponential, division, and multiplication each match correctly rounded 64-bit-significand stages **at this input**. The dynamically linked original `expl@GLIBC_2.2.5` is about **0.6833093254 extended ULP above the exact exponential**, one extended ULP above its correctly rounded value. Subsequent operations propagate that difference across a final-double rounding boundary. This is an observed libm compatibility difference, rather than a HighLow rounding bug at this input. It does not prove correct rounding of the adapter for all inputs, nor mathematical exactness of the enclosing Boost approximation.

The earlier 17,650-input check, with one nonidentical output and maximum one double ULP, remains accurate. Its former unexplained “numerical residual” is now explained by this audit. The 8,640 sampled Spouge outputs and 360 cutoffs were exact, and separately recorded whole-search raw comparisons remain valid for their tested inputs. None establishes universal bitwise numerical parity. Exact reproduction of a particular glibc approximation is a further compatibility requirement; no input-specific patch or production arithmetic change was made.

`erfc-glibc-audit.json` records the cause and scope; `erfc-stages.txt` preserves the measured stages. To reproduce on a GNU x86-64 host with 64-bit-significand long double, fetched pinned sources, `cc`, `rustc`, and Python:

```sh
python upstream/numerical-compatibility/reproduce_erfc.py
```

The script instruments temporary copies of the original polynomial and current Rust helper, prints both paths, and independently checks all three native stages using Python's standard-library Decimal. It changes no production files and adds no search fixtures. Other libm versions may produce a different original result; C exits with status 77 when `LDBL_MANT_DIG` is not 64.
