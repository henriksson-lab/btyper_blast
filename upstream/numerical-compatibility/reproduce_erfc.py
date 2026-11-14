#!/usr/bin/env python3
"""Reproduce one libm compatibility difference, without changing production files."""
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

REPO = Path(__file__).resolve().parents[2]
CORE = REPO / 'vendor/ncbi-blast-2.17.0+-src/c++/src/algo/blast/core'

C_HELPER = r'''
#include <float.h>
#include <stdio.h>
#include <math.h>
static double debug_product(double r, double z) {
    double arg = -z*z;
    long double e = expl(arg), d = e/z, m = (long double)r*d;
    printf("C_PRE %.17g\nC_ARG %.17g\n", r, arg);
    printf("C_EXP %La\nC_DIV %La\nC_MUL %La\n", e, d, m);
    return (double)m;
}
'''
R_HELPER = r'''
fn extended_exp_product(result:f64,z:f64)->f64 {
    let (value,power)=HighLow::exp_reduced(-z*z);
    let e=value.round_extended();
    let d=e.div_scalar(z).round_extended();
    let m=d.mul_scalar(result).round_extended();
    println!("R_EXP {:.17e} {:.17e} {}",e.hi,e.lo,power);
    println!("R_DIV {:.17e} {:.17e} {}",d.hi,d.lo,power);
    println!("R_MUL {:.17e} {:.17e} {}",m.hi,m.lo,power);
    m.round_scaled(power)
}
fn main(){println!("R_FINAL {:016x}",erf_c(6.78671432911633232).to_bits());}
'''


def rounded64(x):
    exponent = math.frexp(float(x))[1] - 1
    quantum = Decimal(2) ** (exponent - 63)
    return (x / quantum).to_integral_value(rounding=ROUND_HALF_EVEN) * quantum


def main():
    source = (CORE / 'boost_erf.c').read_text()
    expression = 'result *= expl(-z * z) / z;'
    assert source.count(expression) == 4
    source = C_HELPER + source.replace(expression, 'result = debug_product(result,z);')
    source += '\nint main(void){if(LDBL_MANT_DIG!=64)return 77;printf("C_FINAL %a\\n",ErfC(6.78671432911633232));}\n'
    rust = (REPO / 'src/core/boost_erf.rs').read_text()
    # Retain the production polynomial and HighLow code; instrument only its
    # final expression wrapper in a temporary standalone compilation unit.
    rust = rust[:rust.index('fn extended_exp_product(')] + R_HELPER
    with tempfile.TemporaryDirectory(prefix='btyper-erfc-audit-') as temporary:
        work = Path(temporary)
        (work / 'stage.c').write_text(source)
        (work / 'stage.rs').write_text(rust)
        subprocess.run(shlex.split(os.environ.get('CC', 'cc')) + [
            '-std=c11', '-D_DEFAULT_SOURCE', '-O2', '-ffp-contract=off',
            '-I' + str(REPO / 'tests/reference/include'),
            '-I' + str(CORE.parents[3] / 'include'), '-I' + str(CORE),
            str(work / 'stage.c'), '-lm', '-o', str(work / 'stage-c')], check=True)
        subprocess.run(['rustc', '--edition=2021', '-O', str(work / 'stage.rs'),
                        '-o', str(work / 'stage-r')], check=True)
        c = subprocess.check_output([str(work / 'stage-c')], text=True)
        rust_output = subprocess.check_output([str(work / 'stage-r')], text=True)
    print(c + rust_output, end='')
    rows = dict(line.split(maxsplit=1) for line in (c + rust_output).splitlines())
    with localcontext() as context:
        context.prec = 160
        arg = Decimal.from_float(float(rows['C_ARG']))
        z = Decimal.from_float(6.78671432911633232)
        coefficient = Decimal.from_float(float(rows['C_PRE']))
        exact_exp = arg.exp()
        stages = [rounded64(exact_exp)]
        stages.append(rounded64(stages[-1] / z))
        stages.append(rounded64(stages[-1] * coefficient))
        for name, expected in zip(('EXP', 'DIV', 'MUL'), stages):
            hi, lo, power = rows['R_' + name].split()
            actual = (Decimal.from_float(float(hi)) + Decimal.from_float(float(lo))) * Decimal(2)**int(power)
            assert actual == expected, (name, actual, expected)
            print(f'Decimal correctly-rounded {name}: {expected}')
        print('All three native stages equal independently rounded p64 stages.')
        print('Correctly-rounded-stage final f64:', float(stages[-1]).hex())
        # Decode the original hexadecimal long-double exp exactly.
        match = re.fullmatch(r'0x([0-9a-f]+)\.([0-9a-f]+)p([+-]?\d+)', rows['C_EXP'])
        integer, fraction, power = match.groups()
        original_exp = Decimal(int(integer + fraction, 16)) * Decimal(2)**(int(power) - 4*len(fraction))
        ulp = Decimal(2)**(math.frexp(float(exact_exp))[1]-64)
        print('Original expl error in p64 ULP:', (original_exp-exact_exp)/ulp)


if __name__ == '__main__':
    main()
