extern crate num_complex;
use self::num_complex::Complex;
use libc::c_double;
use std::f64::consts::PI;
//use libc::c_int;

fn log2_(x: usize) -> usize {
    /* y = log2(x) */
    let mut y = 0;
    let mut x = x;

    while x > 1 {
        x >>= 1;
        y += 1;
    }

    return y;
}

fn pow2(x: usize) -> usize /* y = 2 ^ x */
{
    let y: usize;

    if x == 0 {
        y = 1;
    } else {
        y = 2 << (x - 1);
    }

    return y;
}

#[allow(non_snake_case)]
pub fn safe_FFT_(x: &mut [Complex<f64>]) {
    let N = x.len();
    let number_of_stage = log2_(N); /* FFTの段数 */

    /* バタフライ計算 */
    for stage in 1..=number_of_stage {
        for i in 0..pow2(stage - 1) {
            for j in 0..pow2(number_of_stage - stage) {
                let n = pow2(number_of_stage - stage + 1) * i + j;
                let m = pow2(number_of_stage - stage) + n;
                let r = pow2(stage - 1) * j;
                let a = x[n];
                let b = x[m];
                let cc = Complex::new(((2.0 * PI * r as f64) / N as f64).cos(), -((2.0 * PI * r as f64) / N as f64).sin());
                if stage < number_of_stage {
                    x[n] = a + b;
                    x[m].re = (a.re - b.re) * cc.re - (a.im - b.im) * cc.im;
                    x[m].im = (a.im - b.im) * cc.re + (a.re - b.re) * cc.im;
                } else {
                    x[n] = a + b;
                    x[m] = a - b;
                }
            }
        }
    }

    /* インデックスの並び替えのためのテーブルの作成 */
    let mut index: Vec<usize> = vec![0; N];
    for stage in 1..=number_of_stage {
        for i in 0..pow2(stage - 1) {
            index[pow2(stage - 1) + i] = index[i] + pow2(number_of_stage - stage);
        }
    }

    /* インデックスの並び替え */
    for k in 0..N {
        if index[k] > k {
            let tmp = x[index[k]];
            x[index[k]] = x[k];
            x[k] = tmp;
        }
    }
}

#[allow(unused_variables, non_snake_case)]
pub fn safe_IFFT(x_real: &mut [c_double], x_imag: &mut [c_double]) {
    let N = x_real.len();
    assert_eq!(N, x_imag.len());

    let number_of_stage = log2_(N); /* IFFTの段数 */

    /* バタフライ計算 */
    for stage in 1..=number_of_stage {
        for i in 0..pow2(stage - 1) {
            for j in 0..pow2(number_of_stage - stage) {
                let n = pow2(number_of_stage - stage + 1) * i + j;
                let m = pow2(number_of_stage - stage) + n;
                let r = pow2(stage - 1) * j;
                let a_real = x_real[n];
                let a_imag = x_imag[n];
                let b_real = x_real[m];
                let b_imag = x_imag[m];
                let cc_real = ((2.0 * PI * r as f64) / N as f64).cos();
                let cc_imag = ((2.0 * PI * r as f64) / N as f64).sin();
                if stage < number_of_stage {
                    x_real[n] = a_real + b_real;
                    x_imag[n] = a_imag + b_imag;
                    x_real[m] = (a_real - b_real) * cc_real - (a_imag - b_imag) * cc_imag;
                    x_imag[m] = (a_imag - b_imag) * cc_real + (a_real - b_real) * cc_imag;
                } else {
                    x_real[n] = a_real + b_real;
                    x_imag[n] = a_imag + b_imag;
                    x_real[m] = a_real - b_real;
                    x_imag[m] = a_imag - b_imag;
                }
            }
        }
    }

    /* インデックスの並び替えのためのテーブルの作成 */
    let mut index: Vec<usize> = vec![0; N];
    for stage in 1..=number_of_stage {
        for i in 0..pow2(stage - 1) {
            index[pow2(stage - 1) + i] = index[i] + pow2(number_of_stage - stage);
        }
    }

    /* インデックスの並び替え */
    for k in 0..N {
        if index[k] > k {
            let real = x_real[index[k]];
            let imag = x_imag[index[k]];
            x_real[index[k]] = x_real[k];
            x_imag[index[k]] = x_imag[k];
            x_real[k] = real;
            x_imag[k] = imag;
        }
    }

    /* 計算結果をNで割る */
    for k in 0..N {
        x_real[k] /= N as f64;
        x_imag[k] /= N as f64;
    }
}
