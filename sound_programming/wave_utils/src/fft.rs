extern crate num_complex;
use self::num_complex::Complex;
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
                let c = Complex::new(0.0, -(2.0 * PI * r as f64) / N as f64).exp();
                if stage < number_of_stage {
                    x[n] = a + b;
                    x[m] = (a - b) * c;
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
pub fn safe_IFFT_(x: &mut [Complex<f64>]) {
    let N = x.len();

    let number_of_stage = log2_(N); /* IFFTの段数 */

    /* バタフライ計算 */
    for stage in 1..=number_of_stage {
        for i in 0..pow2(stage - 1) {
            for j in 0..pow2(number_of_stage - stage) {
                let n = pow2(number_of_stage - stage + 1) * i + j;
                let m = pow2(number_of_stage - stage) + n;
                let r = pow2(stage - 1) * j;
                let a = x[n];
                let b = x[m];
                let c = Complex::new(0.0, (2.0 * PI * r as f64) / N as f64).exp();

                if stage < number_of_stage {
                    x[n] = a + b;
                    x[m] = (a - b) * c;
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

    /* 計算結果をNで割る */
    for k in 0..N {
        x[k].re /= N as f64;
        x[k].im /= N as f64;
    }
}
