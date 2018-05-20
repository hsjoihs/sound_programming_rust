pub mod fft;
pub mod filter;
use std::f64::consts::PI;
pub mod wave;

#[allow(non_snake_case)]
pub fn safe_ADSR(
    e: &mut [f64],
    A: usize,
    D: usize,
    S: f64,
    R: usize,
    gate: usize,
    duration: usize,
) {
    if A != 0 {
        for n in 0..A {
            e[n] = 1.0 - (-5.0 * n as f64 / A as f64).exp();
        }
    }

    if D != 0 {
        for n in A..gate {
            e[n] = S + (1.0 - S) * (-5.0 * (n - A) as f64 / D as f64).exp();
        }
    } else {
        for n in A..gate {
            e[n] = S;
        }
    }

    if R != 0 {
        for n in gate..duration {
            e[n] = e[gate - 1] * (-5.0 * (n - gate + 1) as f64 / R as f64).exp();
        }
    }
}

#[derive(Clone)]
pub struct MonoPcm {
    pub fs: usize,
    pub bits: i32,
    pub length: usize,
    pub s: Vec<f64>,
}

impl MonoPcm {
    pub fn new16(fs: usize, length: usize) -> Self {
        MonoPcm {
            fs,
            length,
            bits: 16,
            s: vec![0.0; length],
        }
    }
    pub fn blank_copy(orig: &Self) -> Self {
        MonoPcm {
            s: vec![0.0; orig.length as usize],
            ..*orig
        }
    }

    pub fn mult_constant_gain(&mut self, gain: f64) {
        for n in 0..self.length {
            self.s[n] *= gain;
        }
    }

    pub fn new16_fn(fs: usize, length: usize, mut fun: Box<FnMut(usize) -> f64>) -> Self {
        MonoPcm {
            fs,
            length,
            bits: 16,
            s: (0..length).map(|n| fun(n)).collect(),
        }
    }
}

#[derive(Clone)]
pub struct StereoPcm {
    pub fs: usize,
    pub bits: i32,
    pub length: usize,
    pub s_l: Vec<f64>,
    pub s_r: Vec<f64>,
}

#[allow(non_snake_case)]
pub fn create_Hanning_window(N: usize) -> Vec<f64> {
    (0..N)
        .map(|n| {
            0.5 - 0.5 * (2.0 * PI * (n as f64 + if N % 2 == 0 { 0.0 } else { 0.5 }) / (N as f64))
                .cos()
        })
        .collect()
}

pub fn sinc(x: f64) -> f64 {
    if x == 0.0 {
        1.0
    } else {
        x.sin() / x
    }
}
