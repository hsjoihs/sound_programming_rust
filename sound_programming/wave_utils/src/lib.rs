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

    pub fn mult_varying_gain(&mut self, vec: &[f64], gain: f64) {
        assert!(vec.len() >= self.length);
        for n in 0..self.length {
            self.s[n] *= vec[n] * gain
        }
    }

    pub fn new16_fn(fs: usize, length: usize, mut fun: Box<dyn FnMut(usize) -> f64>) -> Self {
        MonoPcm {
            fs,
            length,
            bits: 16,
            s: (0..length).map(|n| fun(n)).collect(),
        }
    }

    pub fn new16_sawtooth_with_varying_freq(fs: usize, length: usize, f0: &[f64]) -> Self {
        MonoPcm{
            fs,
            length,
            bits: 16,
            s: sawtooth_with_varying_freq(fs, length, f0)
        }
    }
}

pub fn linear(initial_v: f64, final_v: f64, length: usize) -> Vec<f64> {
    (0..length)
        .map(|n| initial_v + (final_v - initial_v) * n as f64 / (length - 1) as f64)
        .collect()
}

pub fn sawtooth_with_varying_freq(pcm0_fs: usize, pcm0_length: usize, f0: &[f64]) -> Vec<f64> {
    assert!(f0.len() >= pcm0_length);
    let mut pcm0_s = vec![0.0; pcm0_length];

    /* ノコギリ波 */
    let mut t0 = (pcm0_fs as f64 / f0[0]) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..pcm0_length {
        pcm0_s[n] = 1.0 - 2.0 * m as f64 / t0 as f64;

        m += 1;
        if m >= t0 {
            t0 = (pcm0_fs as f64 / f0[n]) as usize; /* 基本周期 */
            m = 0;
        }
    }
    pcm0_s
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

#[allow(non_snake_case)]
pub fn determine_J(delta: f64) -> usize {
    let mut J = (3.1 / delta + 0.5) as usize - 1; /* 遅延器の数 */
    if J % 2 == 1 {
        J += 1; /* J+1が奇数になるように調整する */
    }
    return J;
}

pub fn lfo(
    pcm: &MonoPcm,
    center: f64,
    am: f64, /* LFOの振幅 */
    fm: f64, /* LFOの周波数 */
) -> Vec<f64> {
    (0..pcm.length)
        .map(|n| center + am * (2.0 * PI * fm * n as f64 / pcm.fs as f64).sin())
        .collect()
}

pub fn mult(i: usize, d: f64) -> usize {
    ((i as f64) * d) as usize
}
