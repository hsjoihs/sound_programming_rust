use determine_J;
use linear;
use num_complex::Complex;

use mult;
extern crate rand;
use rand::Rng;
use sine_wave;
use std::f64::consts::PI;
use wave_utils::MonoPcm;
use wave_utils::create_Hanning_window;
use wave_utils::fft::safe_FFT_;
use wave_utils::fft::safe_IFFT_;
use wave_utils::filter::safe_FIR_LPF;
use wave_utils::filter::safe_FIR_filtering;
use wave_utils::filter::safe_IIR_LPF;
use wave_utils::filter::safe_IIR_filtering;
use wave_utils::wave::wave_read_16bit_mono_safer3;
use wave_utils::wave::wave_read_16bit_stereo_safer3;
use wave_utils::wave::wave_write_16bit_mono_safer3;
use wave_utils::wave::wave_write_16bit_stereo_safer3;

#[allow(non_snake_case)]
fn verify_(X: Vec<Complex<f64>>) {
    /* 周波数特性 */
    for (k, item) in X.iter().enumerate() {
        assert_close(item.re, 0.0);
        assert_close(
            item.im,
            match k {
                4 => -16.0,
                60 => 16.0,
                _ => 0.0,
            },
        );
    }
}

fn exponential_decay(
    pcm_fs: usize,
    amplitude: f64,
    decay_factor: f64,
    decay_time: f64,
    length: usize,
) -> Vec<f64> {
    (0..length)
        .map(|n| amplitude * (-decay_factor * n as f64 / (pcm_fs as f64 * decay_time)).exp())
        .collect()
}

#[allow(non_snake_case)]
fn dft(N: usize, func: Box<Fn(usize) -> f64>) -> Vec<Complex<f64>> {
    let mut x: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); N];
    let mut X: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); N];

    let pcm_slice = wave_read_16bit_mono_safer3("sine_500hz.wav").s;

    /* 波形 */
    for n in 0..N {
        x[n] = Complex::new(pcm_slice[n] * func(n) /* x(n)の実数部 */, 0.0); /* x(n)の虚数部 */
    }

    /* DFT */
    for k_ in 0..N {
        let k = k_ as f64;
        for n_ in 0..N {
            let n = n_ as f64;
            let N = N as f64;
            let W = Complex::new(0.0, -(2.0 * PI * k * n / N)).exp();
            X[k_] += W * x[n_];
        }
    }

    X
}

fn assert_close(a: f64, b: f64) {
    assert!((a - b).abs() < 1.75e-4);
}

pub fn first() {
    ex1_1();
    ex1_2();
    ex2_1();
    ex2_2();
    ex3_1();
    ex3_2();
    ex3_3();
    ex3_4();
    if false {
        ex3_5(); // slow and random
    }
    ex4_1();
    ex4_2();
    ex4_3();
    ex4_4();
    ex5_1();
    ex5_2();
    ex5_3();
    ex5_4();
    ex5_5();
    ex6_1();
    ex6_2();
    ex6_3();
    ex6_4();
    if false {
        ex6_5(); // slooooow
    }
}

fn ex1_1() {
    let pcm0 = wave_read_16bit_mono_safer3("ex1_1_a.wav"); /* 音データの入力 */
    let pcm1 = pcm0.clone(); /* 音データのコピー */

    wave_write_16bit_mono_safer3("ex1_1_b.wav", &pcm1); /* 音データの出力 */
}

#[allow(non_snake_case)]
fn ex1_2() {
    let pcm0 = wave_read_16bit_stereo_safer3("ex1_2_a.wav");
    let pcm1 = pcm0.clone(); /* 音データのコピー */

    wave_write_16bit_stereo_safer3("ex1_2_b.wav", &pcm1);
}

fn ex2_1() {
    let pcm_fs: usize = 44100; /* 標本化周波数 */
    let pcm_length: usize = pcm_fs * 1; /* 音データの長さ */

    /* サイン波 */
    let pcm = MonoPcm::new16_fn(
        pcm_fs,
        pcm_length,
        Box::new(move |n| {
            let a = 0.1; /* 振幅 */
            let f0 = 500.0; /* 周波数 */
            a * (2.0 * PI * f0 * (n as f64) / (pcm_fs as f64)).sin()
        }),
    );

    wave_write_16bit_mono_safer3("ex2_1.wav", &pcm);
}

fn ex2_2() {
    let pcm_fs: usize = 44100; /* 標本化周波数 */
    let pcm_length: usize = pcm_fs * 2; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    sine_wave(
        &mut pcm,
        261.63,
        0.1,
        mult(pcm_fs, 0.00),
        mult(pcm_fs, 0.25),
    ); /* C4 */
    sine_wave(
        &mut pcm,
        293.66,
        0.1,
        mult(pcm_fs, 0.25),
        mult(pcm_fs, 0.25),
    ); /* D4 */
    sine_wave(
        &mut pcm,
        329.63,
        0.1,
        mult(pcm_fs, 0.50),
        mult(pcm_fs, 0.25),
    ); /* E4 */
    sine_wave(
        &mut pcm,
        349.23,
        0.1,
        mult(pcm_fs, 0.75),
        mult(pcm_fs, 0.25),
    ); /* F4 */
    sine_wave(
        &mut pcm,
        392.00,
        0.1,
        mult(pcm_fs, 1.00),
        mult(pcm_fs, 0.25),
    ); /* G4 */
    sine_wave(
        &mut pcm,
        440.00,
        0.1,
        mult(pcm_fs, 1.25),
        mult(pcm_fs, 0.25),
    ); /* A4 */
    sine_wave(
        &mut pcm,
        493.88,
        0.1,
        mult(pcm_fs, 1.50),
        mult(pcm_fs, 0.25),
    ); /* B4 */
    sine_wave(
        &mut pcm,
        523.25,
        0.1,
        mult(pcm_fs, 1.75),
        mult(pcm_fs, 0.25),
    ); /* C5 */

    wave_write_16bit_mono_safer3("ex2_2.wav", &pcm);
}

fn ex3_1() {
    let f0 = 500.0; /* 基本周波数 */

    let pcm_fs = 44100;
    let pcm_length = pcm_fs * 1;
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    /* ノコギリ波 */
    for i_ in 1..=44 {
        let i = i_ as f64;
        for (n, item) in pcm.s.iter_mut().enumerate() {
            *item += 1.0 / i * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).sin();
        }
    }

    pcm.mult_constant_gain(0.1);

    wave_write_16bit_mono_safer3("ex3_1.wav", &pcm);
}

fn ex3_2() {
    let f0 = 500.0; /* 基本周波数 */

    let pcm_fs = 44100;
    let pcm_length = pcm_fs * 1;
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    /* 矩形波 */
    for j in 0..22 {
        let i = (2 * j + 1) as f64;
        for (n, item) in pcm.s.iter_mut().enumerate() {
            *item += 1.0 / i * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).sin();
        }
    }

    pcm.mult_constant_gain(0.1);

    wave_write_16bit_mono_safer3("ex3_2.wav", &pcm);
}

fn ex3_3() {
    let f0 = 500.0; /* 基本周波数 */

    let pcm_fs = 44100;
    let pcm_length = pcm_fs * 1;
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    /* 三角波 */
    for j in 0..22 {
        let i = (2 * j + 1) as f64;
        for (n, item) in pcm.s.iter_mut().enumerate() {
            *item += 1.0 / i / i * (PI * i / 2.0).sin()
                * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).sin();
        }
    }

    pcm.mult_constant_gain(0.1);

    wave_write_16bit_mono_safer3("ex3_3.wav", &pcm);
}

fn ex3_4() {
    let f0 = 500.0; /* 基本周波数 */

    let pcm_fs = 44100;
    let pcm_length = pcm_fs * 1;
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    /* コサイン波の重ね合わせによるノコギリ波 */
    for i in 1..=44 {
        let i = i as f64;
        for (n, item) in pcm.s.iter_mut().enumerate() {
            *item += 1.0 / i * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).cos();
        }
    }

    pcm.mult_constant_gain(0.1);

    wave_write_16bit_mono_safer3("ex3_4.wav", &pcm);
}

// slow and random; omitted from the test
fn ex3_5() {
    let f0 = 1.0; /* 基本周波数 */

    let pcm_fs = 44100;
    let pcm_length = pcm_fs * 1;
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let mut rng = rand::thread_rng();
    /* 白色雑音 */
    for i in 1..=22050 {
        let theta: f64 = rng.gen_range(0.0, 2.0 * PI);
        if i % 441 == 0 {
            println!("{} / 22050", i);
        }
        let i = i as f64;
        for (n, item) in pcm.s.iter_mut().enumerate() {
            *item += (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64) + theta).sin();
        }
    }

    pcm.mult_constant_gain(0.001);

    wave_write_16bit_mono_safer3("ex3_5.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex4_1() {
    let N = 64;
    let X = dft(N, Box::new(|_| 1.0));
    verify_(X)
}

#[allow(non_snake_case)]
fn ex4_2() {
    let N = 64;
    let w = create_Hanning_window(N); /* ハニング窓 */

    let X = dft(N, Box::new(move |n| w[n]));

    for (k, item) in X.iter().enumerate() {
        assert_close(item.re, 0.0);
        assert_close(
            item.im,
            match k {
                3 => 4.0,
                4 => -8.0,
                5 => 4.0,
                59 => -4.0,
                60 => 8.0,
                61 => -4.0,
                _ => 0.0,
            },
        );
    }
}

#[allow(non_snake_case)]
fn ex4_3() {
    let N = 64;
    let pcm_s = wave_read_16bit_mono_safer3("sine_500hz.wav").s;

    let mut x: Vec<Complex<f64>> = (0..N)
        .map(|n| {
            /* 波形 */
            Complex::new(pcm_s[n], 0.0) /* x(n)の実数部と虚数部 */
        })
        .collect();

    safe_FFT_(&mut x); /* FFTの計算結果はxに上書きされる */
}

#[allow(non_snake_case)]
fn ex4_4() {
    {
        let pcm0_fs = 44100;
        let pcm0_length = pcm0_fs * 1; /* 音データの長さ */
        let mut pcm0 = MonoPcm::new16(pcm0_fs, pcm0_length); /* 音データ */

        let f0 = 500.0; /* 基本周波数 */

        /* 基本音を含む音 */
        for i in 1..=44 {
            for (n, item) in pcm0.s.iter_mut().enumerate() {
                *item += (2.0 * PI * i as f64 * f0 * n as f64 / pcm0_fs as f64).sin();
            }
        }

        pcm0.mult_constant_gain(0.01);
        wave_write_16bit_mono_safer3("ex4_4a.wav", &pcm0);
    }
    {
        let pcm1_fs = 44100;
        let pcm1_length = pcm1_fs * 1; /* 音データの長さ */
        let mut pcm1 = MonoPcm::new16(pcm1_fs, pcm1_length); /* 音データ */
        let f0 = 500.0; /* 基本周波数 */
        /* 基本音を含まない音 */
        for i in 2..=44 {
            for (n, item) in pcm1.s.iter_mut().enumerate() {
                *item += (2.0 * PI * i as f64 * f0 * n as f64 / pcm1_fs as f64).sin();
            }
        }
        pcm1.mult_constant_gain(0.01);
        wave_write_16bit_mono_safer3("ex4_4b.wav", &pcm1);
    }
}

#[allow(non_snake_case)]
fn ex5_1() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 4; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length); /* 音データ */

    /* 振幅の時間エンベロープ */
    let initial_v = 0.5; /* a0[0] = 0.5; */
    let final_v = 0.0; /* a0[pcm_length - 1] = 0.0; */
    let a0 = linear(initial_v, final_v, pcm_length);

    let f0 = 500.0; /* 周波数 */
    for n in 0..pcm_length {
        pcm.s[n] = a0[n] * (2.0 * PI * f0 * n as f64 / pcm_fs as f64).sin();
    }

    wave_write_16bit_mono_safer3("ex5_1.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex5_2() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 4; /* 音データの長さ */

    let a0 = 0.5; /* 振幅 */
    /* 周波数の時間エンベロープ */
    let f0 = linear(2500.0, 1500.0, pcm_length);
    let g0: Vec<f64> = (0..pcm_length)
        .map(|n| {
            f0[0] * n as f64
                + (f0[pcm_length - 1] - f0[0]) * n as f64 * n as f64 / (pcm_length - 1) as f64 / 2.0
        })
        .collect();

    wave_write_16bit_mono_safer3(
        "ex5_2.wav",
        &MonoPcm::new16_fn(
            pcm_fs,
            pcm_length,
            Box::new(move |n| a0 * (2.0 * PI * g0[n] / pcm_fs as f64).sin()),
        ),
    );
}

fn ex5_3() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = (pcm_fs as f64 * 0.2) as usize; /* 音データの長さ */
    let a0 = 0.5; /* 振幅 */

    /* 周波数の時間エンベロープ */
    let f0 = linear(2500.0, 1500.0, pcm_length);
    let mut g0 = vec![0.0; pcm_length];

    for n in 0..pcm_length {
        g0[n] = f0[0] * n as f64
            + (f0[pcm_length - 1] - f0[0]) * n as f64 * n as f64 / (pcm_length - 1) as f64 / 2.0;
    }

    wave_write_16bit_mono_safer3(
        "ex5_3.wav",
        &MonoPcm::new16_fn(
            pcm_fs,
            pcm_length,
            Box::new(move |n| a0 * (2.0 * PI * g0[n] / pcm_fs as f64).sin()),
        ),
    );
}

fn ex5_4() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 4; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let a0 = vec![0.5; pcm_length];
    let a1 = vec![1.0; pcm_length];
    let a2 = vec![0.7; pcm_length];
    let a3 = vec![0.5; pcm_length];
    let a4 = vec![0.3; pcm_length];
    let f0 = vec![440.0; pcm_length];
    let f1 = vec![880.0; pcm_length];
    let f2 = vec![1320.0; pcm_length];
    let f3 = vec![1760.0; pcm_length];
    let f4 = vec![2200.0; pcm_length];

    /* 加算合成 */
    for n in 0..pcm_length {
        pcm.s[n] += a0[n] * (2.0 * PI * f0[n] * n as f64 / pcm_fs as f64).sin();
        pcm.s[n] += a1[n] * (2.0 * PI * f1[n] * n as f64 / pcm_fs as f64).sin();
        pcm.s[n] += a2[n] * (2.0 * PI * f2[n] * n as f64 / pcm_fs as f64).sin();
        pcm.s[n] += a3[n] * (2.0 * PI * f3[n] * n as f64 / pcm_fs as f64).sin();
        pcm.s[n] += a4[n] * (2.0 * PI * f4[n] * n as f64 / pcm_fs as f64).sin();
    }

    pcm.mult_constant_gain(0.1);

    /* フェード処理 */
    for n in 0..(pcm_fs as f64 * 0.01).ceil() as usize {
        pcm.s[n] *= n as f64 / (pcm_fs as f64 * 0.01);
        pcm.s[pcm_length - n - 1] *= n as f64 / (pcm_fs as f64 * 0.01);
    }

    wave_write_16bit_mono_safer3("ex5_4.wav", &pcm);
}

fn ex5_5() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 4; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length); /* 音データ */

    /* 時間エンベロープ */
    let a0 = exponential_decay(pcm_fs, 1.0, 5.0, 4.0, pcm_length);
    let a1 = exponential_decay(pcm_fs, 0.8, 5.0, 2.0, pcm_length);
    let a2 = exponential_decay(pcm_fs, 0.6, 5.0, 1.0, pcm_length);
    let a3 = exponential_decay(pcm_fs, 0.5, 5.0, 0.5, pcm_length);
    let a4 = exponential_decay(pcm_fs, 0.4, 5.0, 0.2, pcm_length);

    let f0 = vec![440.0; pcm_length];
    let f1 = vec![880.0; pcm_length];
    let f2 = vec![1320.0; pcm_length];
    let f3 = vec![1760.0; pcm_length];
    let f4 = vec![2200.0; pcm_length];

    /* 加算合成 */
    for (n, item) in pcm.s.iter_mut().enumerate() {
        *item += a0[n] * (2.0 * PI * f0[n] * n as f64 / pcm_fs as f64).sin();
        *item += a1[n] * (2.0 * PI * f1[n] * n as f64 / pcm_fs as f64).sin();
        *item += a2[n] * (2.0 * PI * f2[n] * n as f64 / pcm_fs as f64).sin();
        *item += a3[n] * (2.0 * PI * f3[n] * n as f64 / pcm_fs as f64).sin();
        *item += a4[n] * (2.0 * PI * f4[n] * n as f64 / pcm_fs as f64).sin();
    }

    pcm.mult_constant_gain(0.1);

    /* フェード処理 */
    for n in 0..(pcm_fs as f64 * 0.01).ceil() as usize {
        pcm.s[n] *= n as f64 / (pcm_fs as f64 * 0.01);
        pcm.s[pcm_length - n - 1] *= n as f64 / (pcm_fs as f64 * 0.01);
    }

    wave_write_16bit_mono_safer3("ex5_5.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex6_1() {
    let pcm0 = wave_read_16bit_mono_safer3("sine_500hz_3500hz.wav");

    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    let fe = 1000.0 / pcm0.fs as f64; /* エッジ周波数 */
    let delta = 1000.0 / pcm0.fs as f64; /* 遷移帯域幅 */
    let J = determine_J(delta); /* 遅延器の数 */

    let mut b = vec![0.0; J + 1];
    let mut w = create_Hanning_window(J + 1); /* ハニング窓 */

    safe_FIR_LPF(fe, J, &mut b, &mut w); /* FIRフィルタの設計 */
    safe_FIR_filtering(&pcm0.s, &mut pcm1.s, pcm1.length, &mut b, J);
    wave_write_16bit_mono_safer3("ex6_1.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex6_2() {
    let pcm0 = wave_read_16bit_mono_safer3("sine_500hz_3500hz.wav");
    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    let fc = 1000.0 / pcm0.fs as f64; /* 遮断周波数 */
    let Q = 1.0 / 2.0f64.sqrt(); /* クオリティファクタ */
    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */

    let mut a = [0.0; 3];
    let mut b = [0.0; 3];

    safe_IIR_LPF(fc, Q, &mut a, &mut b); /* IIRフィルタの設計 */
    safe_IIR_filtering(&pcm0.s, &mut pcm1.s, pcm1.length, &a, &b, I, J);

    wave_write_16bit_mono_safer3("ex6_2.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex6_3() {
    let pcm0 = wave_read_16bit_mono_safer3("sine_500hz_3500hz.wav");
    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    let fe = 1000.0 / pcm0.fs as f64; /* エッジ周波数 */
    let delta = 1000.0 / pcm0.fs as f64; /* 遷移帯域幅 */

    let J = determine_J(delta); /* 遅延器の数 */

    let mut b = vec![0.0; J + 1];
    let mut w = create_Hanning_window(J + 1); /* ハニング窓 */
    safe_FIR_LPF(fe, J, &mut b, &mut w); /* FIRフィルタの設計 */

    let L: usize = 128; /* フレームの長さ */
    let N = 256; /* DFTのサイズ */

    let number_of_frame = pcm0.length as usize / L; /* フレームの数 */
    for frame in 0..number_of_frame {
        let offset = (L * frame) as usize;
        /* X(k) */
        let mut x = vec![Complex::new(0.0, 0.0); N];

        for n in 0..L {
            x[n].re = pcm0.s[offset + n];
        }
        safe_FFT_(&mut x);

        /* B(k) */
        let mut b_ = vec![Complex::new(0.0, 0.0); N];
        for m in 0..=J {
            b_[m].re = b[m];
        }
        safe_FFT_(&mut b_);

        /* フィルタリング */
        let mut y: Vec<_> = (0..N).map(|k| x[k] * b_[k]).collect();
        safe_IFFT_(&mut y);

        /* オーバーラップアド */
        for n in 0..(L * 2) {
            if offset + n < pcm1.length as usize {
                pcm1.s[offset + n] += y[n].re;
            }
        }
    }
    wave_write_16bit_mono_safer3("ex6_3.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex6_4() {
    let pcm0 = wave_read_16bit_mono_safer3("sine_500hz_3500hz.wav");

    let N = 256; /* DFTのサイズ */

    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    let mut b_ = vec![Complex::new(0.0, 0.0); N];

    let w: Vec<f64> = create_Hanning_window(N); /* ハニング窓 */

    let number_of_frame = (pcm0.length as usize - N / 2) / (N / 2); /* フレームの数 */

    for frame in 0..number_of_frame {
        let offset = N / 2 * frame;

        /* X(n) */
        let mut x: Vec<_> = (0..N)
            .map(|n| Complex::new(pcm0.s[offset + n] * w[n], 0.0))
            .collect();
        safe_FFT_(&mut x);

        /* B(k) */
        let fe = 1000.0 / pcm0.fs as f64; /* エッジ周波数 */
        let fe = (fe * N as f64) as usize;
        for k in 0..=fe {
            b_[k] = Complex::new(1.0, 0.0);
        }
        for k in (fe + 1)..=N / 2 {
            b_[k] = Complex::new(0.0, 0.0);
        }
        for k in 1..N / 2 {
            b_[N - k] = b_[k].conj();
        }

        /* フィルタリング */
        let mut y: Vec<_> = (0..N).map(|k| x[k] * b_[k]).collect();
        safe_IFFT_(&mut y);

        /* オーバーラップアド */
        for (n, item) in y.iter().enumerate() {
            pcm1.s[offset + n] += item.re;
        }
    }
    wave_write_16bit_mono_safer3("ex6_4.wav", &pcm1);
}

// slooooow; omitted from the test
#[allow(non_snake_case)]
fn ex6_5() {
    let pcm0 = wave_read_16bit_mono_safer3("drum.wav");
    let pcm1 = wave_read_16bit_mono_safer3("response.wav");

    let J = pcm1.fs; /* 遅延器の数 */

    /* フィルタリング */
    let pcm2 = MonoPcm::new16_fn(
        pcm0.fs,
        pcm0.length,
        Box::new(move |n| {
            let mut a = 0.0;
            for m in 0..=J as usize {
                if n >= m {
                    a += pcm1.s[m] * pcm0.s[n - m];
                }
            }
            if n % 1000 == 0 {
                println!("{} / {}", n, pcm0.length);
            }
            a
        }),
    );

    wave_write_16bit_mono_safer3("ex6_5.wav", &pcm2);
}
