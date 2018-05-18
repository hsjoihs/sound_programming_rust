extern crate num_complex;
extern crate rand;
extern crate sound_programming;
//use std::io::Write;
use num_complex::Complex;
use rand::Rng;
use sound_programming::MonoPcm;
use sound_programming::StereoPcm;
use sound_programming::c_double;
use sound_programming::c_int;
use sound_programming::fft::safe_FFT_;
use sound_programming::fft::safe_IFFT;
use sound_programming::filter::safe_FIR_LPF;
use sound_programming::filter::safe_FIR_filtering;
use sound_programming::filter::safe_IIR_LPF;
use sound_programming::filter::safe_IIR_filtering;
use sound_programming::filter::safe_IIR_resonator;
use sound_programming::safe_ADSR;
use sound_programming::safe_Hanning_window;
use sound_programming::wave::wave_read_16bit_mono_safer3;
use sound_programming::wave::wave_read_16bit_stereo_safer3;
use sound_programming::wave_read_IMA_ADPCM_mono_safer3;
use sound_programming::wave_read_PCMA_mono_safer3;
use sound_programming::wave_read_PCMU_mono_safer3;
use sound_programming::wave_write_16bit_mono_safer3;
use sound_programming::wave_write_16bit_stereo_safer3;
use sound_programming::wave_write_IMA_ADPCM_mono_safer3;
use sound_programming::wave_write_PCMA_mono_safer3;
use sound_programming::wave_write_PCMU_mono_safer3;
use std::f64::consts::PI;
//use std::io;
fn main() {
    ex1_1();
    ex1_2();
    ex2_1();
    ex2_2();
    ex3_1();
    ex3_2();
    ex3_3();
    ex3_4();
    if false {
        ex3_5(); //slow
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
    ex7_1();
    ex7_2();
    ex7_3();
    /*if false*/ {
        ex7_4(); // slow
    }
    ex8_1();
    ex8_2();
    ex8_3();
    if false {
        ex8_4(); // random
        ex8_5(); // random
        ex8_6(); // random
    }
    ex8_7();
    ex8_8();
    ex8_9();
    ex8_10();
    ex8_11();
    ex8_12();
    ex10_4();
    ex11_7();
    ex11_8();
    ex11_9();
}

fn ex1_1() {
    /* 音データの入力 */

    let pcm0 = wave_read_16bit_mono_safer3("ex1_1_a.wav");
    let pcm1 = MonoPcm {
        s: (0..pcm0.length)
	      .map(|n| pcm0.s[n as usize])/* 音データのコピー */
	      .collect(),
        ..pcm0
    };

    wave_write_16bit_mono_safer3("ex1_1_b.wav", &pcm1); /* 音データの出力 */
}

#[allow(non_snake_case)]
fn ex1_2() {
    let pcm0 = wave_read_16bit_stereo_safer3("ex1_2_a.wav");
    let pcm1 = StereoPcm {
        s_l: (0..pcm0.length).map(|n| pcm0.s_l[n as usize]).collect(),
        s_r: (0..pcm0.length).map(|n| pcm0.s_r[n as usize]).collect(),
        ..pcm0
    };

    wave_write_16bit_stereo_safer3("ex1_2_b.wav", &pcm1);
}

fn ex2_1() {
    let pcm_fs: usize = 44100; /* 標本化周波数 */
    let pcm_length: usize = pcm_fs * 1; /* 音データの長さ */

    let a = 0.1; /* 振幅 */
    let f0 = 500.0; /* 周波数 */

    /* サイン波 */
    let pcm = MonoPcm::new16_fn(
        pcm_fs,
        pcm_length,
        Box::new(move |n| a * (2.0 * PI * f0 * (n as f64) / (pcm_fs as f64)).sin()),
    );

    wave_write_16bit_mono_safer3("ex2_1.wav", &pcm);
}

fn sine_wave(pcm: &mut MonoPcm, f0: c_double, a: c_double, offset: c_int, duration: c_int) {
    /* サイン波 */
    let mut s: Vec<c_double> = (0..duration)
        .map(|n| (2.0 * PI * f0 * (n as f64) / (pcm.fs as f64)).sin() * a)
        .collect();

    /* フェード処理 */
    for n in 0..(pcm.fs as f64 * 0.01).ceil() as usize {
        s[n] *= n as c_double / (pcm.fs as f64 * 0.01);
        s[duration as usize - n - 1] *= n as c_double / (pcm.fs as f64 * 0.01);
    }

    for n in 0..duration as usize {
        pcm.s[offset as usize + n] += s[n];
    }
}
fn ex2_2() {
    let pcm_fs: usize = 44100; /* 標本化周波数 */
    let pcm_length: usize = pcm_fs * 2; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    sine_wave(
        &mut pcm,
        261.63,
        0.1,
        itdyi(pcm_fs, 0.00),
        itdyi(pcm_fs, 0.25),
    ); /* C4 */
    sine_wave(
        &mut pcm,
        293.66,
        0.1,
        itdyi(pcm_fs, 0.25),
        itdyi(pcm_fs, 0.25),
    ); /* D4 */
    sine_wave(
        &mut pcm,
        329.63,
        0.1,
        itdyi(pcm_fs, 0.50),
        itdyi(pcm_fs, 0.25),
    ); /* E4 */
    sine_wave(
        &mut pcm,
        349.23,
        0.1,
        itdyi(pcm_fs, 0.75),
        itdyi(pcm_fs, 0.25),
    ); /* F4 */
    sine_wave(
        &mut pcm,
        392.00,
        0.1,
        itdyi(pcm_fs, 1.00),
        itdyi(pcm_fs, 0.25),
    ); /* G4 */
    sine_wave(
        &mut pcm,
        440.00,
        0.1,
        itdyi(pcm_fs, 1.25),
        itdyi(pcm_fs, 0.25),
    ); /* A4 */
    sine_wave(
        &mut pcm,
        493.88,
        0.1,
        itdyi(pcm_fs, 1.50),
        itdyi(pcm_fs, 0.25),
    ); /* B4 */
    sine_wave(
        &mut pcm,
        523.25,
        0.1,
        itdyi(pcm_fs, 1.75),
        itdyi(pcm_fs, 0.25),
    ); /* C5 */

    wave_write_16bit_mono_safer3("ex2_2.wav", &pcm);
}

// int_times_double_yielding_int
fn itdyi(i: usize, d: c_double) -> c_int {
    ((i as c_double) * d) as c_int
}

fn ex3_1() {
    let f0 = 500.0; /* 基本周波数 */

    let pcm_fs = 44100;
    let pcm_length = pcm_fs * 1;
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    /* ノコギリ波 */
    for i_ in 1..=44 {
        let i = i_ as f64;
        for n in 0..pcm_length {
            pcm.s[n] += 1.0 / i * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).sin();
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
        for n in 0..pcm_length {
            pcm.s[n] += 1.0 / i * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).sin();
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
        for n in 0..pcm_length {
            pcm.s[n] += 1.0 / i / i * (PI * i / 2.0).sin()
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
        for n in 0..pcm_length {
            pcm.s[n] += 1.0 / i * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).cos();
        }
    }

    pcm.mult_constant_gain(0.1);

    wave_write_16bit_mono_safer3("ex3_4.wav", &pcm);
}

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
        for n in 0..pcm_length {
            pcm.s[n] += (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64) + theta).sin();
        }
    }

    pcm.mult_constant_gain(0.001);

    wave_write_16bit_mono_safer3("ex3_5.wav", &pcm);
}

#[allow(non_snake_case)]
fn verify_(X: Vec<Complex<f64>>) {
    let N = 64;

    /* 周波数特性 */
    for k in 0..N {
        assert_close(X[k].re, 0.0);
        assert_close(
            X[k].im,
            match k {
                4 => -16.0,
                60 => 16.0,
                _ => 0.0,
            },
        );
    }
}

#[allow(non_snake_case)]
fn ex4_1() {
    let X = foo_(Box::new(|_| 1.0));
    verify_(X)
}

#[allow(non_snake_case)]
fn foo_(func: Box<Fn(usize) -> f64>) -> Vec<Complex<f64>> {
    let N = 64;
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
            let W = Complex::new((2.0 * PI * k * n / N).cos(), -(2.0 * PI * k * n / N).sin());
            X[k_] += W * x[n_];
        }
    }

    X
}

fn assert_close(a: f64, b: f64) {
    assert!((a - b).abs() < 1.75e-4);
}

#[allow(non_snake_case)]
fn ex4_2() {
    let N = 64;
    let mut w: Vec<c_double> = vec![0.0; N];
    safe_Hanning_window(&mut w); /* ハニング窓 */

    let X = foo_(Box::new(move |n| w[n]));

    for k in 0..N {
        assert_close(X[k].re, 0.0);
        assert_close(
            X[k].im,
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
    let mut x: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); N];
    let pcm_s = wave_read_16bit_mono_safer3("sine_500hz.wav").s;

    /* 波形 */
    for n in 0..N {
        x[n] = Complex::new(pcm_s[n], 0.0); /* x(n)の虚数部 */
    }

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
            for n in 0..pcm0_length {
                pcm0.s[n] += (2.0 * PI * i as f64 * f0 * n as f64 / pcm0_fs as f64).sin();
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
            for n in 0..pcm1_length {
                pcm1.s[n] += (2.0 * PI * i as f64 * f0 * n as f64 / pcm1_fs as f64).sin();
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

    let mut a0 = vec![0.0; pcm_length];
    /* 振幅の時間エンベロープ */
    a0[0] = 0.5;
    a0[pcm_length - 1] = 0.0;
    for n in 0..pcm_length {
        a0[n] = a0[0] + (a0[pcm_length - 1] - a0[0]) * n as f64 / (pcm_length - 1) as f64;
    }
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
    let mut f0 = vec![0.0; pcm_length];
    let mut g0 = vec![0.0; pcm_length];
    /* 周波数の時間エンベロープ */
    f0[0] = 2500.0;
    f0[pcm_length - 1] = 1500.0;
    for n in 0..pcm_length {
        f0[n] = f0[0] + (f0[pcm_length - 1] - f0[0]) * n as f64 / (pcm_length - 1) as f64;
    }
    for n in 0..pcm_length {
        g0[n] = f0[0] * n as f64
            + (f0[pcm_length - 1] - f0[0]) * n as f64 * n as f64 / (pcm_length - 1) as f64 / 2.0;
    }

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
    let mut f0 = vec![0.0; pcm_length];
    let mut g0 = vec![0.0; pcm_length];
    /* 周波数の時間エンベロープ */
    f0[0] = 2500.0;
    f0[pcm_length - 1] = 1500.0;
    for n in 0..pcm_length {
        f0[n] = f0[0] + (f0[pcm_length - 1] - f0[0]) * n as f64 / (pcm_length - 1) as f64;
    }
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

    let mut a0 = vec![0.0; pcm_length];
    let mut a1 = vec![0.0; pcm_length];
    let mut a2 = vec![0.0; pcm_length];
    let mut a3 = vec![0.0; pcm_length];
    let mut a4 = vec![0.0; pcm_length];
    let f0 = vec![440.0; pcm_length];
    let f1 = vec![880.0; pcm_length];
    let f2 = vec![1320.0; pcm_length];
    let f3 = vec![1760.0; pcm_length];
    let f4 = vec![2200.0; pcm_length];
    /* 時間エンベロープ */
    for n in 0..pcm_length {
        a0[n] = 1.0 * (-5.0 * n as f64 / (pcm_fs as f64 * 4.0)).exp();
        a1[n] = 0.8 * (-5.0 * n as f64 / (pcm_fs as f64 * 2.0)).exp();
        a2[n] = 0.6 * (-5.0 * n as f64 / (pcm_fs as f64 * 1.0)).exp();
        a3[n] = 0.5 * (-5.0 * n as f64 / (pcm_fs as f64 * 0.5)).exp();
        a4[n] = 0.4 * (-5.0 * n as f64 / (pcm_fs as f64 * 0.2)).exp();
    }
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

    wave_write_16bit_mono_safer3("ex5_5.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex6_1() {
    let pcm0 = wave_read_16bit_mono_safer3("sine_500hz_3500hz.wav");

    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    let fe = 1000.0 / pcm0.fs as f64; /* エッジ周波数 */
    let delta = 1000.0 / pcm0.fs as f64; /* 遷移帯域幅 */

    let mut J = (3.1 / delta + 0.5) as usize - 1; /* 遅延器の数 */
    if J % 2 == 1 {
        J += 1; /* J+1が奇数になるように調整する */
    }

    let mut b = vec![0.0; J + 1];
    let mut w = vec![0.0; J + 1];
    safe_Hanning_window(&mut w); /* ハニング窓 */

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

    let mut J = (3.1 / delta + 0.5) as usize - 1; /* 遅延器の数 */
    if J % 2 == 1 {
        J += 1; /* J+1が奇数になるように調整する */
    }

    let mut b = vec![0.0; J + 1];
    let mut w = vec![0.0; J + 1];
    safe_Hanning_window(&mut w); /* ハニング窓 */
    safe_FIR_LPF(fe, J, &mut b, &mut w); /* FIRフィルタの設計 */

    let L: usize = 128; /* フレームの長さ */
    let N = 256; /* DFTのサイズ */

    let mut y_real = vec![0.0; N];
    let mut y_imag = vec![0.0; N];
    let mut b_ = vec![Complex::new(0.0, 0.0); N];

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
        for m in 0..N {
            b_[m].re = 0.0;
            b_[m].im = 0.0;
        }
        for m in 0..=J {
            b_[m].re = b[m];
        }
        safe_FFT_(&mut b_);

        /* フィルタリング */
        for k in 0..N {
            y_real[k] = x[k].re * b_[k].re - x[k].im * b_[k].im;
            y_imag[k] = x[k].im * b_[k].re + x[k].re * b_[k].im;
        }
        safe_IFFT(&mut y_real, &mut y_imag);

        /* オーバーラップアド */
        for n in 0..(L * 2) {
            if offset + n < pcm1.length as usize {
                pcm1.s[offset + n] += y_real[n];
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

    let mut x = vec![Complex::new(0.0, 0.0); N];
    let mut y_real: Vec<c_double> = vec![0.0; N];
    let mut y_imag: Vec<c_double> = vec![0.0; N];
    let mut b_real: Vec<c_double> = vec![0.0; N];
    let mut b_imag: Vec<c_double> = vec![0.0; N];

    let mut w: Vec<c_double> = vec![0.0; N];
    safe_Hanning_window(&mut w); /* ハニング窓 */

    let number_of_frame = (pcm0.length as usize - N / 2) / (N / 2); /* フレームの数 */

    for frame in 0..number_of_frame {
        let offset = N / 2 * frame;

        /* X(n) */
        for n in 0..N {
            x[n] = Complex::new(pcm0.s[offset + n] * w[n], 0.0);
        }
        safe_FFT_(&mut x);

        /* B(k) */
        let fe = 1000.0 / pcm0.fs as f64; /* エッジ周波数 */
        let fe = (fe * N as f64) as usize;
        for k in 0..=fe {
            b_real[k] = 1.0;
            b_imag[k] = 0.0;
        }
        for k in (fe + 1)..=N / 2 {
            b_real[k] = 0.0;
            b_imag[k] = 0.0;
        }
        for k in 1..N / 2 {
            b_real[N - k] = b_real[k];
            b_imag[N - k] = -b_imag[k];
        }

        /* フィルタリング */
        for k in 0..N {
            y_real[k] = x[k].re * b_real[k] - x[k].im * b_imag[k];
            y_imag[k] = x[k].im * b_real[k] + x[k].re * b_imag[k];
        }
        safe_IFFT(&mut y_real, &mut y_imag);

        /* オーバーラップアド */
        for n in 0..N {
            pcm1.s[offset + n] += y_real[n];
        }
    }
    wave_write_16bit_mono_safer3("ex6_4.wav", &pcm1);
}

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

#[allow(non_snake_case)]
fn ex7_1() {
    let pcm0 = wave_read_16bit_mono_safer3("ex7_1_pulse_train.wav");

    let mut fc = vec![0.0; pcm0.length as usize];
    /* LPFの遮断周波数 */
    for n in 0..pcm0.length as usize {
        fc[n] = 10000.0 * (-5.0 * n as f64 / pcm0.length as f64).exp();
    }
    let Q = 1.0 / 2.0f64.sqrt(); /* クオリティファクタ */
    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */

    let pcm1_length = pcm0.length; /* 音データの長さ */
    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    let mut a = [0.0; 3];
    let mut b = [0.0; 3];

    for n in 0..pcm1_length as usize {
        safe_IIR_LPF(fc[n] / pcm1.fs as f64, Q, &mut a, &mut b); /* IIRフィルタの設計 */

        for m in 0..=J {
            if n >= m {
                pcm1.s[n] += b[m] * pcm0.s[n - m];
            }
        }
        for m in 1..=I {
            if n >= m {
                pcm1.s[n] += -a[m] * pcm1.s[n - m];
            }
        }
    }
    wave_write_16bit_mono_safer3("ex7_1.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex7_2() {
    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    let pcm0 = wave_read_16bit_mono_safer3("white_noise.wav");
    let mut fc = vec![0.0; pcm0.length as usize];
    /* LPFの遮断周波数 */
    for n in 0..pcm0.length as usize {
        fc[n] = 10000.0 * (-5.0 * n as f64 / pcm0.length as f64).exp();
    }
    let Q = 1.0 / 2.0f64.sqrt(); /* クオリティファクタ */
    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */

    let mut pcm1 = MonoPcm::blank_copy(&pcm0); /* 音データ */

    for n in 0..pcm1.length as usize {
        safe_IIR_LPF(fc[n] / pcm1.fs as f64, Q, &mut a, &mut b); /* IIRフィルタの設計 */
        for m in 0..=J {
            if n >= m {
                pcm1.s[n] += b[m] * pcm0.s[n - m];
            }
        }
        for m in 1..=I {
            if n >= m {
                pcm1.s[n] += -a[m] * pcm1.s[n - m];
            }
        }
    }
    wave_write_16bit_mono_safer3("ex7_2.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex7_3() {
    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    let pcm0 = wave_read_16bit_mono_safer3("ex7_3_pulse_train.wav");
    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    let mut s = vec![0.0; pcm1.length as usize];
    let F = [
        800.0 /* F1の周波数 */, 1200.0 /* F2の周波数 */,
        2500.0 /* F3の周波数 */, 3500.0,
    ]; /* F4の周波数 */

    let B = [
        100.0 /* F1の帯域幅 */, 100.0 /* F2の帯域幅 */,
        100.0 /* F3の帯域幅 */, 100.0,
    ]; /* F4の帯域幅 */

    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */

    for num in 0..4 {
        /* IIRフィルタの設計 */
        safe_IIR_resonator(F[num] / pcm0.fs as f64, F[num] / B[num], &mut a, &mut b);
        safe_IIR_filtering(&pcm0.s, &mut s, pcm0.length, &a, &b, I, J); /* フィルタリング */
        for n in 0..pcm1.length as usize {
            pcm1.s[n] += s[n];
            s[n] = 0.0;
        }
    }

    /* ディエンファシス処理 */
    s[0] = pcm1.s[0];
    for n in 1..pcm1.length as usize {
        s[n] = pcm1.s[n] + 0.98 * s[n - 1];
    }
    for n in 0..pcm1.length as usize {
        pcm1.s[n] = s[n];
    }
    wave_write_16bit_mono_safer3("ex7_3.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex7_4() {
    let pcm0 = wave_read_16bit_mono_safer3("synth.wav");
    let mut pcm1 = wave_read_16bit_mono_safer3("vocal.wav");
    let mut pcm2 = MonoPcm::blank_copy(&pcm0);

    let mut s = vec![0.0; pcm0.length]; /* 音データ */
    /* プリエンファシス処理 */
    s[0] = 0.0;
    for n in 1..pcm1.length {
        s[n] = pcm1.s[n] - 0.98 * pcm1.s[n - 1];
    }
    for n in 0..pcm1.length {
        pcm1.s[n] = s[n];
    }

    let N = 1024; /* DFTのサイズ */

    let mut x = vec![Complex::new(0.0, 0.0); N];
    let mut y_real = vec![0.0; N];
    let mut y_imag = vec![0.0; N];
    let mut b_ = vec![Complex::new(0.0, 0.0); N];
    let mut w = vec![0.0; N];
    safe_Hanning_window(&mut w); /* ハニング窓 */
    let number_of_frame = (pcm0.length - N / 2) / (N / 2); /* フレームの数 */
    let band_width = 8;
    let number_of_band = N / 2 / band_width;

    for frame in 0..number_of_frame {
        let offset = N / 2 * frame;
        /* X(n) */
        for n in 0..N {
            x[n] = Complex::new(pcm0.s[offset + n] * w[n], 0.0);
        }
        safe_FFT_(&mut x);

        /* B(k) */
        for n in 0..N {
            b_[n].re = pcm1.s[offset + n] * w[n];
            b_[n].im = 0.0;
        }
        safe_FFT_(&mut b_);

        for k in 0..N {
            b_[k].re = b_[k].norm_sqr().sqrt();
            b_[k].im = 0.0;
        }

        for band in 0..number_of_band {
            let offset = band_width * band;
            let mut a = 0.0;
            for k in 0..band_width {
                a += b_[offset + k].re;
            }
            a /= band_width as f64;
            for k in 0..band_width {
                b_[offset + k].re = a;
            }
        }
        b_[0].re = 0.0;
        b_[N / 2].re = 0.0;
        for k in 1..N / 2 {
            b_[N - k].re = b_[k].re;
        }

        /* フィルタリング */

        for k in 0..N {
            y_real[k] = x[k].re * b_[k].re - x[k].im * b_[k].im;
            y_imag[k] = x[k].im * b_[k].re + x[k].re * b_[k].im;
        }
        safe_IFFT(&mut y_real, &mut y_imag);

        let offset = N / 2 * frame;
        /* オーバーラップアド */
        for n in 0..N {
            pcm2.s[offset + n] += y_real[n];
        }
    }
    wave_write_16bit_mono_safer3("ex7_4.wav", &pcm2);
}

#[allow(non_snake_case)]
fn ex8_1() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 1; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let f0 = 500.0; /* 基本周波数 */
    /* ノコギリ波 */
    let t0 = pcm_fs / f0 as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..pcm_length {
        pcm.s[n] = 1.0 - 2.0 * (m as f64) / (t0 as f64);

        m += 1;
        if m >= t0 {
            m = 0;
        }
    }
    pcm.mult_constant_gain(0.1);

    wave_write_16bit_mono_safer3("ex8_1.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex8_2() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 1; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let f0 = 500.0; /* 基本周波数 */

    /* 矩形波 */
    let t0 = (pcm.fs as f64 / f0) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..pcm.length {
        pcm.s[n] = if (m as f64) < t0 as f64 / 2.0 {
            1.0
        } else {
            -1.0
        };

        m += 1;
        if m >= t0 {
            m = 0;
        }
    }
    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex8_2.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex8_3() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 1; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let f0 = 500.0; /* 基本周波数 */
    /* 三角波 */
    let t0 = (pcm_fs as f64 / f0) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..pcm_length {
        pcm.s[n] = if (m as f64) < t0 as f64 / 2.0 {
            -1.0 + 4.0 * m as f64 / t0 as f64
        } else {
            3.0 - 4.0 * m as f64 / t0 as f64
        };

        m += 1;
        if m >= t0 {
            m = 0;
        }
    }
    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex8_3.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex8_4() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 1; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let mut rng = rand::thread_rng();
    /* 白色雑音 */
    for n in 0..pcm_length {
        pcm.s[n] = rng.gen_range(-1.0, 1.0);
    }
    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex8_4.wav", &pcm);
}

fn ex8_5() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 8; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let mut rng = rand::thread_rng();

    /* 白色雑音 */
    for n in 0..pcm_length {
        pcm.s[n] = rng.gen_range(-1.0, 1.0);
    }

    let mut e = vec![0.0; pcm_length];
    let te = pcm_fs * 2; /* 単調増加または単調減少にかかる時間 */

    /* 時間エンベロープ */
    let mut m = 0;
    for n in 0..pcm_length {
        e[n] = if m < te {
            m as f64 / te as f64
        } else {
            1.0 - (m as f64 - te as f64) / te as f64
        };

        m += 1;
        if m >= te * 2 {
            m = 0;
        }
    }
    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex8_5.wav", &pcm);
}

fn square_wave(
    x: (&mut [c_double], usize),
    f0: c_double,
    gain: c_double,
    offset: usize,
    duration: usize,
) {
    let mut s = vec![0.0; duration];
    let (pcm_s, pcm_fs) = x;
    /* 矩形波 */
    let t0 = (pcm_fs as f64 / f0) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..duration {
        s[n] = if (m as f64) < t0 as f64 / 2.0 {
            1.0
        } else {
            -1.0
        };

        m += 1;
        if m >= t0 {
            m = 0;
        }
    }

    for n in 0..duration {
        s[n] *= gain;
    }

    for n in 0..duration {
        pcm_s[offset + n] += s[n];
    }
}

fn triangle_wave(
    x: (&mut [c_double], usize),
    f0: c_double,
    gain: c_double,
    offset: usize,
    duration: usize,
) {
    let mut s = vec![0.0; duration];
    let (pcm_s, pcm_fs) = x;
    /* 三角波 */
    let t0 = (pcm_fs as f64 / f0) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..duration {
        s[n] = if (m as f64) < t0 as f64 / 2.0 {
            -1.0 + 4.0 * m as f64 / t0 as f64
        } else {
            3.0 - 4.0 * m as f64 / t0 as f64
        };
        m += 1;
        if m >= t0 {
            m = 0;
        }
    }
    for n in 0..duration {
        s[n] *= gain;
    }
    for n in 0..duration {
        pcm_s[offset + n] += s[n];
    }
}

fn white_noise(pcm_s: &mut [c_double], gain: c_double, offset: usize, duration: usize) {
    let mut s = vec![0.0; duration];
    let mut rng = rand::thread_rng();
    /* 白色雑音 */
    for n in 0..duration {
        s[n] = rng.gen_range(-1.0, 1.0);
    }
    for n in 0..duration {
        s[n] *= gain;
    }
    for n in 0..duration {
        pcm_s[offset + n] += s[n];
    }
}

fn ex8_6() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = 7000 * 16; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length); /* 音データ */

    /* メロディパート */
    square_wave((&mut pcm.s, pcm_fs), 659.26, 0.1, 7000 * 0, 6125); /* E5 */
    square_wave((&mut pcm.s, pcm_fs), 659.26, 0.1, 7000 * 1, 6125); /* E5 */
    square_wave((&mut pcm.s, pcm_fs), 659.26, 0.1, 7000 * 3, 6125); /* E5 */
    square_wave((&mut pcm.s, pcm_fs), 523.25, 0.1, 7000 * 5, 6125); /* C5 */
    square_wave((&mut pcm.s, pcm_fs), 659.26, 0.1, 7000 * 6, 6125); /* E5 */
    square_wave((&mut pcm.s, pcm_fs), 783.99, 0.1, 7000 * 8, 6125); /* G5 */
    square_wave((&mut pcm.s, pcm_fs), 392.00, 0.1, 7000 * 12, 6125); /* G4 */

    /* ベースパート */
    triangle_wave((&mut pcm.s, pcm_fs), 146.83, 0.2, 7000 * 0, 6125); /* D3 */
    triangle_wave((&mut pcm.s, pcm_fs), 146.83, 0.2, 7000 * 1, 6125); /* D3 */
    triangle_wave((&mut pcm.s, pcm_fs), 146.83, 0.2, 7000 * 3, 6125); /* D3 */
    triangle_wave((&mut pcm.s, pcm_fs), 146.83, 0.2, 7000 * 5, 6125); /* D3 */
    triangle_wave((&mut pcm.s, pcm_fs), 146.83, 0.2, 7000 * 6, 6125); /* D3 */
    triangle_wave((&mut pcm.s, pcm_fs), 196.00, 0.2, 7000 * 8, 6125); /* G3 */
    triangle_wave((&mut pcm.s, pcm_fs), 196.00, 0.2, 7000 * 12, 6125); /* G3 */

    /* パーカッション */
    white_noise(&mut pcm.s, 0.1, 7000 * 0, 4000);
    white_noise(&mut pcm.s, 0.1, 7000 * 2, 1000);
    white_noise(&mut pcm.s, 0.1, 7000 * 3, 4000);
    white_noise(&mut pcm.s, 0.1, 7000 * 5, 1000);
    white_noise(&mut pcm.s, 0.1, 7000 * 6, 4000);
    white_noise(&mut pcm.s, 0.1, 7000 * 8, 4000);
    white_noise(&mut pcm.s, 0.1, 7000 * 11, 4000);
    white_noise(&mut pcm.s, 0.1, 7000 * 13, 1000);
    white_noise(&mut pcm.s, 0.1, 7000 * 14, 1000);
    white_noise(&mut pcm.s, 0.1, 7000 * 15, 1000);

    wave_write_16bit_mono_safer3("ex8_6.wav", &pcm);
}

fn ex8_7() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = (pcm_fs as f64 * 0.6) as usize - 1; /* 音データの長さ */
    // -1 is here to match the file with the original data
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let mut f0 = vec![0.0; pcm_length];
    /* 基本周波数 */
    for n in 0..itdyi(pcm_fs, 0.1) as usize {
        f0[n] = 987.77; /* B5 */
    }
    for n in itdyi(pcm_fs, 0.1) as usize..pcm_length {
        f0[n] = 1318.51; /* E6 */
    }

    /* 矩形波 */
    let mut t0 = (pcm_fs as f64 / f0[0]) as usize; /* 矩形波の基本周期 */
    let mut m = 0;
    for n in 0..pcm_length {
        pcm.s[n] = if (m as f64) < t0 as f64 / 2.0 {
            1.0
        } else {
            -1.0
        };

        m += 1;
        if m >= t0 {
            t0 = (pcm_fs as f64 / f0[n]) as usize; /* 矩形波の基本周期 */
            m = 0;
        }
    }
    let mut e = vec![0.0; pcm_length];
    let pe = pcm_length; /* 単調減少にかかる時間 */

    /* 時間エンベロープ */
    for n in 0..pcm_length {
        e[n] = 1.0 - n as f64 / pe as f64;
    }
    let gain = 0.1; /* ゲイン */

    for n in 0..pcm_length {
        pcm.s[n] *= e[n] * gain;
    }

    wave_write_16bit_mono_safer3("ex8_7.wav", &pcm);
}

fn ex8_8() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = (pcm_fs as f64 * 0.6) as usize - 1; /* 音データの長さ */
    // -1 is here to match the file with the original data
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    /* 基本周波数 */
    let mut f0 = vec![0.0; pcm_length];
    for n in 0..(pcm_fs as f64 * 0.2) as usize {
        f0[n] = 440.0;
    }
    for n in (pcm_fs as f64 * 0.2) as usize..pcm_length {
        f0[n] = 440.0
            + (880.0 - 440.0) * (n as f64 - pcm_fs as f64 * 0.2)
                / (pcm_length as f64 - 1.0 - pcm_fs as f64 * 0.2);
    }
    /* 矩形波 */
    let mut t0 = (pcm_fs as f64 / f0[0]) as usize; /* 矩形波の基本周期 */
    let mut m = 0;
    for n in 0..pcm_length {
        pcm.s[n] = if (m as f64) < t0 as f64 / 2.0 {
            1.0
        } else {
            -1.0
        };
        m += 1;
        if m >= t0 {
            t0 = (pcm_fs as f64 / f0[n]) as usize; /* 矩形波の基本周期 */
            m = 0;
        }
    }

    let mut e = vec![0.0; pcm_length];
    let pe = pcm_length; /* 単調減少にかかる時間 */
    /* 時間エンベロープ */
    for n in 0..pcm_length {
        e[n] = 1.0 - n as f64 / pe as f64;
    }

    let gain = 0.1; /* ゲイン */
    for n in 0..pcm_length {
        pcm.s[n] *= e[n] * gain;
    }
    wave_write_16bit_mono_safer3("ex8_8.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex8_9() {
    let pcm_fs = 8000; /* 標本化周波数 */
    let pcm_length = pcm_fs * 2; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let mut f0 = vec![0.0; pcm_length];

    /* 基本周波数 */
    f0[0] = 500.0;
    f0[pcm_length - 1] = 3500.0;
    for n in 0..pcm_length {
        f0[n] = f0[0] + (f0[pcm_length - 1] - f0[0]) * n as f64 / (pcm_length - 1) as f64;
    }

    /* ノコギリ波 */
    let mut t0 = (pcm_fs as f64 / f0[0]) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..pcm_length {
        pcm.s[n] = 1.0 - 2.0 * m as f64 / t0 as f64;

        m += 1;
        if m >= t0 {
            t0 = (pcm_fs as f64 / f0[n]) as usize; /* 基本周期 */
            m = 0;
        }
    }
    pcm.mult_constant_gain(0.1);

    wave_write_16bit_mono_safer3("ex8_9.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex8_10() {
    let pcm0_fs = 192000; /* 標本化周波数 */
    let _pcm0_bits = 16; /* 量子化精度 */
    let pcm0_length = pcm0_fs * 2; /* 音データの長さ */
    let mut pcm0_s = vec![0.0; pcm0_length];
    let mut f0 = vec![0.0; pcm0_length];
    /* 基本周波数 */
    f0[0] = 500.0;
    f0[pcm0_length - 1] = 3500.0;
    for n in 0..pcm0_length {
        f0[n] = f0[0] + (f0[pcm0_length - 1] - f0[0]) * n as f64 / (pcm0_length - 1) as f64;
    }

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

    let pcm1_fs = 8000; /* 標本化周波数 */
    let pcm1_length = pcm1_fs * 2; /* 音データの長さ */
    let mut pcm1 = MonoPcm::new16(pcm1_fs, pcm1_length);
    let ratio = pcm0_fs / pcm1_fs; /* ダウンサンプリングのレシオ */
    let fe = 0.45 / ratio as f64; /* エッジ周波数 */
    let delta = 0.1 / ratio as f64; /* 遷移帯域幅 */
    let mut J = (3.1 / delta + 0.5) as usize - 1; /* 遅延器の数 */
    if J % 2 == 1 {
        J += 1; /* J+1が奇数になるように調整する */
    }
    let mut b = vec![0.0; J + 1];
    let mut w = vec![0.0; J + 1];

    safe_Hanning_window(&mut w); /* ハニング窓 */
    safe_FIR_LPF(fe, J, &mut b, &mut w); /* FIRフィルタの設計 */

    /* フィルタリング */
    for n in 0..pcm1_length {
        for m in 0..=J {
            if n * ratio + J / 2 >= m && n * ratio + J / 2 < pcm0_length + m {
                pcm1.s[n] += b[m] * pcm0_s[n * ratio + J / 2 - m];
            }
        }
    }
    pcm1.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex8_10.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex8_11() {
    let pcm0_fs = 192000; /* 標本化周波数 */
    let pcm0_length = pcm0_fs * 2; /* 音データの長さ */

    let mut pcm0 = MonoPcm::new16(pcm0_fs, pcm0_length);

    let mut f0 = vec![0.0; pcm0_length];

    /* 基本周波数 */
    f0[0] = 500.0;
    f0[pcm0_length - 1] = 3500.0;
    for n in 0..pcm0_length {
        f0[n] = f0[0] + (f0[pcm0_length - 1] - f0[0]) * n as f64 / (pcm0_length - 1) as f64;
    }
    /* ノコギリ波 */
    let mut t0 = (pcm0_fs as f64 / f0[0]) as usize; /* 基本周期 */

    let mut m = 0;
    for n in 0..pcm0_length {
        pcm0.s[n] = 1.0 - 2.0 * m as f64 / t0 as f64;

        m += 1;
        if m >= t0 {
            t0 = (pcm0_fs as f64 / f0[n]) as usize; /* 基本周期 */
            m = 0;
        }
    }
    let pcm1_fs = 8000; /* 標本化周波数 */
    let pcm1_length = pcm1_fs * 2; /* 音データの長さ */
    let mut pcm1 = MonoPcm::new16(pcm1_fs, pcm1_length);

    let ratio = pcm0_fs / pcm1_fs; /* ダウンサンプリングのレシオ */

    for n in 0..pcm1_length {
        pcm1.s[n] = pcm0.s[n * ratio];
    }

    pcm1.mult_constant_gain(0.1);

    wave_write_16bit_mono_safer3("ex8_11.wav", &pcm1);
}

#[allow(non_snake_case, unused_variables, unused_mut)]
fn ex8_12() {
    let mut pcm0 = MonoPcm::new16(192000, 192000 * 2);

    let mut f0 = vec![0.0; pcm0.length];
    /* 基本周波数 */
    f0[0] = 500.0;
    f0[pcm0.length - 1] = 3500.0;
    for n in 0..pcm0.length {
        f0[n] = f0[0] + (f0[pcm0.length - 1] - f0[0]) * n as f64 / (pcm0.length - 1) as f64;
    }
    {
        /* ノコギリ波 */
        let mut t0 = (pcm0.fs as f64 / f0[0]) as usize; /* 基本周期 */
        let mut m = 0;
        for n in 0..pcm0.length {
            pcm0.s[n] = 1.0 - 2.0 * m as f64 / t0 as f64;

            m += 1;
            if m >= t0 {
                t0 = (pcm0.fs as f64 / f0[n]) as usize; /* 基本周期 */
                m = 0;
            }
        }
    }
    let mut pcm1 = MonoPcm::new16(8000, 8000 * 2);
    let ratio = pcm0.fs / pcm1.fs; /* ダウンサンプリングのレシオ */
    let fe = 0.45 / ratio as f64; /* エッジ周波数 */
    let delta = 0.1 / ratio as f64; /* 遷移帯域幅 */
    let mut J = (3.1 / delta + 0.5) as usize - 1; /* 遅延器の数 */
    if J % 2 == 1 {
        J += 1; /* J+1が奇数になるように調整する */
    }
    let mut b = vec![0.0; J + 1];
    let mut w = vec![0.0; J + 1];
    safe_Hanning_window(&mut w); /* ハニング窓 */
    safe_FIR_LPF(fe, J, &mut b, &mut w); /* FIRフィルタの設計 */
    for n in 0..pcm1.length {
        for m in 0..=J {
            if n * ratio + J / 2 >= m && n * ratio + J / 2 < pcm0.length + m {
                pcm1.s[n] += b[m] * pcm0.s[n * ratio + J / 2 - m];
            }
        }
    }
    let mut pcm2 = MonoPcm::new16(192000, 192000 * 2);
    /* 0を挿入する */
    for n in 0..pcm1.length {
        pcm2.s[n * ratio] = pcm1.s[n];
    }

    let mut pcm3 = MonoPcm::new16(192000, 192000 * 2);
    for n in 0..pcm3.length {
        for m in 0..=J {
            if n + J / 2 >= m && n + J / 2 < pcm2.length + m {
                pcm3.s[n] += b[m] * pcm2.s[n + J / 2 - m];
            }
        }
    }

    pcm3.mult_constant_gain(ratio as f64 * 0.1);

    wave_write_16bit_mono_safer3("ex8_12.wav", &pcm3);
}

#[allow(non_snake_case)]
fn ex10_4() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 4; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let mut ac: Vec<c_double> = vec![0.0; pcm_length];
    let mut am: Vec<c_double> = vec![0.0; pcm_length];

    /* キャリア振幅 */
    let gate = pcm_fs * 4;
    let duration = pcm_fs * 4;
    let A = 0;
    let D = pcm_fs * 4;
    let S = 0.0;
    let R = pcm_fs * 4;
    safe_ADSR(&mut ac, A, D, S, R, gate, duration);

    let fc = 440.0; /* キャリア周波数 */

    /* モジュレータ振幅 */
    let gate = pcm_fs * 4;
    let duration = pcm_fs * 4;
    let A = 0;
    let D = pcm_fs * 2;
    let S = 0.0;
    let R = pcm_fs * 2;
    safe_ADSR(&mut am, A, D, S, R, gate, duration);

    let ratio = 3.5;
    let fm = fc * ratio; /* モジュレータ周波数 */

    /* FM音源 */
    for n in 0..pcm_length {
        pcm.s[n] = ac[n]
            * (2.0 * PI * fc * n as f64 / pcm_fs as f64
                + am[n] * (2.0 * PI * fm * n as f64 / pcm_fs as f64).sin())
                .sin();
    }

    pcm.mult_constant_gain(0.1);

    wave_write_16bit_mono_safer3("ex10_4.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex11_7() {
    let pcm0 = wave_read_16bit_mono_safer3("vocal.wav");
    wave_write_PCMU_mono_safer3("pcmu.wav", &pcm0);
    let pcm1 = wave_read_PCMU_mono_safer3("pcmu.wav");
    wave_write_16bit_mono_safer3("ex11_7_pcm.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex11_8() {
    let pcm0 = wave_read_16bit_mono_safer3("vocal.wav");
    wave_write_PCMA_mono_safer3("pcma.wav", &pcm0);
    let pcm1 = wave_read_PCMA_mono_safer3("pcma.wav");
    wave_write_16bit_mono_safer3("ex11_8_pcm.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex11_9() {
    let pcm0 = wave_read_16bit_mono_safer3("vocal.wav");
    wave_write_IMA_ADPCM_mono_safer3("ima_adpcm.wav", &pcm0);
    let pcm1 = wave_read_IMA_ADPCM_mono_safer3("ima_adpcm.wav");
    wave_write_16bit_mono_safer3("ex11_9_pcm.wav", &pcm1);
}
