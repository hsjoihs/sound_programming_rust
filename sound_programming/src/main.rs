extern crate num_complex;
extern crate rand;
extern crate sound_programming;
extern crate wave_utils;
//use std::io::Write;
use sound_programming::first::first;
use sound_programming::determine_J;
use sound_programming::mult;
use num_complex::Complex;
use rand::Rng;
use wave_utils::MonoPcm;
use wave_utils::c_double;
use wave_utils::fft::safe_FFT_;
use wave_utils::fft::safe_IFFT_;
use wave_utils::filter::safe_FIR_LPF;
use wave_utils::filter::safe_IIR_LPF;
use wave_utils::filter::safe_IIR_filtering;
use wave_utils::filter::safe_IIR_resonator;
use wave_utils::safe_ADSR;
use wave_utils::safe_Hanning_window;
use wave_utils::sinc;
use wave_utils::wave::wave_read_16bit_mono_safer3;
use wave_utils::wave_read_IMA_ADPCM_mono_safer3;
use wave_utils::wave_read_PCMA_mono_safer3;
use wave_utils::wave_read_PCMU_mono_safer3;
use wave_utils::wave_write_16bit_mono_safer3;
use wave_utils::wave_write_IMA_ADPCM_mono_safer3;
use wave_utils::wave_write_PCMA_mono_safer3;
use wave_utils::wave_write_PCMU_mono_safer3;
use std::f64::consts::PI;
//use std::io;


fn second() {
    ex7_1();
    ex7_2();
    ex7_3();
    ex7_4();
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
}

fn third(){
    ex9_1();
    ex9_2();
    ex9_3();
    ex9_4();
    ex9_5();
    ex9_6();
    ex9_7();
    ex9_8();
    ex9_9();
    ex9_10();
    ex9_11();
    ex10_1();
    ex10_2();
    ex10_3();
    ex10_4();
    ex10_5();
    ex10_6();
    ex11_7();
    ex11_8();
    ex11_9();
}

fn main() {
    if true {
        first();
    }
    if false {
        second();
    }
    if false{
        third();
    }
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

    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    let mut a = [0.0; 3];
    let mut b = [0.0; 3];

    filter_with_IIR_LPF(&pcm0, &mut pcm1, Q, I, J, &fc, &mut a, &mut b);

    wave_write_16bit_mono_safer3("ex7_1.wav", &pcm1);
}

#[allow(non_snake_case)]
fn filter_with_IIR_LPF(
    pcm0: &MonoPcm,
    pcm1: &mut MonoPcm,
    Q: f64,
    I: usize,
    J: usize,
    fc: &[f64],
    mut a: &mut [f64],
    mut b: &mut [f64],
) {
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

    filter_with_IIR_LPF(&pcm0, &mut pcm1, Q, I, J, &fc, &mut a, &mut b);

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

    let mut w = vec![0.0; N];
    safe_Hanning_window(&mut w); /* ハニング窓 */
    let number_of_frame = (pcm0.length - N / 2) / (N / 2); /* フレームの数 */
    let band_width = 8;
    let number_of_band = N / 2 / band_width;

    for frame in 0..number_of_frame {
        let offset = N / 2 * frame;
        /* X(n) */
        let mut x: Vec<_> = (0..N)
            .map(|n| Complex::new(pcm0.s[offset + n] * w[n], 0.0))
            .collect();
        safe_FFT_(&mut x);

        /* B(k) */
        let mut b_: Vec<_> = (0..N)
            .map(|n| Complex::new(pcm1.s[offset + n] * w[n], 0.0))
            .collect();
        safe_FFT_(&mut b_);

        for item in b_.iter_mut() {
            *item = Complex::new(item.norm_sqr().sqrt(), 0.0);
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
        let mut y: Vec<_> = (0..N).map(|k| x[k] * b_[k]).collect();
        safe_IFFT_(&mut y);

        let offset = N / 2 * frame;
        /* オーバーラップアド */
        for (n,item) in y.iter().enumerate() {
            pcm2.s[offset + n] += item.re;
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

// random; omitted from the test
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

// random; omitted from the test
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

fn square_wave(pcm: &mut MonoPcm, f0: c_double, gain: c_double, offset: usize, duration: usize) {
    let mut s = vec![0.0; duration];
    /* 矩形波 */
    let t0 = (pcm.fs as f64 / f0) as usize; /* 基本周期 */
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
        pcm.s[offset + n] += s[n];
    }
}

fn triangle_wave(pcm: &mut MonoPcm, f0: c_double, gain: c_double, offset: usize, duration: usize) {
    let mut s = vec![0.0; duration];
    /* 三角波 */
    let t0 = (pcm.fs as f64 / f0) as usize; /* 基本周期 */
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
        pcm.s[offset + n] += s[n];
    }
}

fn white_noise(pcm: &mut MonoPcm, gain: c_double, offset: usize, duration: usize) {
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
        pcm.s[offset + n] += s[n];
    }
}

// random; omitted from the test
fn ex8_6() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = 7000 * 16; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length); /* 音データ */

    /* メロディパート */
    square_wave(&mut pcm, 659.26, 0.1, 7000 * 0, 6125); /* E5 */
    square_wave(&mut pcm, 659.26, 0.1, 7000 * 1, 6125); /* E5 */
    square_wave(&mut pcm, 659.26, 0.1, 7000 * 3, 6125); /* E5 */
    square_wave(&mut pcm, 523.25, 0.1, 7000 * 5, 6125); /* C5 */
    square_wave(&mut pcm, 659.26, 0.1, 7000 * 6, 6125); /* E5 */
    square_wave(&mut pcm, 783.99, 0.1, 7000 * 8, 6125); /* G5 */
    square_wave(&mut pcm, 392.00, 0.1, 7000 * 12, 6125); /* G4 */

    /* ベースパート */
    triangle_wave(&mut pcm, 146.83, 0.2, 7000 * 0, 6125); /* D3 */
    triangle_wave(&mut pcm, 146.83, 0.2, 7000 * 1, 6125); /* D3 */
    triangle_wave(&mut pcm, 146.83, 0.2, 7000 * 3, 6125); /* D3 */
    triangle_wave(&mut pcm, 146.83, 0.2, 7000 * 5, 6125); /* D3 */
    triangle_wave(&mut pcm, 146.83, 0.2, 7000 * 6, 6125); /* D3 */
    triangle_wave(&mut pcm, 196.00, 0.2, 7000 * 8, 6125); /* G3 */
    triangle_wave(&mut pcm, 196.00, 0.2, 7000 * 12, 6125); /* G3 */

    /* パーカッション */
    white_noise(&mut pcm, 0.1, 7000 * 0, 4000);
    white_noise(&mut pcm, 0.1, 7000 * 2, 1000);
    white_noise(&mut pcm, 0.1, 7000 * 3, 4000);
    white_noise(&mut pcm, 0.1, 7000 * 5, 1000);
    white_noise(&mut pcm, 0.1, 7000 * 6, 4000);
    white_noise(&mut pcm, 0.1, 7000 * 8, 4000);
    white_noise(&mut pcm, 0.1, 7000 * 11, 4000);
    white_noise(&mut pcm, 0.1, 7000 * 13, 1000);
    white_noise(&mut pcm, 0.1, 7000 * 14, 1000);
    white_noise(&mut pcm, 0.1, 7000 * 15, 1000);

    wave_write_16bit_mono_safer3("ex8_6.wav", &pcm);
}

fn ex8_7() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = (pcm_fs as f64 * 0.6) as usize - 1; /* 音データの長さ */
    // -1 is here to match the file with the original data
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let mut f0 = vec![0.0; pcm_length];
    /* 基本周波数 */
    for n in 0..mult(pcm_fs, 0.1) as usize {
        f0[n] = 987.77; /* B5 */
    }
    for n in mult(pcm_fs, 0.1) as usize..pcm_length {
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
    let J = determine_J(delta); /* 遅延器の数 */
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

#[allow(non_snake_case)]
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
    let J = determine_J(delta); /* 遅延器の数 */
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

fn ex9_1() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 2; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);
    let vco = 500.0; /* 基本周波数 */

    /* ノコギリ波 */
    let t0 = (pcm.fs as f64 / vco) as usize; /* 基本周期 */
    {
        let mut m = 0;
        for n in 0..pcm_length {
            pcm.s[n] = 1.0 - 2.0 * m as f64 / t0 as f64;

            m += 1;
            if m >= t0 {
                m = 0;
            }
        }
    }

    let mut vca = vec![0.0; pcm_length];

    /* 時間エンベロープ */
    vca[0] = 1.0;
    let am = 0.2; /* LFOの振幅 */
    let fm = 2.0; /* LFOの周波数 */
    /* LFO */
    for n in 0..pcm_length {
        vca[n] = vca[0] + am * (2.0 * PI * fm * n as f64 / pcm.fs as f64).sin();
    }

    let gain = 0.1; /* ゲイン */

    for n in 0..pcm.length {
        pcm.s[n] *= vca[n] * gain;
    }

    wave_write_16bit_mono_safer3("ex9_1.wav", &pcm);
}

fn ex9_2() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 2; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);
    let mut vco = vec![0.0; pcm_length];
    /* 時間エンベロープ */
    vco[0] = 500.0; /* Hz */
    let am = 100.0; /* LFOの振幅 */
    let fm = 2.0; /* LFOの周波数 */
    /* LFO */
    for n in 0..pcm.length {
        vco[n] = vco[0] + am * (2.0 * PI * fm * n as f64 / pcm.fs as f64).sin();
    }

    /* ノコギリ波 */
    let mut t0 = (pcm.fs as f64 / vco[0]) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..pcm.length {
        pcm.s[n] = 1.0 - 2.0 * m as f64 / t0 as f64;

        m += 1;
        if m >= t0 {
            t0 = (pcm.fs as f64 / vco[n]) as usize; /* 基本周期 */
            m = 0;
        }
    }

    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex9_2.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex9_3() {
    let pcm0_fs = 44100;
    let pcm0_length = 44100 * 2;
    let mut pcm0 = MonoPcm::new16(pcm0_fs, pcm0_length);

    let vco = 500.0; /* 基本周波数 */

    /* ノコギリ波 */
    let t0 = (pcm0.fs as f64 / vco) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..pcm0.length {
        pcm0.s[n] = 1.0 - 2.0 * m as f64 / t0 as f64;

        m += 1;
        if m >= t0 {
            m = 0;
        }
    }

    let mut vcf = vec![0.0; pcm0.length];

    /* 時間エンベロープ */
    vcf[0] = 1000.0; /* Hz */
    let am = 800.0; /* LFOの振幅 */
    let fm = 2.0; /* LFOの周波数 */
    /* LFO */
    for n in 0..pcm0.length {
        vcf[n] = vcf[0] + am * (2.0 * PI * fm * n as f64 / pcm0.fs as f64).sin();
    }

    let Q = 5.0; /* レゾナンス */
    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */

    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    /* フィルタリング */
    for n in 0..pcm1.length {
        safe_IIR_LPF(vcf[n] / pcm1.fs as f64, Q, &mut a, &mut b); /* IIRフィルタの設計 */

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

    pcm1.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex9_3.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex9_4() {
    let pcm0_fs = 44100; /* 標本化周波数 */
    let pcm0_length = pcm0_fs * 4; /* 音データの長さ */
    let mut pcm0 = MonoPcm::new16(pcm0_fs, pcm0_length);

    let vco = 440.0; /* 基本周波数 */

    /* ノコギリ波 */
    let t0 = (pcm0.fs as f64 / vco) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..pcm0.length {
        pcm0.s[n] = 1.0 - 2.0 * m as f64 / t0 as f64;

        m += 1;
        if m >= t0 {
            m = 0;
        }
    }

    let vcf = 1500.0; /* 遮断周波数 */
    let Q = 5.0; /* レゾナンス */
    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */
    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    safe_IIR_LPF(vcf / pcm0.fs as f64, Q, &mut a, &mut b); /* IIRフィルタの設計 */

    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    /* フィルタリング */
    for n in 0..pcm1.length {
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

    let mut vca = vec![0.0; pcm0.length]; /* 振幅 */
    let gate = pcm1.fs * 4;
    let duration = pcm1.fs * 4;
    let A = 0;
    let D = 0;
    let S = 1.0;
    let R = 0;
    safe_ADSR(&mut vca, A, D, S, R, gate, duration);
    let gain = 0.1; /* ゲイン */
    for n in 0..pcm1.length {
        pcm1.s[n] *= vca[n] * gain;
    }
    /* フェード処理 */
    for n in 0..(pcm1.fs as f64 * 0.01).ceil() as usize {
        pcm1.s[n] *= n as f64 / (pcm1.fs as f64 * 0.01);
        pcm1.s[pcm1.length - n - 1] *= n as f64 / (pcm1.fs as f64 * 0.01);
    }
    wave_write_16bit_mono_safer3("ex9_4.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex9_5() {
    let pcm0_fs = 44100; /* 標本化周波数 */
    let pcm0_length = pcm0_fs * 4; /* 音データの長さ */
    let mut pcm0 = MonoPcm::new16(pcm0_fs, pcm0_length);

    let vco = 440.0; /* 基本周波数 */

    /* ノコギリ波 */
    let t0 = (pcm0.fs as f64 / vco) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..pcm0.length {
        pcm0.s[n] = 1.0 - 2.0 * m as f64 / t0 as f64;

        m += 1;
        if m >= t0 {
            m = 0;
        }
    }

    let vcf = 4000.0; /* 遮断周波数 */
    let Q = 1.0 / 2.0f64.sqrt(); /* レゾナンス */
    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */
    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    safe_IIR_LPF(vcf / pcm0.fs as f64, Q, &mut a, &mut b); /* IIRフィルタの設計 */

    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    /* フィルタリング */
    for n in 0..pcm1.length {
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
    let mut vca = vec![0.0; pcm0.length]; /* 振幅 */
    let gate = pcm1.fs * 3;
    let duration = pcm1.fs * 4;
    let A = pcm1.fs * 1;
    let D = 0;
    let S = 1.0;
    let R = pcm1.fs * 1;

    safe_ADSR(&mut vca, A, D, S, R, gate, duration);

    let gain = 0.1; /* ゲイン */
    for n in 0..pcm1.length {
        pcm1.s[n] *= vca[n] * gain;
    }

    /* フェード処理 */
    for n in 0..(pcm1.fs as f64 * 0.01).ceil() as usize {
        pcm1.s[n] *= n as f64 / (pcm1.fs as f64 * 0.01);
        pcm1.s[pcm1.length - n - 1] *= n as f64 / (pcm1.fs as f64 * 0.01);
    }
    wave_write_16bit_mono_safer3("ex9_5.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex9_6() {
    let pcm0_fs = 44100; /* 標本化周波数 */
    let pcm0_length = pcm0_fs * 4; /* 音データの長さ */

    let mut pcm0 = MonoPcm::new16(pcm0_fs, pcm0_length);

    let vco = 440.0; /* 基本周波数 */

    /* ノコギリ波 */
    let t0 = (pcm0.fs as f64 / vco) as usize; /* 基本周期 */
    let mut m = 0;
    for n in 0..pcm0.length {
        pcm0.s[n] = 1.0 - 2.0 * m as f64 / t0 as f64;

        m += 1;
        if m >= t0 {
            m = 0;
        }
    }

    let mut vcf = vec![0.0; pcm0.length]; /* 遮断周波数 */
    let gate = pcm0.fs * 1;
    let duration = pcm0.fs * 4;
    let A = 0;
    let D = pcm0.fs * 1;
    let S = 0.0;
    let R = pcm0.fs * 1;
    safe_ADSR(&mut vcf, A, D, S, R, gate, duration);
    let offset = 500.0; /* 時間エンベロープのオフセット */
    let depth = 500.0; /* 時間エンベロープのデプス */

    for n in 0..pcm0.length {
        vcf[n] = offset + vcf[n] * depth;
    }
    let Q = 1.0 / 2.0f64.sqrt(); /* レゾナンス */
    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */
    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    let mut pcm1 = MonoPcm::blank_copy(&pcm0);
    /* フィルタリング */
    for n in 0..pcm1.length {
        safe_IIR_LPF(vcf[n] / pcm1.fs as f64, Q, &mut a, &mut b); /* IIRフィルタの設計 */

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
    let mut vca = vec![0.0; pcm1.length]; /* 振幅 */

    let gate = pcm0.fs * 4;
    let duration = pcm0.fs * 4;
    let A = 0;
    let D = pcm0.fs * 4;
    let S = 0.0;
    let R = pcm0.fs * 1;
    safe_ADSR(&mut vca, A, D, S, R, gate, duration);

    let gain = 0.1; /* ゲイン */

    for n in 0..pcm1.length {
        pcm1.s[n] *= vca[n] * gain;
    }
    /* フェード処理 */
    for n in 0..(pcm1.fs as f64 * 0.01).ceil() as usize {
        pcm1.s[n] *= n as f64 / (pcm1.fs as f64 * 0.01);
        pcm1.s[pcm1.length - n - 1] *= n as f64 / (pcm1.fs as f64 * 0.01);
    }
    wave_write_16bit_mono_safer3("ex9_6.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex9_7() {
    let pcm0_fs = 44100; /* 標本化周波数 */
    let pcm0_length = pcm0_fs * 1; /* 音データの長さ */
    let mut pcm0 = MonoPcm::new16(pcm0_fs, pcm0_length);

    let mut vco = vec![0.0; pcm0.length]; /* 基本周波数 */
    let gate = pcm0.fs * 1;
    let duration = pcm0.fs * 1;
    let A = 0;
    let D = (pcm0.fs as f64 * 0.4) as usize;
    let S = 0.0;
    let R = (pcm0.fs as f64 * 0.4) as usize;
    safe_ADSR(&mut vco, A, D, S, R, gate, duration);
    let offset = 40.0; /* 時間エンベロープのオフセット */
    let depth = 120.0; /* 時間エンベロープのデプス */
    for n in 0..pcm0.length {
        vco[n] = offset + vco[n] * depth;
    }
    {
        /* 三角波 */
        let mut t0 = (pcm0.fs as f64 / vco[0]) as usize; /* 基本周期 */
        let mut m = 0;
        for n in 0..pcm0.length {
            pcm0.s[n] = if (m as f64) < t0 as f64 / 2.0 {
                -1.0 + 4.0 * m as f64 / t0 as f64
            } else {
                3.0 - 4.0 * m as f64 / t0 as f64
            };

            m += 1;
            if m >= t0 {
                t0 = (pcm0.fs as f64 / vco[n]) as usize; /* 基本周期 */
                m = 0;
            }
        }
    }

    let mut vcf = vec![0.0; pcm0.length]; /* 遮断周波数 */
    let gate = pcm0.fs * 1;
    let duration = pcm0.fs * 1;
    let A = 0;
    let D = (pcm0.fs as f64 * 0.4) as usize;
    let S = 0.0;
    let R = (pcm0.fs as f64 * 0.4) as usize;
    safe_ADSR(&mut vcf, A, D, S, R, gate, duration);

    let offset = 80.0; /* 時間エンベロープのオフセット */
    let depth = 240.0; /* 時間エンベロープのデプス */
    for n in 0..pcm0.length {
        vcf[n] = offset + vcf[n] * depth;
    }
    let Q = 5.0; /* レゾナンス */
    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */

    let mut pcm1 = MonoPcm::blank_copy(&pcm0);
    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    /* フィルタリング */
    for n in 0..pcm1.length {
        safe_IIR_LPF(vcf[n] / pcm1.fs as f64, Q, &mut a, &mut b); /* IIRフィルタの設計 */

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

    let mut vca = vec![0.0; pcm1.length]; /* 振幅 */
    let gate = pcm1.fs * 1;
    let duration = pcm1.fs * 1;
    let A = 0;
    let D = (pcm0.fs as f64 * 0.4) as usize;
    let S = 0.0;
    let R = (pcm0.fs as f64 * 0.4) as usize;
    safe_ADSR(&mut vca, A, D, S, R, gate, duration);

    let gain = 0.9; /* ゲイン */
    for n in 0..pcm1.length {
        pcm1.s[n] *= vca[n] * gain;
    }

    /* フェード処理 */
    for n in 0..(pcm1.fs as f64 * 0.01).ceil() as usize {
        pcm1.s[n] *= n as f64 / (pcm1.fs as f64 * 0.01);
        pcm1.s[pcm1.length - n - 1] *= n as f64 / (pcm1.fs as f64 * 0.01);
    }
    wave_write_16bit_mono_safer3("ex9_7.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex9_8() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 1; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let vco = 500.0; /* 基本周波数 */

    /* パルス列 */
    let N = 128;
    generate_pulse_sequence(&mut pcm, vco, N);

    let mut s = vec![0.0; pcm.length];

    /* 積分フィルタ */
    s[0] = pcm.s[0] - 0.5;
    for n in 1..pcm.length {
        s[n] = pcm.s[n] + 0.98 * s[n - 1];
    }
    for n in 0..pcm.length {
        pcm.s[n] = s[n] * 2.0;
    }

    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex9_8.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex9_9() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 1; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let vco = 500.0; /* 基本周波数 */

    /* 双極性パルス列 */
    let t0 = pcm.fs as f64 / vco; /* 基本周期 */
    let mut t = 0.0;
    let mut sign = 1.0;
    let N = 128;
    while t < pcm.length as f64 {
        let ta = t as i32;

        let tb = if t == ta as f64 { ta } else { ta + 1 };

        for n in (tb - N / 2)..=(ta + N / 2) {
            if n >= 0 && n < pcm.length as i32 {
                pcm.s[n as usize] += sign * sinc(PI * (t - n as f64))
                    * (0.5 + 0.5 * (2.0 * PI * (t - n as f64) / (N * 2 + 1) as f64).cos());
            }
        }

        t += t0 / 2.0;
        sign *= -1.0;
    }
    let mut s = vec![0.0; pcm.length];

    /* 積分フィルタ */
    s[0] = pcm.s[0] - 0.5;
    for n in 1..pcm.length {
        s[n] = pcm.s[n] + 0.98 * s[n - 1];
    }

    for n in 0..pcm.length {
        pcm.s[n] = s[n] * 2.0;
    }

    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex9_9.wav", &pcm);
}

#[allow(non_snake_case, unused_mut, unused_variables)]
fn ex9_10() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 1; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);
    let vco = 500.0; /* 基本周波数 */
    /* 双極性パルス列 */
    let t0 = pcm.fs as f64 / vco; /* 基本周期 */
    let mut t = 0.0;
    let mut sign = 1.0;
    let N = 128;

    while t < pcm.length as f64 {
        let ta = t as i32;

        let tb = if t == ta as f64 { ta } else { ta + 1 };

        for n in (tb - N / 2)..=(ta + N / 2) {
            if n >= 0 && n < pcm.length as i32 {
                pcm.s[n as usize] += sign * sinc(PI * (t - n as f64))
                    * (0.5 + 0.5 * (2.0 * PI * (t - n as f64) / (N * 2 + 1) as f64).cos());
            }
        }

        t += t0 / 2.0;
        sign *= -1.0;
    }
    let mut s = vec![0.0; pcm.length];

    /* 積分フィルタ */
    s[0] = pcm.s[0] - 0.5;
    for n in 1..pcm.length {
        s[n] = pcm.s[n] + 0.98 * s[n - 1];
    }

    for n in 0..pcm.length {
        pcm.s[n] = s[n] * 2.0;
    }

    for n in 0..pcm.length {
        pcm.s[n] *= 2.0 / t0;
    }

    /* 積分フィルタ */
    s[0] = pcm.s[0] - 0.5;
    for n in 1..pcm.length {
        s[n] = pcm.s[n] + 0.98 * s[n - 1];
    }

    for n in 0..pcm.length {
        pcm.s[n] = s[n] * 2.0;
    }
    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex9_10.wav", &pcm);
}

#[allow(non_snake_case)]
fn generate_pulse_sequence(pcm: &mut MonoPcm, vco: f64, N: usize) {
    let t0 = pcm.fs as f64 / vco; /* 基本周期 */
    let mut t = 0.0;
    while t < pcm.length as f64 {
        let ta = t as usize;

        let tb = if t == ta as f64 { ta } else { ta + 1 };

        for n in (tb as i32 - N as i32 / 2)..=(ta as i32 + N as i32 / 2) {
            if n >= 0 && n < pcm.length as i32 {
                pcm.s[n as usize] += sinc(PI * (t - n as f64))
                    * (0.5 + 0.5 * (2.0 * PI * (t - n as f64) / (N * 2 + 1) as f64).cos());
            }
        }

        t += t0;
    }

    for n in 0..pcm.length {
        pcm.s[n] -= 1.0 / t0 as f64;
    }
}

#[allow(non_snake_case)]
fn ex9_11() {
    let pcm0_fs = 44100; /* 標本化周波数 */
    let pcm0_length = pcm0_fs * 4; /* 音データの長さ */
    let mut pcm0 = MonoPcm::new16(pcm0_fs, pcm0_length);

    let mut s = vec![0.0; pcm0.length];
    let vco = 440.0; /* 基本周波数 */

    /* パルス列 */
    let N = 128;
    generate_pulse_sequence(&mut pcm0, vco, N);

    /* 積分フィルタ */
    s[0] = pcm0.s[0] - 0.5;
    for n in 1..pcm0.length {
        s[n] = pcm0.s[n] + 0.98 * s[n - 1];
    }

    for n in 0..pcm0.length {
        pcm0.s[n] = s[n] * 2.0;
    }
    let vcf = 4000.0; /* 遮断周波数 */
    let Q = 1.0 / 2.0f64.sqrt(); /* レゾナンス */
    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */
    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    safe_IIR_LPF(vcf / pcm0.fs as f64, Q, &mut a, &mut b); /* IIRフィルタの設計 */

    let mut pcm1 = MonoPcm::blank_copy(&pcm0);
    /* フィルタリング */
    for n in 0..pcm1.length {
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

    let mut vca = vec![0.0; pcm1.length]; /* 振幅 */
    let gate = pcm1.fs * 3;
    let duration = pcm1.fs * 4;
    let A = pcm1.fs * 1;
    let D = 0;
    let S = 1.0;
    let R = pcm1.fs * 1;
    safe_ADSR(&mut vca, A, D, S, R, gate, duration);
    let gain = 0.1; /* ゲイン */

    for n in 0..pcm1.length {
        pcm1.s[n] *= vca[n] * gain;
    }

    let pcm2_fs = 44100; /* 標本化周波数 */
    let pcm2_length = pcm2_fs * 4; /* 音データの長さ */
    let mut pcm2 = MonoPcm::new16(pcm2_fs, pcm2_length);
    let vco = 440.5; /* 基本周波数 */

    /* パルス列 */
    generate_pulse_sequence(&mut pcm2, vco, N);

    /* 積分フィルタ */
    s[0] = pcm2.s[0] - 0.5;
    for n in 1..pcm2.length {
        s[n] = pcm2.s[n] + 0.98 * s[n - 1];
    }

    for n in 0..pcm2.length {
        pcm2.s[n] = s[n] * 2.0;
    }

    let vcf = 4000.0; /* 遮断周波数 */
    let Q = 1.0 / 2.0f64.sqrt(); /* レゾナンス */
    let I = 2; /* 遅延器の数 */
    let J = 2; /* 遅延器の数 */
    safe_IIR_LPF(vcf / pcm2.fs as f64, Q, &mut a, &mut b); /* IIRフィルタの設計 */

    let mut pcm3 = MonoPcm::blank_copy(&pcm2);
    /* フィルタリング */
    for n in 0..pcm3.length {
        for m in 0..=J {
            if n >= m {
                pcm3.s[n] += b[m] * pcm2.s[n - m];
            }
        }
        for m in 1..=I {
            if n >= m {
                pcm3.s[n] += -a[m] * pcm3.s[n - m];
            }
        }
    }
    let mut vca = vec![0.0; pcm3.length]; /* 振幅 */
    let gate = pcm3.fs * 3;
    let duration = pcm3.fs * 4;
    let A = pcm3.fs * 1;
    let D = 0;
    let S = 1.0;
    let R = pcm3.fs * 1;
    safe_ADSR(&mut vca, A, D, S, R, gate, duration);
    let gain = 0.1; /* ゲイン */

    for n in 0..pcm3.length {
        pcm3.s[n] *= vca[n] * gain;
    }

    /* デチューン */
    for n in 0..pcm3.length {
        pcm1.s[n] += pcm3.s[n];
    }
    wave_write_16bit_mono_safer3("ex9_11.wav", &pcm1);
}

#[allow(non_snake_case, unused_mut, unused_variables)]
fn ex10_1() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 1; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);
    let ac = 1.0; /* キャリア振幅 */
    let fc = 500.0; /* キャリア周波数 */

    let am = 1.0; /* モジュレータ振幅 */
    let ratio = 1.0; /* 周波数比 */
    let fm = fc * ratio; /* モジュレータ周波数 */

    /* FM音源 */
    for n in 0..pcm.length {
        pcm.s[n] = ac * (2.0 * PI * fc * n as f64 / pcm.fs as f64
            + am * (2.0 * PI * fm * n as f64 / pcm.fs as f64).sin())
            .sin();
    }

    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex10_1.wav", &pcm);
}

fn ex10_2() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 1; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);
    let ac = 1.0; /* キャリア振幅 */
    let fc = 500.0; /* キャリア周波数 */

    let am = 1.0; /* モジュレータ振幅 */
    let ratio = 2.0; /* 周波数比 */
    let fm = fc * ratio; /* モジュレータ周波数 */
    /* FM音源 */
    for n in 0..pcm.length {
        pcm.s[n] = ac * (2.0 * PI * fc * n as f64 / pcm.fs as f64
            + am * (2.0 * PI * fm * n as f64 / pcm.fs as f64).sin())
            .sin();
    }
    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex10_2.wav", &pcm);
}

#[allow(non_snake_case)]
fn ex10_3() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 1; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let ac = 1.0; /* キャリア振幅 */
    let fc = 500.0; /* キャリア周波数 */

    let am = 1.0; /* モジュレータ振幅 */
    let ratio = 3.5; /* 周波数比 */
    let fm = fc * ratio; /* モジュレータ周波数 */
    /* FM音源 */
    for n in 0..pcm.length {
        pcm.s[n] = ac * (2.0 * PI * fc * n as f64 / pcm.fs as f64
            + am * (2.0 * PI * fm * n as f64 / pcm.fs as f64).sin())
            .sin();
    }
    pcm.mult_constant_gain(0.1);
    wave_write_16bit_mono_safer3("ex10_3.wav", &pcm);
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
fn ex10_5() {
    let pcm0_fs = 44100; /* 標本化周波数 */
    let pcm0_length = pcm0_fs * 4; /* 音データの長さ */
    let mut pcm0 = MonoPcm::new16(pcm0_fs, pcm0_length);

    let mut pcm1 = MonoPcm::blank_copy(&pcm0);

    let mut ac = vec![0.0; pcm0.length];
    let mut am = vec![0.0; pcm0.length];
    {
        /* キャリア振幅 */
        let gate = pcm0.fs * 4;
        let duration = pcm0.fs * 4;
        let A = 0;
        let D = pcm0.fs * 4;
        let S = 0.0;
        let R = pcm0.fs * 1;
        safe_ADSR(&mut ac, A, D, S, R, gate, duration);
    }
    let fc = 440.0; /* キャリア周波数 */
    {
        /* モジュレータ振幅 */
        let gate = pcm0.fs * 4;
        let duration = pcm0.fs * 4;
        let A = 0;
        let D = pcm0.fs * 2;
        let S = 0.0;
        let R = pcm0.fs * 2;
        safe_ADSR(&mut am, A, D, S, R, gate, duration);
    }

    let ratio = 1.0;
    let fm = 440.0 * ratio; /* モジュレータ周波数 */
    /* FM音源 */
    for n in 0..pcm0.length {
        pcm0.s[n] = ac[n]
            * (2.0 * PI * fc * n as f64 / pcm0.fs as f64
                + am[n] * (2.0 * PI * fm * n as f64 / pcm0.fs as f64).sin())
                .sin();
    }
    {
        /* キャリア振幅 */
        let gate = pcm1.fs * 4;
        let duration = pcm1.fs * 4;
        let A = 0;
        let D = pcm1.fs * 1;
        let S = 0.0;
        let R = pcm1.fs * 1;
        safe_ADSR(&mut ac, A, D, S, R, gate, duration);
    }
    let fc = 440.0; /* キャリア周波数 */
    {
        /* モジュレータ振幅 */
        let gate = pcm1.fs * 4;
        let duration = pcm1.fs * 4;
        let A = 0;
        let D = pcm1.fs * 1;
        let S = 0.0;
        let R = pcm1.fs * 1;
        safe_ADSR(&mut am, A, D, S, R, gate, duration);
    }
    let ratio = 14.0;
    let fm = fc * ratio; /* モジュレータ周波数 */
    /* FM音源 */
    for n in 0..pcm1.length {
        pcm1.s[n] = ac[n]
            * (2.0 * PI * fc * n as f64 / pcm1.fs as f64
                + am[n] * (2.0 * PI * fm * n as f64 / pcm1.fs as f64).sin())
                .sin();
    }

    let gain = 0.1; /* ゲイン */

    for n in 0..pcm1.length {
        pcm1.s[n] += pcm0.s[n];
        pcm1.s[n] *= gain;
    }

    wave_write_16bit_mono_safer3("ex10_5.wav", &pcm1);
}

#[allow(non_snake_case, unused_mut, unused_variables)]
fn ex10_6() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = pcm_fs * 4; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length);

    let mut ac = vec![0.0; pcm.length];
    let mut am = vec![0.0; pcm.length];

    {
        /* キャリア振幅 */
        let gate = pcm.fs * 4;
        let duration = pcm.fs * 4;
        let A = 0;
        let D = pcm.fs * 4;
        let S = 0.0;
        let R = pcm.fs * 4;
        safe_ADSR(&mut ac, A, D, S, R, gate, duration);
    }
    let fc = 440.0; /* キャリア周波数 */
    {
        /* モジュレータ振幅 */
        let gate = pcm.fs * 4;
        let duration = pcm.fs * 4;
        let A = 0;
        let D = pcm.fs * 2;
        let S = 0.0;
        let R = pcm.fs * 2;
        safe_ADSR(&mut am, A, D, S, R, gate, duration);
    }
    let ratio = 3.5;
    let fm = fc * ratio; /* モジュレータ周波数 */
    /* AM変調 */
    for n in 0..pcm.length {
        pcm.s[n] = ac[n] * (2.0 * PI * fc * n as f64 / pcm.fs as f64).sin()
            * (1.0 + am[n] * (2.0 * PI * fm * n as f64 / pcm.fs as f64).sin());
    }

    pcm.mult_constant_gain(0.1);

    wave_write_16bit_mono_safer3("ex10_6.wav", &pcm);
}
/*



*/

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
