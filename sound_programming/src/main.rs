extern crate num_complex;
extern crate rand;
extern crate sound_programming;
extern crate wave_utils;
//use std::io::Write;
use sound_programming::first::first;
use sound_programming::mult;
use sound_programming::second::second;
use std::f64::consts::PI;
use wave_utils::MonoPcm;
use wave_utils::c_double;
use wave_utils::filter::safe_IIR_LPF;
use wave_utils::safe_ADSR;
use wave_utils::sinc;
use wave_utils::wave::wave_read_16bit_mono_safer3;
use wave_utils::wave::wave_read_16bit_stereo_safer3;
use wave_utils::wave::wave_read_8bit_mono_safer3;
use wave_utils::wave::wave_read_8bit_stereo_safer3;
use wave_utils::wave::wave_read_IMA_ADPCM_mono_safer3;
use wave_utils::wave::wave_read_PCMA_mono_safer3;
use wave_utils::wave::wave_read_PCMU_mono_safer3;
use wave_utils::wave::wave_write_16bit_mono_safer3;
use wave_utils::wave::wave_write_16bit_stereo_safer3;
use wave_utils::wave::wave_write_8bit_mono_safer3;
use wave_utils::wave::wave_write_8bit_stereo_safer3;
use wave_utils::wave::wave_write_IMA_ADPCM_mono_safer3;
use wave_utils::wave::wave_write_PCMA_mono_safer3;
use wave_utils::wave::wave_write_PCMU_mono_safer3;
//use std::io;

fn third() {
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
    ex11_1();
    ex11_2();
    ex11_3();
    ex11_4();
    ex11_5();
    ex11_6();
}

fn main() {
    if true {
        first();
    }
    if false {
        second();
    }
    if false {
        third();
    }
    ex11_7();
    ex11_8();
    ex11_9();
    eightbit();
}

fn eightbit() {
    {
        let pcm0 = wave_read_16bit_mono_safer3("ex1_1_a.wav"); /* 音データの入力 */
        let mut pcm1 = pcm0.clone(); /* 音データのコピー */
        pcm1.bits = 8;

        wave_write_8bit_mono_safer3("ex1_1_8bit.wav", &pcm1); /* 音データの出力 */
    }
    {
        let pcm1 = wave_read_8bit_mono_safer3("ex1_1_8bit.wav");
        let mut pcm2 = pcm1.clone();
        pcm2.bits = 16;
        wave_write_16bit_mono_safer3("ex1_1_c.wav", &pcm2); /* 音データの出力 */
    }
    {
        let pcm0 = wave_read_16bit_stereo_safer3("ex1_2_a.wav");
        let mut pcm1 = pcm0.clone(); /* 音データのコピー */
        pcm1.bits = 8;

        wave_write_8bit_stereo_safer3("ex1_2_8bit.wav", &pcm1);
    }
    {
        let pcm1 = wave_read_8bit_stereo_safer3("ex1_2_8bit.wav");
        let mut pcm2 = pcm1.clone();
        pcm2.bits = 16;
        wave_write_16bit_stereo_safer3("ex1_2_c.wav", &pcm2); /* 音データの出力 */
    }
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

#[allow(non_snake_case)]
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

#[allow(non_snake_case)]
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

#[allow(non_snake_case)]
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

#[allow(non_snake_case, unused_mut, unused_variables)]
fn ex11_1() {
    let pcm0 = wave_read_16bit_mono_safer3("sine_2s.wav");
    let rate = 2.0;
    assert!(1.0 < rate);

    let pcm1_fs = pcm0.fs; /* 標本化周波数 */
    let pcm1_length = (pcm0.length as f64 / rate) as usize + 1; /* 音データの長さ */
    let mut pcm1 = MonoPcm::new16(pcm1_fs, pcm1_length);

    let template_size = mult(pcm1.fs, 0.01); /* 相関関数のサイズ */
    let pmin = mult(pcm1.fs, 0.005); /* ピークの探索範囲の下限 */
    let pmax = mult(pcm1.fs, 0.02); /* ピークの探索範囲の上限 */

    let mut x = vec![0.0; template_size];
    let mut y = vec![0.0; template_size];
    let mut r = vec![0.0; pmax + 1];

    let mut offset0 = 0;
    let mut offset1 = 0;

    while offset0 + pmax * 2 < pcm0.length {
        for n in 0..template_size {
            x[n] = pcm0.s[offset0 + n]; /* 本来の音データ */
        }

        let mut rmax = 0.0;
        let mut p = pmin;
        for m in pmin..=pmax {
            for n in 0..template_size {
                y[n] = pcm0.s[offset0 + m + n]; /* mサンプルずらした音データ */
            }
            r[m] = 0.0;
            for n in 0..template_size {
                r[m] += x[n] * y[n]; /* 相関関数 */
            }
            if r[m] > rmax {
                rmax = r[m]; /* 相関関数のピーク */
                p = m; /* 波形の周期 */
            }
        }

        for n in 0..p {
            pcm1.s[offset1 + n] = pcm0.s[offset0 + n] * (p - n) as f64 / p as f64; /* 単調減少の重みづけ */
            pcm1.s[offset1 + n] += pcm0.s[offset0 + p + n] * n as f64 / p as f64; /* 単調増加の重みづけ */
        }

        let q = (p as f64 / (rate - 1.0) + 0.5) as usize;
        for n in p..q {
            if offset0 + p + n >= pcm0.length {
                break;
            }
            pcm1.s[offset1 + n] = pcm0.s[offset0 + p + n];
        }

        offset0 += p + q; /* offset0の更新 */
        offset1 += q; /* offset1の更新 */
    }
    wave_write_16bit_mono_safer3("ex11_1.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex11_2() {
    let pcm0 = wave_read_16bit_mono_safer3("sine_1s.wav");
    let rate = 0.5;
    assert!(0.5 <= rate && rate < 1.0);

    let pcm1_fs = pcm0.fs; /* 標本化周波数 */
    let pcm1_length = (pcm0.length as f64 / rate) as usize + 1; /* 音データの長さ */
    let mut pcm1 = MonoPcm::new16(pcm1_fs, pcm1_length);

    let template_size = mult(pcm1.fs, 0.01); /* 相関関数のサイズ */
    let pmin = mult(pcm1.fs, 0.005); /* ピークの探索範囲の下限 */
    let pmax = mult(pcm1.fs, 0.02); /* ピークの探索範囲の上限 */

    let mut x = vec![0.0; template_size];
    let mut y = vec![0.0; template_size];
    let mut r = vec![0.0; pmax + 1];

    let mut offset0 = 0;
    let mut offset1 = 0;

    while offset0 + pmax * 2 < pcm0.length {
        for n in 0..template_size {
            x[n] = pcm0.s[offset0 + n]; /* 本来の音データ */
        }

        let mut rmax = 0.0;
        let mut p = pmin;
        for m in pmin..=pmax {
            for n in 0..template_size {
                y[n] = pcm0.s[offset0 + m + n]; /* mサンプルずらした音データ */
            }
            r[m] = 0.0;
            for n in 0..template_size {
                r[m] += x[n] * y[n]; /* 相関関数 */
            }
            if r[m] > rmax {
                rmax = r[m]; /* 相関関数のピーク */
                p = m; /* 波形の周期 */
            }
        }

        for n in 0..p {
            pcm1.s[offset1 + n] = pcm0.s[offset0 + n];
        }
        for n in 0..p {
            /* 単調減少の重みづけ */
            pcm1.s[offset1 + p + n] = pcm0.s[offset0 + p + n] * (p - n) as f64 / p as f64;

            /* 単調増加の重みづけ */
            pcm1.s[offset1 + p + n] += pcm0.s[offset0 + n] * n as f64 / p as f64;
        }

        let q = (p as f64 * rate / (1.0 - rate) + 0.5) as usize;
        for n in p..q {
            if offset0 + n >= pcm0.length {
                break;
            }
            pcm1.s[offset1 + p + n] = pcm0.s[offset0 + n];
        }

        offset0 += q; /* offset0の更新 */
        offset1 += p + q; /* offset1の更新 */
    }
    wave_write_16bit_mono_safer3("ex11_2.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex11_3() {
    let pcm0 = wave_read_16bit_mono_safer3("ex11_sine_500hz.wav");
    let pitch = 2.0; /* 音の高さを2倍にする */

    let pcm1_fs = pcm0.fs; /* 標本化周波数 */
    let pcm1_length = (pcm0.length as f64 / pitch) as usize; /* 音データの長さ */
    let mut pcm1 = MonoPcm::new16(pcm1_fs, pcm1_length);

    let N = 128; /* ハニング窓のサイズ */

    Hanning_something(&mut pcm1, &pcm0, pitch, N);

    wave_write_16bit_mono_safer3("ex11_3.wav", &pcm1);
}

#[allow(non_snake_case)]
fn something(t: f64, n: i32, N: usize) -> f64 {
    sinc(PI * (t - n as f64)) * (0.5 + 0.5 * (2.0 * PI * (t - n as f64) / (N * 2 + 1) as f64).cos())
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
                pcm.s[n as usize] += something(t, n, N);
            }
        }

        t += t0;
    }

    for n in 0..pcm.length {
        pcm.s[n] -= 1.0 / t0 as f64;
    }
}

#[allow(non_snake_case, unused_mut, unused_variables)]
fn ex11_4() {
    let pcm0 = wave_read_16bit_mono_safer3("sine_1000hz.wav");
    let rate = 2.0;
    assert!(1.0 < rate);

    let pcm1_fs = pcm0.fs; /* 標本化周波数 */
    let pcm1_length = (pcm0.length as f64 / rate) as usize + 1; /* 音データの長さ */
    let mut pcm1 = MonoPcm::new16(pcm1_fs, pcm1_length);

    let template_size = mult(pcm1.fs, 0.01); /* 相関関数のサイズ */
    let pmin = mult(pcm1.fs, 0.005); /* ピークの探索範囲の下限 */
    let pmax = mult(pcm1.fs, 0.02); /* ピークの探索範囲の上限 */

    let mut x = vec![0.0; template_size];
    let mut y = vec![0.0; template_size];
    let mut r = vec![0.0; pmax + 1];

    let mut offset0 = 0;
    let mut offset1 = 0;

    while offset0 + pmax * 2 < pcm0.length {
        for n in 0..template_size {
            x[n] = pcm0.s[offset0 + n]; /* 本来の音データ */
        }

        let mut rmax = 0.0;
        let mut p = pmin;
        for m in pmin..=pmax {
            for n in 0..template_size {
                y[n] = pcm0.s[offset0 + m + n]; /* mサンプルずらした音データ */
            }
            r[m] = 0.0;
            for n in 0..template_size {
                r[m] += x[n] * y[n]; /* 相関関数 */
            }
            if r[m] > rmax {
                rmax = r[m]; /* 相関関数のピーク */
                p = m; /* 波形の周期 */
            }
        }

        for n in 0..p {
            pcm1.s[offset1 + n] = pcm0.s[offset0 + n] * (p - n) as f64 / p as f64; /* 単調減少の重み付け */
            pcm1.s[offset1 + n] += pcm0.s[offset0 + p + n] * n as f64 / p as f64; /* 単調増加の重み付け */
        }

        let q = (p as f64 / (rate - 1.0) + 0.5) as usize;
        for n in p..q {
            if offset0 + p + n >= pcm0.length {
                break;
            }
            pcm1.s[offset1 + n] = pcm0.s[offset0 + p + n];
        }

        offset0 += p + q; /* offset0の更新 */
        offset1 += q; /* offset1の更新 */
    }

    let pitch = 1.0 / rate;

    let mut pcm2 = MonoPcm::blank_copy(&pcm0);

    let N = 128; /* ハニング窓のサイズ */
    Hanning_something(&mut pcm2, &pcm1, pitch, N);

    wave_write_16bit_mono_safer3("ex11_4.wav", &pcm2);
}

#[allow(non_snake_case, unused_mut, unused_variables)]
fn ex11_5() {
    let pcm0 = wave_read_16bit_mono_safer3("ex11_sine_500hz.wav");
    let rate = 0.5;
    assert!(0.5 <= rate && rate < 1.0);

    let pcm1_fs = pcm0.fs; /* 標本化周波数 */
    let pcm1_length = (pcm0.length as f64 / rate) as usize + 1; /* 音データの長さ */
    let mut pcm1 = MonoPcm::new16(pcm1_fs, pcm1_length);

    let template_size = mult(pcm1.fs, 0.01); /* 相関関数のサイズ */
    let pmin = mult(pcm1.fs, 0.005); /* ピークの探索範囲の下限 */
    let pmax = mult(pcm1.fs, 0.02); /* ピークの探索範囲の上限 */

    let mut x = vec![0.0; template_size];
    let mut y = vec![0.0; template_size];
    let mut r = vec![0.0; pmax + 1];

    let mut offset0 = 0;
    let mut offset1 = 0;

    while offset0 + pmax * 2 < pcm0.length {
        for n in 0..template_size {
            x[n] = pcm0.s[offset0 + n]; /* 本来の音データ */
        }

        let mut rmax = 0.0;
        let mut p = pmin;
        for m in pmin..=pmax {
            for n in 0..template_size {
                y[n] = pcm0.s[offset0 + m + n]; /* mサンプルずらした音データ */
            }
            r[m] = 0.0;
            for n in 0..template_size {
                r[m] += x[n] * y[n]; /* 相関関数 */
            }
            if r[m] > rmax {
                rmax = r[m]; /* 相関関数のピーク */
                p = m; /* 波形の周期 */
            }
        }

        for n in 0..p {
            pcm1.s[offset1 + n] = pcm0.s[offset0 + n];
        }
        for n in 0..p {
            /* 単調減少の重み付け */
            pcm1.s[offset1 + p + n] = pcm0.s[offset0 + p + n] * (p - n) as f64 / p as f64;

            /* 単調増加の重み付け */
            pcm1.s[offset1 + p + n] += pcm0.s[offset0 + n] * n as f64 / p as f64;
        }

        let q = (p as f64 * rate / (1.0 - rate) + 0.5) as usize;
        for n in p..q {
            if offset0 + n >= pcm0.length {
                break;
            }
            pcm1.s[offset1 + p + n] = pcm0.s[offset0 + n];
        }

        offset0 += q; /* offset0の更新 */
        offset1 += p + q; /* offset1の更新 */
    }

    let pitch = 1.0 / rate;
    let mut pcm2 = MonoPcm::blank_copy(&pcm0);
    let N = 128; /* ハニング窓のサイズ */
    Hanning_something(&mut pcm2, &pcm1, pitch, N);
    wave_write_16bit_mono_safer3("ex11_5.wav", &pcm2);
}

#[allow(non_snake_case)]
fn Hanning_something(pcm1: &mut MonoPcm, pcm0: &MonoPcm, pitch: f64, N: usize) {
    for o in 0..pcm1.length {
        let mut t = pitch * o as f64;
        let mut tmp = 0.0;

        let ta = t as usize;

        let tb = if t == ta as f64 { ta } else { ta + 1 };
        for m in (tb as i32 - N as i32 / 2)..=(ta as i32 + N as i32 / 2) {
            if m >= 0 && m < pcm0.length as i32 {
                tmp += pcm0.s[m as usize] * something(t, m, N);
            }
        }

        pcm1.s[o as usize] += tmp;
    }
}

#[allow(non_snake_case, unused_mut, unused_variables)]
fn ex11_6() {
    let pcm0 = wave_read_16bit_mono_safer3("vocal.wav");
    let rate = 0.5;
    assert!(0.5 <= rate && rate < 1.0);

    let pcm1_fs = pcm0.fs; /* 標本化周波数 */
    let pcm1_length = (pcm0.length as f64 / rate) as usize + 1; /* 音データの長さ */
    let mut pcm1 = MonoPcm::new16(pcm1_fs, pcm1_length);

    let template_size = mult(pcm1.fs, 0.01); /* 相関関数のサイズ */
    let pmin = mult(pcm1.fs, 0.005); /* ピークの探索範囲の下限 */
    let pmax = mult(pcm1.fs, 0.02); /* ピークの探索範囲の上限 */

    let mut x = vec![0.0; template_size];
    let mut y = vec![0.0; template_size];
    let mut r = vec![0.0; pmax + 1];

    let mut offset0 = 0;
    let mut offset1 = 0;

    while offset0 + pmax * 2 < pcm0.length {
        for n in 0..template_size {
            x[n] = pcm0.s[offset0 + n]; /* 本来の音データ */
        }

        let mut rmax = 0.0;
        let mut p = pmin;
        for m in pmin..=pmax {
            for n in 0..template_size {
                y[n] = pcm0.s[offset0 + m + n]; /* mサンプルずらした音データ */
            }
            r[m] = 0.0;
            for n in 0..template_size {
                r[m] += x[n] * y[n]; /* 相関関数 */
            }
            if r[m] > rmax {
                rmax = r[m]; /* 相関関数のピーク */
                p = m; /* 波形の周期 */
            }
        }

        for n in 0..p {
            pcm1.s[offset1 + n] = pcm0.s[offset0 + n];
        }
        for n in 0..p {
            /* 単調減少の重み付け */
            pcm1.s[offset1 + p + n] = pcm0.s[offset0 + p + n] * (p - n) as f64 / p as f64;

            /* 単調増加の重み付け */
            pcm1.s[offset1 + p + n] += pcm0.s[offset0 + n] * n as f64 / p as f64;
        }

        let q = (p as f64 * rate / (1.0 - rate) + 0.5) as usize;
        for n in p..q {
            if offset0 + n >= pcm0.length {
                break;
            }
            pcm1.s[offset1 + p + n] = pcm0.s[offset0 + n];
        }

        offset0 += q; /* offset0の更新 */
        offset1 += p + q; /* offset1の更新 */
    }

    let pitch = 1.0 / rate;

    let mut pcm2 = MonoPcm::blank_copy(&pcm0);

    let N = 128;

    Hanning_something(&mut pcm2, &pcm1, pitch, N);

    wave_write_16bit_mono_safer3("ex11_6.wav", &pcm2);
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
    let mut pcm0 = wave_read_16bit_mono_safer3("vocal.wav");
    pcm0.bits = 8;
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
