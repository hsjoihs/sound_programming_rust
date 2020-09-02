extern crate num_complex;
extern crate rand;

use rand::Rng;
use wave_utils::MonoPcm;
use wave_utils::wave::wave_write_16bit_mono_safer3;

pub fn second() {
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = 7000 * 16; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length); /* 音データ */

    let mut counter = 0;
    for pitch in "e5e5__e5__c5e5__g5______g4".chars().collect::<Vec<_>>().chunks(2) {
        let freq = match pitch {
            ['e', '5'] => Some(440. * 2f64.powf(1. / 12.).powi(7)),
            ['c', '5'] => Some(440. * 2f64.powf(1. / 12.).powi(3)),
            ['g', '5'] => Some(440. * 2f64.powf(1. / 12.).powi(10)),
            ['g', '4'] => Some(440. * 2f64.powf(1. / 12.).powi(-2)),
            ['_', '_'] => None,
            _ => panic!("unrecognized pitch {:?}", pitch)
        };
        match freq {
            None => {},
            Some(f) => {
                square_wave(&mut pcm, f, 0.1, 7000 * counter, 6125);
            }
        }

        counter += 1;
    }

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

fn square_wave(pcm: &mut MonoPcm, f0: f64, gain: f64, offset: usize, duration: usize) {
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

fn triangle_wave(pcm: &mut MonoPcm, f0: f64, gain: f64, offset: usize, duration: usize) {
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

fn white_noise(pcm: &mut MonoPcm, gain: f64, offset: usize, duration: usize) {
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
