extern crate num_complex;
extern crate rand;

use rand::Rng;
use wave_utils::MonoPcm;
use wave_utils::wave::wave_write_16bit_mono_safer3;
extern crate fraction;
type F = fraction::GenericFraction<i32>;

#[derive(Copy, Clone, Debug)]
enum NoteName {
    C = 3, 
    CSh = 4, 
    D = 5, 
    DSh = 6,
    E = 7,
    F = 8,
    FSh = 9,
    G = 10,
    GSh = 11,
    A = 12,
    ASh = 13,
    B = 14
}

#[derive(Copy, Clone, Debug)]
struct Pitch(NoteName, i32);
impl Pitch {
    pub fn get_pitch(self) -> f64 {
        let Pitch(name, octave) = self;
        440. * 2f64.powf(1. / 12.).powi(name as i32 + 12 * (octave - 5))
    }
}

pub fn second() {
    let eighth_note_length = 8000; // number of samples per an eighth note
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = eighth_note_length * 16; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length); /* 音データ */

    let mut counter = 0;
    use self::NoteName::*;
    for pitch in vec![
        Some(Pitch(E, 5)),
        Some(Pitch(E, 5)),
        None,
        Some(Pitch(E, 5)),
        None,
        Some(Pitch(C, 5)),
        Some(Pitch(E, 5)),
        None,
        Some(Pitch(G, 5)),
        None,
        None,
        None,
        Some(Pitch(G, 4))
    ] {
        let freq = match pitch {
            Some(a) => Some(a.get_pitch()),
            None => None,
        };
        match freq {
            None => {},
            Some(f) => {
                square_wave(&mut pcm, f, 0.1, eighth_note_length * counter, 7000);
            }
        }

        counter += 1;
    }

    /* ベースパート */
    triangle_wave(&mut pcm, 146.83, 0.2, eighth_note_length * 0, 7000); /* D3 */
    triangle_wave(&mut pcm, 146.83, 0.2, eighth_note_length * 1, 7000); /* D3 */
    triangle_wave(&mut pcm, 146.83, 0.2, eighth_note_length * 3, 7000); /* D3 */
    triangle_wave(&mut pcm, 146.83, 0.2, eighth_note_length * 5, 7000); /* D3 */
    triangle_wave(&mut pcm, 146.83, 0.2, eighth_note_length * 6, 7000); /* D3 */
    triangle_wave(&mut pcm, 196.00, 0.2, eighth_note_length * 8, 7000); /* G3 */
    triangle_wave(&mut pcm, 196.00, 0.2, eighth_note_length * 12, 7000); /* G3 */

    /* パーカッション */
    white_noise(&mut pcm, 0.1, eighth_note_length * 0, 4000);
    white_noise(&mut pcm, 0.1, eighth_note_length * 2, 1000);
    white_noise(&mut pcm, 0.1, eighth_note_length * 3, 4000);
    white_noise(&mut pcm, 0.1, eighth_note_length * 5, 1000);
    white_noise(&mut pcm, 0.1, eighth_note_length * 6, 4000);
    white_noise(&mut pcm, 0.1, eighth_note_length * 8, 4000);
    white_noise(&mut pcm, 0.1, eighth_note_length * 11, 4000);
    white_noise(&mut pcm, 0.1, eighth_note_length * 13, 1000);
    white_noise(&mut pcm, 0.1, eighth_note_length * 14, 1000);
    white_noise(&mut pcm, 0.1, eighth_note_length * 15, 1000);

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
