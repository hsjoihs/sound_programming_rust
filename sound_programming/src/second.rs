extern crate num_complex;
extern crate rand;

use rand::Rng;
use wave_utils::wave::wave_write_16bit_mono_safer3;
use wave_utils::MonoPcm;

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
    B = 14,
}

#[derive(Copy, Clone, Debug)]
struct Pitch(NoteName, i32);
impl Pitch {
    pub fn to_hertz(self) -> f64 {
        let Pitch(name, octave) = self;
        504.997 * 2f64.powf(1. / 12.).powi(name as i32 + 12 * (octave - 5))
    }
}

fn render_pitched<F>(
    mut pcm: &mut wave_utils::MonoPcm,
    unit_length: usize,
    font: F,
    melody: &[(usize, Option<Pitch>)],
) where
    F: Fn(&mut wave_utils::MonoPcm, usize, usize, f64, usize) -> (),
{
    let mut counter = 0;
    for (len, pitch) in melody {
        if let Some(p) = pitch {
            font(&mut pcm, unit_length, counter, p.to_hertz(), *len);
        }
        counter += len;
    }
}

fn render_pitchless<F, A>(
    mut pcm: &mut wave_utils::MonoPcm,
    unit_length: usize,
    font: F,
    melody: &[(usize, Option<A>)],
) where
    F: Fn(&mut wave_utils::MonoPcm, usize, usize, A) -> (),
    A: Clone,
{
    let mut counter = 0;
    for (len, addendum) in melody {
        if let Some(s) = addendum {
            font(&mut pcm, unit_length, counter, s.clone());
        }
        counter += len;
    }
}

fn melody_font(
    mut pcm: &mut wave_utils::MonoPcm,
    unit_length: usize,
    counter: usize,
    pitch: f64,
    len: usize,
) {
    square_wave(
        &mut pcm,
        pitch,
        0.1,
        unit_length * counter,
        unit_length * len * 7 / 8,
    )
}

fn base_font(
    mut pcm: &mut wave_utils::MonoPcm,
    unit_length: usize,
    counter: usize,
    pitch: f64,
    len: usize,
) {
    triangle_wave(
        &mut pcm,
        pitch,
        0.2,
        unit_length * counter,
        unit_length * len * 7 / 8,
    )
}

#[derive(Copy, Clone, Debug)]
enum Addendum {
    Short,
    Long,
}

fn percussion_font(
    mut pcm: &mut wave_utils::MonoPcm,
    unit_length: usize,
    counter: usize,
    addendum: Addendum,
) {
    white_noise(
        &mut pcm,
        0.1,
        unit_length * counter,
        match addendum {
            Addendum::Short => unit_length / 8,
            Addendum::Long => unit_length / 2,
        },
    );
}

pub fn second() {
    let unit_length = 9320; // number of samples per an eighth note
    let pcm_fs = 44100; /* 標本化周波数 */
    let pcm_length = unit_length * 128; /* 音データの長さ */
    let mut pcm = MonoPcm::new16(pcm_fs, pcm_length); /* 音データ */

    use self::NoteName::*;
    let melody_1 = vec![
        (1, None),
        (1, Some(Pitch(E, 5))),
        (1, None),
        (1, Some(Pitch(E, 5))),
        (1, None),
        (1, Some(Pitch(C, 5))),
        (1, None),
        (1, Some(Pitch(C, 5))),
        // --------------------
        (1, None),
        (1, Some(Pitch(A, 4))),
        (1, None),
        (1, Some(Pitch(A, 4))),
        (3, Some(Pitch(B, 4))),
        (1, None),
        // --------------------
        (1, None),
        (1, Some(Pitch(E, 5))),
        (1, None),
        (1, Some(Pitch(E, 5))),
        (1, None),
        (1, Some(Pitch(C, 5))),
        (1, None),
        (1, Some(Pitch(C, 5))),
        // --------------------
        (1, None),
        (1, Some(Pitch(F, 5))),
        (1, Some(Pitch(E, 5))),
        (1, Some(Pitch(D, 5))),
        (1, Some(Pitch(C, 5))),
        (1, None),
        (1, Some(Pitch(B, 4))),
        (1, None),
        // --------------------
        (1, None),
        (1, Some(Pitch(E, 5))),
        (1, None),
        (1, Some(Pitch(E, 5))),
        (1, None),
        (1, Some(Pitch(C, 5))),
        (1, None),
        (1, Some(Pitch(C, 5))),
        // --------------------
        (1, None),
        (1, Some(Pitch(A, 4))),
        (1, None),
        (1, Some(Pitch(A, 4))),
        (3, Some(Pitch(B, 4))),
        (1, None),
        // --------------------
        (1, None),
        (1, Some(Pitch(E, 5))),
        (1, None),
        (1, Some(Pitch(E, 5))),
        (1, Some(Pitch(F, 5))),
        (1, Some(Pitch(E, 5))),
        (1, Some(Pitch(D, 5))),
        (1, Some(Pitch(C, 5))),
        // --------------------
        (1, Some(Pitch(B, 4))),
        (1, Some(Pitch(G, 4))),
        (1, None),
        (1, Some(Pitch(G, 4))),
        (1, Some(Pitch(A, 4))),
        (1, Some(Pitch(B, 4))),
        (1, Some(Pitch(C, 5))),
        (1, None),
    ];

    let melody_2 = melody_1.iter().map(|(len, pitch)| {
        (*len, pitch.map(|Pitch(name, octave)| Pitch(name, octave - 1)))
    }).collect::<Vec<_>>();
    render_pitched(
        &mut pcm,
        unit_length,
        melody_font,
        &[&melody_1[..], &melody_2[..]].concat(),
    );

    render_pitched(
        &mut pcm,
        unit_length,
        base_font,
        &vec![
            (8, Some(Pitch(G, 3))),
            (8, Some(Pitch(F, 3))),
            (8, Some(Pitch(G, 3))),
            (8, Some(Pitch(D, 3))),
            // ===================
            (8, Some(Pitch(E, 3))),
            (8, Some(Pitch(A, 2))),
            (8, Some(Pitch(G, 2))),
            (2, Some(Pitch(A, 2))),
            (2, Some(Pitch(B, 2))),
            (4, Some(Pitch(C, 3))),
            // ===================
            (4, Some(Pitch(G, 3))),
            (4, Some(Pitch(G, 3))),
            (4, Some(Pitch(F, 3))),
            (4, Some(Pitch(F, 3))),
            (4, Some(Pitch(G, 3))),
            (4, Some(Pitch(G, 3))),
            (4, Some(Pitch(D, 3))),
            (4, Some(Pitch(D, 3))),
            // ===================
            (4, Some(Pitch(E, 3))),
            (4, Some(Pitch(E, 3))),
            (4, Some(Pitch(A, 2))),
            (4, Some(Pitch(A, 2))),
            (4, Some(Pitch(G, 2))),
            (4, Some(Pitch(G, 2))),
            (2, Some(Pitch(A, 2))),
            (2, Some(Pitch(B, 2))),
            (4, Some(Pitch(C, 3))),
        ],
    );

    let rhythm_1 = vec![
        (1, Some(Addendum::Long)),
        (1, None),
        (1, Some(Addendum::Short)),
        (1, Some(Addendum::Long)),
        (1, None),
        (1, Some(Addendum::Short)),
        (1, Some(Addendum::Long)),
        (1, None),
        // --------------------
        (1, Some(Addendum::Long)),
        (1, None),
        (1, None),
        (1, Some(Addendum::Long)),
        (1, None),
        (1, Some(Addendum::Short)),
        (1, Some(Addendum::Short)),
        (1, Some(Addendum::Short)),
    ];

    let rhythm_2 = vec![
        (1, Some(Addendum::Long)),
        (1, None),
        (1, Some(Addendum::Short)),
        (1, Some(Addendum::Long)),
        (1, None),
        (1, Some(Addendum::Short)),
        (1, Some(Addendum::Long)),
        (1, None),
        // --------------------
        (1, Some(Addendum::Long)),
        (1, Some(Addendum::Short)),
        (1, Some(Addendum::Long)),
        (1, Some(Addendum::Short)),
        (1, None),
        (1, Some(Addendum::Short)),
        (1, Some(Addendum::Long)),
        (1, None),
    ];

    render_pitchless(
        &mut pcm,
        unit_length,
        percussion_font,
        &[
            &rhythm_1[..],
            &rhythm_1[..],
            &rhythm_1[..],
            &rhythm_2[..],
            &rhythm_1[..],
            &rhythm_1[..],
            &rhythm_1[..],
            &rhythm_2[..],
        ]
        .concat(),
    );
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
