extern crate num_complex;
extern crate rand;

use rand::Rng;
use wave_utils::MonoPcm;

#[derive(Copy, Clone, Debug)]
pub enum NoteName {
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
pub struct Pitch(pub NoteName, pub i32);
impl Pitch {
    pub fn to_hertz(self) -> f64 {
        let Pitch(name, octave) = self;
        504.997 * 2f64.powf(1. / 12.).powi(name as i32 + 12 * (octave - 5))
    }
}

pub fn render_pitched<F>(
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

pub fn render_pitchless<F, A>(
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

pub fn melody_font(
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

pub fn base_font(
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
pub enum Addendum {
    Short,
    Long,
}

pub fn percussion_font(
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
