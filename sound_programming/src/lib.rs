extern crate num_complex;
extern crate rand;
extern crate wave_utils;
use std::f64::consts::PI;
use wave_utils::MonoPcm;

pub mod first;
pub mod second;
pub mod third;

fn sine_wave(pcm: &mut MonoPcm, f0: f64, a: f64, offset: usize, duration: usize) {
    /* サイン波 */
    let mut s: Vec<f64> = (0..duration)
        .map(|n| (2.0 * PI * f0 * (n as f64) / (pcm.fs as f64)).sin() * a)
        .collect();

    /* フェード処理 */
    for n in 0..(pcm.fs as f64 * 0.01).ceil() as usize {
        s[n] *= n as f64 / (pcm.fs as f64 * 0.01);
        s[duration - n - 1] *= n as f64 / (pcm.fs as f64 * 0.01);
    }

    for n in 0..duration as usize {
        pcm.s[offset + n] += s[n];
    }
}
