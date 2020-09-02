extern crate num_complex;
extern crate rand;
extern crate sound_programming;
extern crate wave_utils;
use wave_utils::wave::wave_write_16bit_mono_safer3;
use wave_utils::MonoPcm;
use sound_programming::second::*;
fn main() {
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
    wave_write_16bit_mono_safer3("rebaku.wav", &pcm);
}