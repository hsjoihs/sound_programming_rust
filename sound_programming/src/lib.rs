extern crate libc;
pub use libc::c_char;
pub use libc::c_double;
pub use libc::c_int;
use std::f64::consts::PI;
use std::ffi::CString;
use std::mem;
use std::slice::from_raw_parts_mut;
pub mod fft;
pub mod filter;
pub mod wave;

#[allow(non_snake_case)]
pub fn safe_ADSR(
    e: &mut [c_double],
    A: usize,
    D: usize,
    S: c_double,
    R: usize,
    gate: usize,
    duration: usize,
) {
    if A != 0 {
        for n in 0..A {
            e[n] = 1.0 - (-5.0 * n as f64 / A as f64).exp();
        }
    }

    if D != 0 {
        for n in A..gate {
            e[n] = S + (1.0 - S) * (-5.0 * (n - A) as f64 / D as f64).exp();
        }
    } else {
        for n in A..gate {
            e[n] = S;
        }
    }

    if R != 0 {
        for n in gate..duration {
            e[n] = e[gate - 1] * (-5.0 * (n - gate + 1) as f64 / R as f64).exp();
        }
    }
}

pub fn sinc(x: c_double) -> c_double {
    if x == 0.0 {
        1.0
    } else {
        x.sin() / x
    }
}

pub struct MonoPcm {
    pub fs: usize,
    pub bits: i32,
    pub length: usize,
    pub s: Vec<f64>,
}

impl MonoPcm {
    pub fn new16(fs: usize, length: usize) -> Self {
        MonoPcm {
            fs,
            length,
            bits: 16,
            s: vec![0.0; length],
        }
    }
    pub fn blank_copy(orig: &Self) -> Self {
        MonoPcm {
            s: vec![0.0; orig.length as usize],
            ..*orig
        }
    }

    pub fn mult_constant_gain(&mut self, gain: f64) {
        for n in 0..self.length {
            self.s[n] *= gain;
        }
    }

    pub fn new16_fn(fs: usize, length: usize, mut fun: Box<FnMut(usize) -> f64>) -> Self {
        MonoPcm {
            fs,
            length,
            bits: 16,
            s: (0..length).map(|n| fun(n)).collect(),
        }
    }
}

pub struct StereoPcm {
    pub fs: usize,
    pub bits: i32,
    pub length: usize,
    pub s_l: Vec<f64>,
    pub s_r: Vec<f64>,
}

#[repr(C)]
pub struct MONO_PCM {
    pub fs: c_int,        /* 標本化周波数 */
    pub bits: c_int,      /* 量子化精度 */
    pub length: c_int,    /* 音データの長さ */
    pub s: *mut c_double, /* 音データ */
}

#[repr(C)]
pub struct MONO_PCM_CONST {
    pub fs: c_int,          /* 標本化周波数 */
    pub bits: c_int,        /* 量子化精度 */
    pub length: c_int,      /* 音データの長さ */
    pub s: *const c_double, /* 音データ */
}

#[allow(non_snake_case)]
#[repr(C)]
pub struct STEREO_PCM {
    pub fs: c_int,         /* 標本化周波数 */
    pub bits: c_int,       /* 量子化精度 */
    pub length: c_int,     /* 音データの長さ */
    pub sL: *mut c_double, /* 音データ（Lチャンネル） */
    pub sR: *mut c_double, /* 音データ（Rチャンネル） */
}

#[allow(non_snake_case)]
#[repr(C)]
pub struct STEREO_PCM_CONST {
    pub fs: c_int,           /* 標本化周波数 */
    pub bits: c_int,         /* 量子化精度 */
    pub length: c_int,       /* 音データの長さ */
    pub sL: *const c_double, /* 音データ（Lチャンネル） */
    pub sR: *const c_double, /* 音データ（Rチャンネル） */
}

#[link(name = "wave")]
extern "C" {
    pub fn wave_write_8bit_mono(pcm: *const MONO_PCM_CONST, file_name: *const c_char);
    pub fn wave_write_8bit_stereo(pcm: *const STEREO_PCM_CONST, file_name: *const c_char);

    fn wave_write_16bit_mono(pcm: *const MONO_PCM_CONST, file_name: *const c_char);
    fn wave_write_16bit_stereo(pcm: *const STEREO_PCM_CONST, file_name: *const c_char);

    fn wave_read_PCMA_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
    fn wave_write_PCMA_mono(pcm: *const MONO_PCM_CONST, file_name: *const c_char);
    fn wave_read_IMA_ADPCM_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
    fn wave_write_IMA_ADPCM_mono(pcm: *const MONO_PCM_CONST, file_name: *const c_char);
    fn wave_read_PCMU_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
    fn wave_write_PCMU_mono(pcm: *const MONO_PCM_CONST, file_name: *const c_char);

}

#[allow(non_snake_case)]
pub unsafe fn Hanning_window(w: *mut c_double, N: c_int) {
    let w_slice = from_raw_parts_mut(w, N as usize);
    safe_Hanning_window(w_slice);
}

#[allow(non_snake_case)]
pub fn safe_Hanning_window(w_slice: &mut [c_double]) {
    let N = w_slice.len();
    for n in 0..N as usize {
        w_slice[n] = 0.5 - 0.5 * (2.0 * PI * (n as f64 + if N % 2 == 0 { 0.0 } else { 0.5 })
            / (N as f64))
            .cos();
    }
}

pub fn to_c_str(a: &str) -> *mut i8 {
    CString::new(a).unwrap().into_raw()
}

#[allow(non_snake_case)]
pub fn wave_write_16bit_mono_safer2(path: &str, x: (&[f64], usize, i32, usize)) {
    let pcm1: MONO_PCM_CONST = MONO_PCM_CONST {
        fs: x.1 as i32,
        bits: x.2,
        length: x.3 as i32,
        s: x.0.as_ptr(),
    };
    unsafe {
        wave_write_16bit_mono(&pcm1, to_c_str(path));
    }
}

#[allow(non_snake_case)]
pub fn wave_write_16bit_mono_safer3(path: &str, pcm: &MonoPcm) {
    let pcm1: MONO_PCM_CONST = MONO_PCM_CONST {
        fs: pcm.fs as i32,
        bits: pcm.bits,
        length: pcm.length as i32,
        s: pcm.s.as_ptr(),
    };
    unsafe {
        wave_write_16bit_mono(&pcm1, to_c_str(path));
    }
}

#[allow(non_snake_case)]
pub fn wave_read_IMA_ADPCM_mono_safer3(path: &str) -> MonoPcm {
    unsafe {
        let mut pcm: MONO_PCM = mem::uninitialized();
        wave_read_IMA_ADPCM_mono(&mut pcm, to_c_str(path));
        MonoPcm {
            fs: pcm.fs as usize,
            bits: pcm.bits,
            length: pcm.length as usize,
            s: from_raw_parts_mut(pcm.s, pcm.length as usize).to_vec(),
        }
    }
}

#[allow(non_snake_case)]
pub fn wave_write_IMA_ADPCM_mono_safer3(path: &str, pcm: &MonoPcm) {
    let pcm1: MONO_PCM_CONST = MONO_PCM_CONST {
        fs: pcm.fs as i32,
        bits: pcm.bits,
        length: pcm.length as i32,
        s: pcm.s.as_ptr(),
    };
    unsafe {
        wave_write_IMA_ADPCM_mono(&pcm1, to_c_str(path));
    }
}

#[allow(non_snake_case)]
pub fn wave_read_PCMU_mono_safer3(path: &str) -> MonoPcm {
    unsafe {
        let mut pcm: MONO_PCM = mem::uninitialized();
        wave_read_PCMU_mono(&mut pcm, to_c_str(path));
        MonoPcm {
            fs: pcm.fs as usize,
            bits: pcm.bits,
            length: pcm.length as usize,
            s: from_raw_parts_mut(pcm.s, pcm.length as usize).to_vec(),
        }
    }
}

#[allow(non_snake_case)]
pub fn wave_write_PCMU_mono_safer3(path: &str, pcm: &MonoPcm) {
    let pcm1: MONO_PCM_CONST = MONO_PCM_CONST {
        fs: pcm.fs as i32,
        bits: pcm.bits,
        length: pcm.length as i32,
        s: pcm.s.as_ptr(),
    };
    unsafe {
        wave_write_PCMU_mono(&pcm1, to_c_str(path));
    }
}

#[allow(non_snake_case)]
pub fn wave_read_PCMA_mono_safer3(path: &str) -> MonoPcm {
    unsafe {
        let mut pcm: MONO_PCM = mem::uninitialized();
        wave_read_PCMA_mono(&mut pcm, to_c_str(path));
        MonoPcm {
            fs: pcm.fs as usize,
            bits: pcm.bits,
            length: pcm.length as usize,
            s: from_raw_parts_mut(pcm.s, pcm.length as usize).to_vec(),
        }
    }
}

#[allow(non_snake_case)]
pub fn wave_write_PCMA_mono_safer3(path: &str, pcm: &MonoPcm) {
    let pcm1: MONO_PCM_CONST = MONO_PCM_CONST {
        fs: pcm.fs as i32,
        bits: pcm.bits,
        length: pcm.length as i32,
        s: pcm.s.as_ptr(),
    };
    unsafe {
        wave_write_PCMA_mono(&pcm1, to_c_str(path));
    }
}

#[allow(non_snake_case)]
pub fn wave_write_16bit_stereo_safer2(path: &str, x: (&[f64], &[f64], usize, i32, usize)) {
    let pcm1: STEREO_PCM_CONST = STEREO_PCM_CONST {
        fs: x.2 as i32,
        bits: x.3,
        length: x.4 as i32,
        sL: x.0.as_ptr(),
        sR: x.1.as_ptr(),
    };
    unsafe {
        wave_write_16bit_stereo(&pcm1, to_c_str(path));
    }
}

#[allow(non_snake_case)]
pub fn wave_write_16bit_stereo_safer3(path: &str, pcm: &StereoPcm) {
    let pcm1: STEREO_PCM_CONST = STEREO_PCM_CONST {
        fs: pcm.fs as i32,
        bits: pcm.bits,
        length: pcm.length as i32,
        sL: pcm.s_l.as_ptr(),
        sR: pcm.s_r.as_ptr(),
    };
    unsafe {
        wave_write_16bit_stereo(&pcm1, to_c_str(path));
    }
}
