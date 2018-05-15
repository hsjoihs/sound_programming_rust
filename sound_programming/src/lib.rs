extern crate libc;
use std::slice::from_raw_parts_mut;
use std::f64::consts::PI;
pub use libc::c_char;
pub use libc::c_int;
pub use libc::c_double;

#[link(name = "adsr")]
extern {
    pub fn ADSR(e: *mut c_double, A: c_int, D: c_int, S: c_double, R: c_int, gate: c_int, duration: c_int);
}

#[link(name = "fft")]
extern {
	pub fn FFT(x_real: *mut c_double, x_imag: *mut c_double, N: c_int);
	pub fn IFFT(x_real: *mut c_double, x_imag: *mut c_double, N: c_int);
}

#[link(name = "fir_filter")]
extern {
	pub fn FIR_LPF(fe: c_double, J: c_int, b: *mut c_double, w: *mut c_double);
	pub fn FIR_HPF(fe: c_double, J: c_int, b: *mut c_double, w: *mut c_double);
	pub fn FIR_BPF(fe1: c_double, fe2: c_double, J: c_int, b: *mut c_double, w: *mut c_double);
	pub fn FIR_BEF(fe1: c_double, fe2: c_double, J: c_int, b: *mut c_double, w: *mut c_double);
	pub fn FIR_filtering(x: *mut c_double, y: *mut c_double, L: c_int, b: *mut c_double, J: c_int);
}

#[link(name = "iir_filter")]
extern {
	 pub fn IIR_LPF(fc  : c_double, Q:  c_double , a: *mut c_double, b: *mut c_double)
	;pub fn IIR_HPF(fc  : c_double, Q:  c_double , a: *mut c_double, b: *mut c_double)
	;pub fn IIR_BPF(fc1 : c_double, fc2: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_BEF(fc1 : c_double, fc2: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_resonator(fc: c_double, Q: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_notch(fc: c_double, Q: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_low_shelving(fc: c_double, Q: c_double,g: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_high_shelving(fc: c_double, Q: c_double,g: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_peaking(fc: c_double, Q: c_double,g: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_filtering(x: *mut c_double, y: *mut c_double, L: c_int, a: *mut c_double, b: *mut c_double, I: c_int, J: c_int);
}

pub fn sinc(x: c_double) -> c_double {  
  if x == 0.0 {
    1.0
  } else {
     x.sin() / x
  }
}

#[repr(C)]
pub struct MONO_PCM {
  pub fs: c_int, /* 標本化周波数 */
  pub bits: c_int, /* 量子化精度 */
  pub length: c_int, /* 音データの長さ */
  pub s: *mut c_double /* 音データ */
}

#[allow(non_snake_case)]
#[repr(C)]
pub struct STEREO_PCM {
  pub fs: c_int, /* 標本化周波数 */
  pub bits: c_int, /* 量子化精度 */
  pub length: c_int, /* 音データの長さ */
  pub sL: *mut c_double, /* 音データ（Lチャンネル） */
  pub sR: *mut c_double /* 音データ（Rチャンネル） */
}



#[link(name = "wave")]
extern {
	pub fn wave_read_8bit_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
	pub fn wave_write_8bit_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
	pub fn wave_read_8bit_stereo(pcm: *mut STEREO_PCM, file_name: *const c_char);
	pub fn wave_write_8bit_stereo(pcm: *mut STEREO_PCM, file_name: *const c_char);

	pub fn wave_read_16bit_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
	pub fn wave_write_16bit_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
	pub fn wave_read_16bit_stereo(pcm: *mut STEREO_PCM, file_name: *const c_char);
	pub fn wave_write_16bit_stereo(pcm: *mut STEREO_PCM, file_name: *const c_char);
	
	pub fn wave_read_IMA_ADPCM_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
	pub fn wave_write_IMA_ADPCM_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
}
/*
#[link(name = "window_function")]
extern {
	pub fn Hanning_window(w: *mut c_double, N: c_int);
}
*/
#[allow(non_snake_case)]
pub unsafe fn Hanning_window(w: *mut c_double, N: c_int)
{ 
  let w_slice = from_raw_parts_mut(w, N as usize);
  if N % 2 == 0 {/* Nが偶数のとき */
    for n in 0..N as usize {
      w_slice[n] = 0.5 - 0.5 * (2.0 * PI * (n as f64) / (N as f64)).cos();
    }
  } else { /* Nが奇数のとき */
    for n in 0..N as usize {
      w_slice[n] = 0.5 - 0.5 * (2.0 * PI * (n as f64 + 0.5) / (N as f64)).cos();
    }
  }
}


