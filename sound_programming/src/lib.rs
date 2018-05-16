extern crate libc;
use std::slice::from_raw_parts_mut;
use std::f64::consts::PI;
pub use libc::c_char;
pub use libc::c_int;
pub use libc::c_double;
use std::ffi::CString;
pub mod fft;
pub mod filter;
pub mod wave;

#[allow(non_snake_case)]
pub fn safe_ADSR(e: &mut [c_double], A: usize, D: usize, S: c_double, R: usize, gate: usize, duration: usize){
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
  
    if R != 0
    {
        for n in gate..duration
        {
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

	/*pub*/ //fn wave_read_16bit_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
	/*pub*/ fn wave_write_16bit_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
	/*pub*/ //fn wave_read_16bit_stereo(pcm: *mut STEREO_PCM, file_name: *const c_char);
	/*pub*/ fn wave_write_16bit_stereo(pcm: *mut STEREO_PCM, file_name: *const c_char);
	
	pub fn wave_read_IMA_ADPCM_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
	pub fn wave_write_IMA_ADPCM_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
}



#[allow(non_snake_case)]
pub unsafe fn Hanning_window(w: *mut c_double, N: c_int)
{ 
  let w_slice = from_raw_parts_mut(w, N as usize);
  safe_Hanning_window(w_slice);
}


#[allow(non_snake_case)]
pub fn safe_Hanning_window(w_slice: &mut [c_double])
{ 
	let N = w_slice.len();
    for n in 0..N as usize {
	  	w_slice[n] = 0.5 - 0.5 * (
	  		2.0 * PI * 
	  		(n as f64 + if N % 2 == 0 { 0.0 } else { 0.5 })
	  		/ (N as f64)
	  	).cos();
  	}
}






pub fn to_c_str(a: &str) -> *mut i8 {
	CString::new(a).unwrap().into_raw()
}

#[allow(non_snake_case)]
pub fn wave_write_16bit_mono_safer2(path: &str, x : (&mut [f64], usize, i32, usize) ){
	let mut pcm1 : MONO_PCM = MONO_PCM{
			fs : x.1 as i32,
			bits : x.2,
			length : x.3 as i32,
			s : x.0.as_mut_ptr()
		};
	unsafe{
  		wave_write_16bit_mono(&mut pcm1, to_c_str(path)); 
	}
}

#[allow(non_snake_case)]
pub fn wave_write_16bit_stereo_safer2(path: &str, x : (&mut [f64], &mut [f64], usize, i32, usize) ){
	let mut pcm1 : STEREO_PCM = STEREO_PCM{
			fs : x.2 as i32,
			bits : x.3,
			length : x.4 as i32,
			sL : x.0.as_mut_ptr(),
			sR : x.1.as_mut_ptr()
		};
	unsafe{
  		wave_write_16bit_stereo(&mut pcm1, to_c_str(path)); 
	}
}
