use libc::c_double;
use libc::c_int;

#[link(name = "fft")]
extern {
	pub fn FFT(x_real: *mut c_double, x_imag: *mut c_double, N: c_int);
	pub fn IFFT(x_real: *mut c_double, x_imag: *mut c_double, N: c_int);
}

#[allow(unused_variables, non_snake_case)]
pub fn safe_FFT (x_real: &mut [c_double], x_imag: &mut [c_double]){
	let N = x_real.len() ;
	assert_eq!(N, x_imag.len());
	unsafe{
		FFT(x_real.as_mut_ptr(), x_imag.as_mut_ptr(), N as c_int);
	}
}
#[allow(unused_variables, non_snake_case)]
pub fn safe_IFFT(x_real: &mut [c_double], x_imag: &mut [c_double]){
	let N = x_real.len() ;
	assert_eq!(N, x_imag.len());
	unsafe{
		IFFT(x_real.as_mut_ptr(), x_imag.as_mut_ptr(), N as c_int);
	}
}