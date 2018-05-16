use sinc;
use std::f64::consts::PI;
use libc::c_double;
use libc::c_int;

#[link(name = "fir_filter")]
extern {
	//pub fn FIR_LPF(fe: c_double, J: c_int, b: *mut c_double, w: *mut c_double);
	pub fn FIR_HPF(fe: c_double, J: c_int, b: *mut c_double, w: *mut c_double);
	pub fn FIR_BPF(fe1: c_double, fe2: c_double, J: c_int, b: *mut c_double, w: *mut c_double);
	pub fn FIR_BEF(fe1: c_double, fe2: c_double, J: c_int, b: *mut c_double, w: *mut c_double);
	//*pub*/ fn FIR_filtering(x: *const c_double, y: *mut c_double, L: c_int, b: *mut c_double, J: c_int);
}

#[allow(non_snake_case)]
pub fn safe_FIR_LPF(fe: c_double, J: usize, b: &mut [c_double], w: &mut [c_double]){
	assert_eq!(J%2, 0);
	assert_eq!(J+1, b.len());
	assert_eq!(J+1, w.len());

	// k = m + J/2
	for (k, item) in b.iter_mut().enumerate() {
		*item = 2.0 * fe * sinc(2.0 * PI * fe * (k as f64 - (J/2) as f64));
	}

	for (m, item) in b.iter_mut().enumerate(){
		*item *= w[m];
	}

}

#[allow(non_snake_case)]
pub fn safe_FIR_filtering(x: &[c_double], y: &mut [c_double], L: usize, b: &mut [c_double], J: usize){
	// check index here
	assert_eq!(J+1, b.len());
	
	for n in 0..L {
		for m in 0..=J {
			if n >= m {
				y[n] += b[m] * x[n-m];
			}
		}
	}
}


#[link(name = "iir_filter")]
extern {
	 /*pub*/ fn IIR_LPF(fc  : c_double, Q:  c_double , a: *mut c_double, b: *mut c_double)
	;pub fn IIR_HPF(fc  : c_double, Q:  c_double , a: *mut c_double, b: *mut c_double)
	;pub fn IIR_BPF(fc1 : c_double, fc2: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_BEF(fc1 : c_double, fc2: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_resonator(fc: c_double, Q: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_notch(fc: c_double, Q: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_low_shelving(fc: c_double, Q: c_double,g: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_high_shelving(fc: c_double, Q: c_double,g: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_peaking(fc: c_double, Q: c_double,g: c_double, a: *mut c_double, b: *mut c_double)
	;pub fn IIR_filtering(x: *const c_double, y: *mut c_double, L: c_int, a: *const c_double, b: *const c_double, I: c_int, J: c_int);
}

#[allow(non_snake_case)]
pub fn safe_IIR_LPF(fc  : c_double, Q:  c_double , a: &mut [c_double], b: &mut [c_double]){
	unsafe{
		IIR_LPF(fc, Q, a.as_mut_ptr(), b.as_mut_ptr()); /* IIRフィルタの設計 */
	}
}