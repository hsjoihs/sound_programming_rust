extern crate sound_programming;
extern crate rand;
//use std::io::Write;
use sound_programming::to_c_str;
use sound_programming::wave_read_16bit_mono_safer2;
use sound_programming::fft::safe_IFFT;
use sound_programming::fft::safe_FFT;
use sound_programming::safe_Hanning_window;
use std::slice::from_raw_parts;
use std::slice::from_raw_parts_mut;
use std::f64::consts::PI;
use std::mem;
use sound_programming::wave_write_16bit_mono;
use sound_programming::wave_write_16bit_stereo;
use sound_programming::wave_read_16bit_stereo;
use sound_programming::STEREO_PCM;
use sound_programming::c_int;
use sound_programming::c_double;
use sound_programming::MONO_PCM;
use sound_programming::sinc;
use rand::Rng;
//use std::io;
fn main() {
	ex1_1();
	ex1_2();
	ex2_1();
	ex2_2();
	ex3_1();
	ex3_2();
	ex3_3();
	ex3_4();
	// ex3_5(); //slow
	ex4_1();
	ex4_2();
	ex4_3();
	ex6_4();
    assert_eq!(sinc(2.1),2.1f64.sin()/2.1 );
}



fn ex1_1(){

	
		/* 音データの入力 */	

		let (pcm0_slice, pcm0_fs, pcm0_bits, pcm0_length) = wave_read_16bit_mono_safer2("ex1_1_a.wav");

		let mut pcm1_s : Vec<c_double> = (0..pcm0_length)
		  .map(|n| pcm0_slice[n as usize])/* 音データのコピー */
		  .collect();

		let mut pcm1 : MONO_PCM = MONO_PCM{
			fs : pcm0_fs, /* 標本化周波数 */
			bits : pcm0_bits, /* 量子化精度 */
			length : pcm0_length, /* 音データの長さ */
			s : pcm1_s.as_mut_ptr()  /* 音データ */
		};
	unsafe{
		wave_write_16bit_mono(&mut pcm1, to_c_str("ex1_1_b.wav")); /* 音データの出力 */
	}	

	// pcm0.s possibly leaks? (may be handled nicely by drop of pcm0_slice)

}

#[allow(non_snake_case)]
fn ex1_2(){
	unsafe{
		let mut pcm0 : STEREO_PCM = mem::uninitialized();
		wave_read_16bit_stereo(&mut pcm0, to_c_str("ex1_2_a.wav")); /* 音データの入力 */


		let pcm0_sliceL = from_raw_parts(pcm0.sL, pcm0.length as usize);
		let pcm0_sliceR = from_raw_parts(pcm0.sR, pcm0.length as usize);

		let mut pcm1_sL : Vec<c_double> = (0..pcm0.length)
		 .map(|n| pcm0_sliceL[n as usize])
		 .collect();

		let mut pcm1_sR : Vec<c_double> = (0..pcm0.length)
		 .map(|n| pcm0_sliceR[n as usize])
		 .collect();
		
		let mut pcm1 : STEREO_PCM = STEREO_PCM {
			fs: pcm0.fs,
			bits: pcm0.bits,
			length: pcm0.length,
			sL : pcm1_sL.as_mut_ptr(),
			sR : pcm1_sR.as_mut_ptr()
		};

		wave_write_16bit_stereo(&mut pcm1, to_c_str("ex1_2_b.wav")); /* 音データの出力 */
	}
	// pcm0.sL and pcm0.sR possibly leaks?
}


fn ex2_1(){

	let pcm_fs : usize = 44100; /* 標本化周波数 */
	let pcm_length : usize = pcm_fs * 1; /* 音データの長さ */
	
	
	let a = 0.1; /* 振幅 */
	let f0 = 500.0; /* 周波数 */
	
	/* サイン波 */
	let mut pcm_s : Vec<c_double> = (0..pcm_length)
		.map(|n| a * (2.0 * PI * f0 * (n as f64) / (pcm_fs as f64)).sin())
		.collect();

	let mut pcm : MONO_PCM = MONO_PCM{ 
		fs : pcm_fs as i32, /* 標本化周波数 */
		bits : 16, /* 量子化精度 */
		length : pcm_length as i32, /* 音データの長さ */
		s: pcm_s.as_mut_ptr()
	};
	unsafe {
		wave_write_16bit_mono(&mut pcm, to_c_str("ex2_1.wav"));
	} 	
}





unsafe fn sine_wave(pcm : *mut MONO_PCM, f0: c_double, a: c_double, offset: c_int, duration: c_int) {
	/* サイン波 */
	let mut s : Vec<c_double> = (0..duration)
	 .map(|n| (2.0 * PI * f0 * (n as f64) / ((*pcm).fs as f64)).sin() * a)
	 .collect();


	/* フェード処理 */
	for n in 0..((*pcm).fs as f64*0.01).ceil() as usize {
		s[n] *= n as c_double / ((*pcm).fs as f64 * 0.01);
    	s[duration as usize - n - 1] *= n as c_double / ((*pcm).fs as f64 * 0.01);
	}

	for n in 0..duration as usize {
		let mut slice = from_raw_parts_mut((*pcm).s, (*pcm).length as usize);
		slice[offset as usize + n] += s[n];
	}
}
fn ex2_2(){
	let pcm_fs : usize = 44100; /* 標本化周波数 */
	let pcm_length : usize = pcm_fs * 2; /* 音データの長さ */
	let mut pcm_s : Vec<c_double> = vec![0.0  ; pcm_length];

	let mut pcm : MONO_PCM = MONO_PCM{ 
		fs : pcm_fs as i32, /* 標本化周波数 */
		bits : 16, /* 量子化精度 */
		length : pcm_length as i32, /* 音データの長さ */
		s: pcm_s.as_mut_ptr()
	};
unsafe{
  sine_wave(&mut pcm, 261.63, 0.1, itdyi(pcm.fs, 0.00), itdyi(pcm.fs, 0.25)); /* C4 */
  sine_wave(&mut pcm, 293.66, 0.1, itdyi(pcm.fs, 0.25), itdyi(pcm.fs, 0.25)); /* D4 */
  sine_wave(&mut pcm, 329.63, 0.1, itdyi(pcm.fs, 0.50), itdyi(pcm.fs, 0.25)); /* E4 */
  sine_wave(&mut pcm, 349.23, 0.1, itdyi(pcm.fs, 0.75), itdyi(pcm.fs, 0.25)); /* F4 */
  sine_wave(&mut pcm, 392.00, 0.1, itdyi(pcm.fs, 1.00), itdyi(pcm.fs, 0.25)); /* G4 */
  sine_wave(&mut pcm, 440.00, 0.1, itdyi(pcm.fs, 1.25), itdyi(pcm.fs, 0.25)); /* A4 */
  sine_wave(&mut pcm, 493.88, 0.1, itdyi(pcm.fs, 1.50), itdyi(pcm.fs, 0.25)); /* B4 */
  sine_wave(&mut pcm, 523.25, 0.1, itdyi(pcm.fs, 1.75), itdyi(pcm.fs, 0.25)); /* C5 */
  
  
	wave_write_16bit_mono(&mut pcm, to_c_str("ex2_2.wav"));
  }

}

// int_times_double_yielding_int
fn itdyi (i : c_int, d: c_double) -> c_int {
	((i as c_double) * d) as c_int
}


fn ex3_1(){
	let f0 = 500.0; /* 基本周波数 */

	let pcm_fs = 44100;
	let pcm_length = pcm_fs * 1;
	let mut pcm_s : Vec<c_double> = vec![0.0  ; pcm_length];

	/* ノコギリ波 */
  	for i_ in 1..=44 {
  		let i = i_ as f64;
   		for n in 0..pcm_length {
      		pcm_s[n] += 1.0 / i * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).sin();
    	}
  	}
  
    let gain = 0.1; /* ゲイン */
    
    for n in 0..pcm_length{
        pcm_s[n] *= gain;
    }


	let mut pcm : MONO_PCM = MONO_PCM{
		fs : pcm_fs as i32,
		bits : 16,
		length : pcm_length as i32,
		s : pcm_s.as_mut_ptr()
	};

	unsafe{
		wave_write_16bit_mono(&mut pcm, to_c_str("ex3_1.wav"));
	}
}

fn ex3_2(){
	let f0 = 500.0; /* 基本周波数 */

	let pcm_fs = 44100;
	let pcm_length = pcm_fs * 1;
	let mut pcm_s : Vec<c_double> = vec![0.0  ; pcm_length];

	/* 矩形波 */
	for j in 0..22 {
		let i = (2*j+1) as f64;
	  	for n in 0..pcm_length{
	      pcm_s[n] += 1.0 / i * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).sin();
	    }
	}
  
    let gain = 0.1; /* ゲイン */
    
    for n in 0..pcm_length{
        pcm_s[n] *= gain;
    }


	let mut pcm : MONO_PCM = MONO_PCM{
		fs : pcm_fs as i32,
		bits : 16,
		length : pcm_length as i32,
		s : pcm_s.as_mut_ptr()
	};

	unsafe{
		wave_write_16bit_mono(&mut pcm, to_c_str("ex3_2.wav"));
	}
}


fn ex3_3(){
	let f0 = 500.0; /* 基本周波数 */

	let pcm_fs = 44100;
	let pcm_length = pcm_fs * 1;
	let mut pcm_s : Vec<c_double> = vec![0.0  ; pcm_length];

	/* 三角波 */
  	for j in 0..22 {
		let i = (2*j+1) as f64;
   		for n in 0..pcm_length {
      		pcm_s[n] += 1.0 / i / i * (PI * i / 2.0).sin() * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).sin();
  		}
  	}
  
    let gain = 0.1; /* ゲイン */
    
    for n in 0..pcm_length{
        pcm_s[n] *= gain;
    }


	let mut pcm : MONO_PCM = MONO_PCM{
		fs : pcm_fs as i32,
		bits : 16,
		length : pcm_length as i32,
		s : pcm_s.as_mut_ptr()
	};

	unsafe{
		wave_write_16bit_mono(&mut pcm, to_c_str("ex3_3.wav"));
	}
}


fn ex3_4(){
	let f0 = 500.0; /* 基本周波数 */

	let pcm_fs = 44100;
	let pcm_length = pcm_fs * 1;
	let mut pcm_s : Vec<c_double> = vec![0.0  ; pcm_length];

	/* コサイン波の重ね合わせによるノコギリ波 */
	for i in 1..=44{
		let i = i as f64;
	    for n in 0..pcm_length{
	      pcm_s[n] += 1.0 / i * (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64)).cos();
	    }
	 }
  
    let gain = 0.1; /* ゲイン */
    
    for n in 0..pcm_length{
        pcm_s[n] *= gain;
    }


	let mut pcm : MONO_PCM = MONO_PCM{
		fs : pcm_fs as i32,
		bits : 16,
		length : pcm_length as i32,
		s : pcm_s.as_mut_ptr()
	};

	unsafe{
		wave_write_16bit_mono(&mut pcm, to_c_str("ex3_4.wav"));
	}
}


fn ex3_5(){
	let f0 = 1.0; /* 基本周波数 */
  
  

	let pcm_fs = 44100;
	let pcm_length = pcm_fs * 1;
	let mut pcm_s : Vec<c_double> = vec![0.0  ; pcm_length];

	let mut rng = rand::thread_rng();
	/* 白色雑音 */
  for i in 1..=22050{
  	let theta: f64 = rng.gen_range(0.0, 2.0 * PI); 
  	if i % 441 == 0 {
  		println!("{} / 22050", i);
  	}
  	let i = i as f64;   
    for n in 0..pcm_length{
      pcm_s[n] += (2.0 * PI * i * f0 * (n as f64) / (pcm_fs as f64) + theta).sin();
    }
  }
  
    let gain = 0.001; /* ゲイン */
    
    for n in 0..pcm_length{
        pcm_s[n] *= gain;
    }


	let mut pcm : MONO_PCM = MONO_PCM{
		fs : pcm_fs as i32,
		bits : 16,
		length : pcm_length as i32,
		s : pcm_s.as_mut_ptr()
	};

	unsafe{
		wave_write_16bit_mono(&mut pcm, to_c_str("ex3_5.wav"));
	}
}

#[allow(non_snake_case)]
fn verify(X_real: Vec<c_double>, X_imag: Vec<c_double>){
	let N = 64;

	/* 周波数特性 */
	for k in 0..N {
		assert_close(X_real[k], 0.0);
		assert_close(X_imag[k], match k {
			4 => -16.0,
			60 => 16.0,
			_ => 0.0
		});
	}
}

#[allow(non_snake_case)]
fn ex4_1(){
	let (X_real, X_imag) = foo(Box::new(|_| 1.0));	
	verify(X_real, X_imag)
}

#[allow(non_snake_case)]
fn foo(func : Box<Fn(usize) -> f64>) -> (Vec<c_double>, Vec<c_double>){
	let N = 64;
	let mut x_real : Vec<c_double> = vec![0.0; N];
	let mut x_imag : Vec<c_double> = vec![0.0; N];
	let mut X_real : Vec<c_double> = vec![0.0; N];
	let mut X_imag : Vec<c_double> = vec![0.0; N];

		let (pcm_slice, _, _, _) = wave_read_16bit_mono_safer2("sine_500hz.wav");

		/* 波形 */
		for n in 0..N {
 		   x_real[n] = pcm_slice[n] * func(n); /* x(n)の実数部 */
 		   x_imag[n] = 0.0; /* x(n)の虚数部 */
		}
	
		/* DFT */
		for k_ in 0..N {
			let k = k_ as f64;
			for n_ in 0..N {
				let n = n_ as f64;
				let N = N as f64;
                let W_real = (2.0 * PI * k * n / N).cos();
                let W_imag = -(2.0 * PI * k * n / N).sin();
                X_real[k_] += W_real * x_real[n_] - W_imag * x_imag[n_]; /* X(k)の実数部 */
                X_imag[k_] += W_real * x_imag[n_] + W_imag * x_real[n_]; /* X(k)の虚数部 */
			}
		}

		
		(X_real, X_imag)
	
}

fn assert_close(a : f64, b: f64){
	assert!((a-b).abs() < 1.75e-4);
}

#[allow(non_snake_case)]
fn ex4_2(){
	let N = 64;
	let mut w : Vec<c_double> = vec![0.0; N];
  	safe_Hanning_window(&mut w); /* ハニング窓 */

  	let (X_real, X_imag) = foo(Box::new(move |n| w[n]));

  	for k in 0..N{
  		assert_close(X_real[k], 0.0);
  		assert_close(X_imag[k], match k {
			3 => 4.0,
			4 => -8.0,
			5 => 4.0,
			59 => -4.0,
			60 => 8.0,
			61 => -4.0,
			_ => 0.0
		});
  	}
}




#[allow(non_snake_case)]
fn ex4_3(){
	let N = 64;
	let mut x_real : Vec<c_double> = vec![0.0; N];
	let mut x_imag : Vec<c_double> = vec![0.0; N];
	let (pcm_slice, _, _, _) = wave_read_16bit_mono_safer2("sine_500hz.wav");
		

	/* 波形 */
	for n in 0..N {
 		x_real[n] = pcm_slice[n]; /* x(n)の実数部 */
 		x_imag[n] = 0.0; /* x(n)の虚数部 */
	}

	safe_FFT(&mut x_real, &mut x_imag);  /* FFTの計算結果はx_realとx_imagに上書きされる */
}


#[allow(non_snake_case)]
fn ex6_4(){
	
		
    	let (pcm0_s, pcm0_fs, pcm0_bits, pcm0_length) = wave_read_16bit_mono_safer2("sine_500hz_3500hz.wav");

		let N = 256; /* DFTのサイズ */

		let mut pcm1_s : Vec<c_double> = vec![0.0; pcm0_length as usize];

		
		let mut x_real : Vec<c_double> = vec![0.0; N];
    	let mut x_imag : Vec<c_double> = vec![0.0; N];
    	let mut y_real : Vec<c_double> = vec![0.0; N];
    	let mut y_imag : Vec<c_double> = vec![0.0; N];
    	let mut b_real : Vec<c_double> = vec![0.0; N];
    	let mut b_imag : Vec<c_double> = vec![0.0; N];


    	let mut w : Vec<c_double> = vec![0.0; N];
    	safe_Hanning_window(&mut w); /* ハニング窓 */

    	let number_of_frame = (pcm0_length as usize - N / 2) / (N / 2); /* フレームの数 */

  		for frame  in 0..number_of_frame {
    		let offset = N / 2 * frame;
    
    		/* X(n) */
    		for n in 0..N {
      			x_real[n] = pcm0_s[offset + n] * w[n];
      			x_imag[n] = 0.0;
    		}
   			safe_FFT(&mut x_real, &mut x_imag);
    
    		/* B(k) */
    		let fe = 1000.0 / pcm0_fs as f64; /* エッジ周波数 */
    		let fe = (fe * N as f64) as usize;
    		for k in 0..=fe {
      			b_real[k] = 1.0;
      			b_imag[k] = 0.0;
    		}
    		for k in (fe+1)..= N / 2 {
      			b_real[k] = 0.0;
      			b_imag[k] = 0.0;
    		}
    		for k in 1..N / 2   {
      			b_real[N - k] = b_real[k];
      			b_imag[N - k] = -b_imag[k];
    		}
    
    		/* フィルタリング */
    		for k in 0..N {
      			y_real[k] = x_real[k] * b_real[k] - x_imag[k] * b_imag[k];
      			y_imag[k] = x_imag[k] * b_real[k] + x_real[k] * b_imag[k];
    		}
    		safe_IFFT(&mut y_real, &mut y_imag);
    
    		/* オーバーラップアド */
    		for n in 0..N {
      			pcm1_s[offset + n] += y_real[n];
    		}
  		}
  	
  		let mut pcm1 : MONO_PCM = MONO_PCM{
			fs : pcm0_fs,
			bits : pcm0_bits,
			length : pcm0_length,
			s :pcm1_s.as_mut_ptr()
		};
	unsafe{
  		wave_write_16bit_mono(&mut pcm1, to_c_str("ex6_4.wav")); 
  
	}

}

