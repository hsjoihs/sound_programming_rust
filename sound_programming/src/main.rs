extern crate sound_programming;
use sound_programming::wave_write_16bit_mono;
use std::f64::consts::PI;
use sound_programming::c_double;
use sound_programming::MONO_PCM;
use sound_programming::sinc;
fn main() {
	ex2_1();
	ex2_2();
	unsafe{
    	assert_eq!(sinc(2.1),2.1f64.sin()/2.1 );
    }
}


fn ex2_1(){

	let pcm_fs : usize = 44100;
	let pcm_length : usize = pcm_fs * 1;
	let mut pcm_s : Vec<c_double> = vec![0.0  ; pcm_length];
	
	
	let a = 0.1; /* 振幅 */
	let f0 = 500.0; /* 周波数 */
	
	/* サイン波 */
	for n in 0..pcm_length {
		pcm_s[n] = a * (2.0 * PI * f0 * (n as f64) / (pcm_fs as f64)).sin();
	}

	let mut pcm : MONO_PCM = MONO_PCM{ 
		fs : pcm_fs as i32, /* 標本化周波数 */
		bits : 16, /* 量子化精度 */
		length : pcm_length as i32, /* 音データの長さ */
		s: pcm_s.as_mut_ptr()
	};
	unsafe {
		wave_write_16bit_mono(&mut pcm, "ex2_1.wav".as_ptr() as *const i8);
	} 	
}

// ex2_1.c:
/*
int main(void)
{
  MONO_PCM pcm;
  int n;
  double a, f0;
  
  pcm.fs = 44100; /* 標本化周波数 */
  pcm.bits = 16; /* 量子化精度 */
  pcm.length = pcm.fs * 1; /* 音データの長さ */
  pcm.s = calloc(pcm.length, sizeof(double)); /* 音データ */
  
  a = 0.1; /* 振幅 */
  f0 = 500.0; /* 周波数 */
  
  /* サイン波 */
  for (n = 0; n < pcm.length; n++)
  {
    pcm.s[n] = a * sin(2.0 * M_PI * f0 * n / pcm.fs);
  }
  
  wave_write_16bit_mono(&pcm, "ex2_1.wav");
  
  free(pcm.s);
  
  return 0;
}
*/

fn ex2_2(){
unsafe {

}
}

// ex2_2.c:
/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wave.h"

void sine_wave(MONO_PCM *pcm, double f0, double a, int offset, int duration)
{
  int n;
  double *s;
  
  s = calloc(duration, sizeof(double));
  
  /* サイン波 */
  for (n = 0; n < duration; n++)
  {
    s[n] = a * sin(2.0 * M_PI * f0 * n / pcm->fs);
  }
  
  /* フェード処理 */
  for (n = 0; n < pcm->fs*0.01; n++)
  {
    s[n] *= (double)n / (pcm->fs * 0.01);
    s[duration - n - 1] *= (double)n / (pcm->fs * 0.01);
  }
  
  for (n = 0; n < duration; n++)
  {
    pcm->s[offset + n] += s[n];
  }
  
  free(s);
}

int main(void)
{
  MONO_PCM pcm;
  
  pcm.fs = 44100; /* 標本化周波数 */
  pcm.bits = 16; /* 量子化精度 */
  pcm.length = pcm.fs * 2; /* 音データの長さ */
  pcm.s = calloc(pcm.length, sizeof(double)); /* 音データ */
  
  sine_wave(&pcm, 261.63, 0.1, pcm.fs * 0.00, pcm.fs * 0.25); /* C4 */
  sine_wave(&pcm, 293.66, 0.1, pcm.fs * 0.25, pcm.fs * 0.25); /* D4 */
  sine_wave(&pcm, 329.63, 0.1, pcm.fs * 0.50, pcm.fs * 0.25); /* E4 */
  sine_wave(&pcm, 349.23, 0.1, pcm.fs * 0.75, pcm.fs * 0.25); /* F4 */
  sine_wave(&pcm, 392.00, 0.1, pcm.fs * 1.00, pcm.fs * 0.25); /* G4 */
  sine_wave(&pcm, 440.00, 0.1, pcm.fs * 1.25, pcm.fs * 0.25); /* A4 */
  sine_wave(&pcm, 493.88, 0.1, pcm.fs * 1.50, pcm.fs * 0.25); /* B4 */
  sine_wave(&pcm, 523.25, 0.1, pcm.fs * 1.75, pcm.fs * 0.25); /* C5 */
  
  wave_write_16bit_mono(&pcm, "ex2_2.wav");
  
  free(pcm.s);
  
  return 0;
}
*/
