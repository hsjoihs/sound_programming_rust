extern crate num_complex;
extern crate rand;
extern crate sound_programming;
extern crate wave_utils;
use sound_programming::first::first;
use sound_programming::second::second;
use sound_programming::third::third;
use wave_utils::wave::wave_read_16bit_mono_safer3;
use wave_utils::wave::wave_read_16bit_stereo_safer3;
use wave_utils::wave::wave_read_8bit_mono_safer3;
use wave_utils::wave::wave_read_8bit_stereo_safer3;
use wave_utils::wave::wave_read_IMA_ADPCM_mono_safer3;
use wave_utils::wave::wave_read_PCMA_mono_safer3;
use wave_utils::wave::wave_read_PCMU_mono_safer3;
use wave_utils::wave::wave_write_16bit_mono_safer3;
use wave_utils::wave::wave_write_16bit_stereo_safer3;
use wave_utils::wave::wave_write_8bit_mono_safer3;
use wave_utils::wave::wave_write_8bit_stereo_safer3;
use wave_utils::wave::wave_write_IMA_ADPCM_mono_safer3;
use wave_utils::wave::wave_write_PCMA_mono_safer3;
use wave_utils::wave::wave_write_PCMU_mono_safer3;

fn main() {
    if false {
        first();
    }
    if true {
        second();
    }
    if false {
        third();
    }
    if false {
        ex11_7();
        ex11_8();
        ex11_9();
        eightbit();
    }
}

fn eightbit() {
    {
        let pcm0 = wave_read_16bit_mono_safer3("ex1_1_a.wav"); /* 音データの入力 */
        let mut pcm1 = pcm0; /* 音データのムーブ */
        pcm1.bits = 8;

        wave_write_8bit_mono_safer3("ex1_1_8bit.wav", &pcm1); /* 音データの出力 */
    }
    {
        let pcm1 = wave_read_8bit_mono_safer3("ex1_1_8bit.wav");
        let mut pcm2 = pcm1;
        pcm2.bits = 16;
        wave_write_16bit_mono_safer3("ex1_1_c.wav", &pcm2); /* 音データの出力 */
    }
    {
        let pcm0 = wave_read_16bit_stereo_safer3("ex1_2_a.wav");
        let mut pcm1 = pcm0; /* 音データのムーブ */
        pcm1.bits = 8;

        wave_write_8bit_stereo_safer3("ex1_2_8bit.wav", &pcm1);
    }
    {
        let pcm1 = wave_read_8bit_stereo_safer3("ex1_2_8bit.wav");
        let mut pcm2 = pcm1;
        pcm2.bits = 16;
        wave_write_16bit_stereo_safer3("ex1_2_c.wav", &pcm2); /* 音データの出力 */
    }
}

#[allow(non_snake_case)]
fn ex11_7() {
    let pcm0 = wave_read_16bit_mono_safer3("vocal.wav");
    wave_write_PCMU_mono_safer3("pcmu.wav", &pcm0);
    let pcm1 = wave_read_PCMU_mono_safer3("pcmu.wav");
    wave_write_16bit_mono_safer3("ex11_7_pcm.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex11_8() {
    let pcm0 = wave_read_16bit_mono_safer3("vocal.wav");
    wave_write_PCMA_mono_safer3("pcma.wav", &pcm0);
    let pcm1 = wave_read_PCMA_mono_safer3("pcma.wav");
    wave_write_16bit_mono_safer3("ex11_8_pcm.wav", &pcm1);
}

#[allow(non_snake_case)]
fn ex11_9() {
    let pcm0 = wave_read_16bit_mono_safer3("vocal.wav");
    wave_write_IMA_ADPCM_mono_safer3("ima_adpcm.wav", &pcm0);
    let pcm1 = wave_read_IMA_ADPCM_mono_safer3("ima_adpcm.wav");
    wave_write_16bit_mono_safer3("ex11_9_pcm.wav", &pcm1);
}
