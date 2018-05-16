use std::slice::from_raw_parts;
use std::fs::File;
use std::io::Read;
use std::mem;


macro_rules! READ_ARR {
    ($fp: expr, $i : ident, $t: ty; $len: expr) => (
        unsafe{
            let mut buf = [0; mem::size_of::<$t>()*$len];
            $fp.read(&mut buf).unwrap();
            let slice = from_raw_parts(buf.as_ptr() as *const $t, $len);
            for i in 0..$len{
                $i[i] = slice[i];
            }
            
            
        }
    )
}

#[allow(non_snake_case)]
fn foo(file_name : &str) -> (File, i32, i16, i32) {
    let mut fp = File::open(file_name).expect("file not found");
    let mut riff_chunk_ID    = [0;4]; READ_ARR!(fp, riff_chunk_ID,   i8 ; 4);
    let mut riff_chunk_size  = [0;1]; READ_ARR!(fp, riff_chunk_size, i32; 1);
    let mut file_format_type = [0;4]; READ_ARR!(fp, file_format_type,i8 ; 4);
    let mut fmt_chunk_ID     = [0;4]; READ_ARR!(fp, fmt_chunk_ID    ,i8 ; 4);
    let mut fmt_chunk_size   = [0;1]; READ_ARR!(fp, fmt_chunk_size,  i32; 1);
    let mut wave_format_type = [0;1]; READ_ARR!(fp, wave_format_type,i16; 1);
    let mut channel          = [0;1]; READ_ARR!(fp, channel,         i16; 1);
    let mut samples_per_sec  = [0;1]; READ_ARR!(fp, samples_per_sec, i32; 1);
    let mut bytes_per_sec    = [0;1]; READ_ARR!(fp, bytes_per_sec,   i32; 1);
    let mut block_size       = [0;1]; READ_ARR!(fp, block_size,      i16; 1);
    let mut bits_per_sample  = [0;1]; READ_ARR!(fp, bits_per_sample, i16; 1);
    let mut data_chunk_ID    = [0;4]; READ_ARR!(fp, data_chunk_ID,   i8 ; 4);
    let mut data_chunk_size  = [0;1]; READ_ARR!(fp, data_chunk_size, i32; 1);



    return (fp, samples_per_sec[0], bits_per_sample[0], data_chunk_size[0]);
}

//not tested
pub fn wave_read_8bit_mono_safer2(path : &str) -> (Vec<f64>, usize, i32, usize){
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = foo(path);
    let pcm_length = data_chunk_size as usize;
    let mut pcm_s = vec![0.0; pcm_length];

    for n in 0..pcm_length {
        let mut data = [0;1]; READ_ARR!(fp, data, u8; 1);
        pcm_s[n] = (data[0] as f64 - 128.0) / 128.0; /* 音データを-1以上1未満の範囲に正規化する */
    }

    return (pcm_s, pcm_fs as usize, pcm_bits as i32, pcm_length as usize);
}

// not tested
#[allow(non_snake_case)]
pub fn wave_read_8bit_stereo_safer2(path : &str) -> (Vec<f64>, Vec<f64>, usize, i32, usize){
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = foo(path);
    let pcm_length = (data_chunk_size / 2) as usize; /* 音データの長さ */

    let mut pcm_sL = vec![0.0; pcm_length];
    let mut pcm_sR = vec![0.0; pcm_length];
    for n in 0..pcm_length {
        let mut data = [0;1]; READ_ARR!(fp, data, u8; 1);
        pcm_sL[n] = (data[0] as f64 - 128.0) / 128.0; /* 音データを-1以上1未満の範囲に正規化する */
        let mut data = [0;1]; READ_ARR!(fp, data, i16; 1);
        pcm_sR[n] = (data[0] as f64 - 128.0) / 128.0; /* 音データを-1以上1未満の範囲に正規化する */
    }
    return (pcm_sL, pcm_sR, pcm_fs as usize, pcm_bits as i32, pcm_length as usize);
}


pub fn wave_read_16bit_mono_safer2(path : &str) -> (Vec<f64>, usize, i32, usize){
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = foo(path);
    let pcm_length = (data_chunk_size / 2) as usize;
    let mut pcm_s = vec![0.0; pcm_length];

    for n in 0..pcm_length {
        let mut data = [0;1]; READ_ARR!(fp, data, i16; 1);
        pcm_s[n] = (data[0] as f64) / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
    }

    return (pcm_s, pcm_fs as usize, pcm_bits as i32, pcm_length as usize);
}



#[allow(non_snake_case)]
pub fn wave_read_16bit_stereo_safer2(path : &str) -> (Vec<f64>, Vec<f64>, usize, i32, usize){
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = foo(path);
    let pcm_length = (data_chunk_size / 4) as usize; /* 音データの長さ */

    let mut pcm_sL = vec![0.0; pcm_length];
    let mut pcm_sR = vec![0.0; pcm_length];
    for n in 0..pcm_length {
        let mut data = [0;1]; READ_ARR!(fp, data, i16; 1);
        pcm_sL[n] = data[0] as f64 / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
        let mut data = [0;1]; READ_ARR!(fp, data, i16; 1);
        pcm_sR[n] = data[0] as f64 / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
    }
    
    return (pcm_sL, pcm_sR, pcm_fs as usize, pcm_bits as i32, pcm_length as usize);
}

