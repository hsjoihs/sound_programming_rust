extern crate byteorder;
use MonoPcm;
use StereoPcm;
use std::fs::File;
use self::byteorder::{LittleEndian, ReadBytesExt};



fn read_i8x4<T>(mut fp: T) -> [i8; 4] where T : byteorder::ReadBytesExt{
    let mut arr = [0; 4];
    for item in arr.iter_mut() {
       *item = fp.read_i8().unwrap();
    }
    return arr
}

#[allow(non_snake_case)]
fn read_header(file_name: &str) -> (File, i32, i16, i32) {
    let mut fp = File::open(file_name).expect("file not found");

    let _riff_chunk_ID = read_i8x4(&mut fp);
    let _riff_chunk_size = fp.read_i32::<LittleEndian>().unwrap();
    let _file_format_type = read_i8x4(&mut fp);
    let _fmt_chunk_ID = read_i8x4(&mut fp);
    let _fmt_chunk_size = fp.read_i32::<LittleEndian>().unwrap();
    let _wave_format_type = fp.read_i16::<LittleEndian>().unwrap();
    let _channel = fp.read_i16::<LittleEndian>().unwrap();
    let samples_per_sec = fp.read_i32::<LittleEndian>().unwrap();
    let _bytes_per_sec = fp.read_i32::<LittleEndian>().unwrap();
    let _block_size = fp.read_i16::<LittleEndian>().unwrap();
    let bits_per_sample = fp.read_i16::<LittleEndian>().unwrap();
    let _data_chunk_ID = read_i8x4(&mut fp);
    let data_chunk_size = fp.read_i32::<LittleEndian>().unwrap();

    return (
        fp,
        samples_per_sec,
        bits_per_sample,
        data_chunk_size,
    );
}

//not tested
pub fn wave_read_8bit_mono_safer3(path: &str) -> MonoPcm {
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = read_header(path);
    let pcm_length = data_chunk_size as usize;
    let mut pcm_s = vec![0.0; pcm_length];

    for n in 0..pcm_length {
        let data = fp.read_u8().unwrap();
        pcm_s[n] = (data as f64 - 128.0) / 128.0; /* 音データを-1以上1未満の範囲に正規化する */
    }

    return MonoPcm {
        s: pcm_s,
        fs: pcm_fs as usize,
        bits: pcm_bits as i32,
        length: pcm_length as usize,
    };
}

// not tested
#[allow(non_snake_case)]
pub fn wave_read_8bit_stereo_safer3(path: &str) -> StereoPcm {
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = read_header(path);
    let pcm_length = (data_chunk_size / 2) as usize; /* 音データの長さ */

    let mut pcm_sL = vec![0.0; pcm_length];
    let mut pcm_sR = vec![0.0; pcm_length];
    for n in 0..pcm_length {
        let data = fp.read_u8().unwrap();
        pcm_sL[n] = (data as f64 - 128.0) / 128.0; /* 音データを-1以上1未満の範囲に正規化する */
        let data = fp.read_u8().unwrap();
        pcm_sR[n] = (data as f64 - 128.0) / 128.0; /* 音データを-1以上1未満の範囲に正規化する */
    }
    return StereoPcm {
        s_l: pcm_sL,
        s_r: pcm_sR,
        fs: pcm_fs as usize,
        bits: pcm_bits as i32,
        length: pcm_length as usize,
    };
}

pub fn wave_read_16bit_mono_safer3(path: &str) -> MonoPcm {
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = read_header(path);
    let pcm_length = (data_chunk_size / 2) as usize;
    let mut pcm_s = vec![0.0; pcm_length];

    for n in 0..pcm_length {
        let data = fp.read_i16::<LittleEndian>().unwrap();
        pcm_s[n] = (data as f64) / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
    }

    return MonoPcm {
        s: pcm_s,
        fs: pcm_fs as usize,
        bits: pcm_bits as i32,
        length: pcm_length as usize,
    };
}

#[allow(non_snake_case)]
pub fn wave_read_16bit_stereo_safer3(path: &str) -> StereoPcm {
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = read_header(path);
    let pcm_length = (data_chunk_size / 4) as usize; /* 音データの長さ */

    let mut pcm_sL = vec![0.0; pcm_length];
    let mut pcm_sR = vec![0.0; pcm_length];
    for n in 0..pcm_length {
        let data = fp.read_i16::<LittleEndian>().unwrap();
        pcm_sL[n] = data as f64 / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
        let data = fp.read_i16::<LittleEndian>().unwrap();
        pcm_sR[n] = data as f64 / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
    }

    return StereoPcm {
        s_l: pcm_sL,
        s_r: pcm_sR,
        fs: pcm_fs as usize,
        bits: pcm_bits as i32,
        length: pcm_length as usize,
    };
}
