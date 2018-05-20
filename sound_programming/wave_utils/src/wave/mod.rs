extern crate byteorder;
use self::byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use MONO_PCM;
use MONO_PCM_CONST;
use MonoPcm;
use StereoPcm;
use libc::c_char;
use std::mem;
use std::slice::from_raw_parts_mut;
use to_c_str;
use wave::write::Pcm;
use wave::write::wave_write_header;

mod read;
mod write;

pub fn wave_read_8bit_mono_safer3(path: &str) -> MonoPcm {
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = read::read_header(path);
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

#[allow(non_snake_case)]
pub fn wave_read_8bit_stereo_safer3(path: &str) -> StereoPcm {
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = read::read_header(path);
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
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = read::read_header(path);
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
    let (mut fp, pcm_fs, pcm_bits, data_chunk_size) = read::read_header(path);
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

#[link(name = "wave")]
extern "C" {

    fn wave_read_PCMA_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
    fn wave_read_IMA_ADPCM_mono(pcm: *mut MONO_PCM, file_name: *const c_char);
    fn wave_write_IMA_ADPCM_mono(pcm: *const MONO_PCM_CONST, file_name: *const c_char);
    fn wave_read_PCMU_mono(pcm: *mut MONO_PCM, file_name: *const c_char);

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
pub fn wave_write_8bit_mono_safer3(path: &str, pcm: &MonoPcm) {
    let mut fp = wave_write_header::<MonoPcm, u8>(path, pcm);
    for n in 0..pcm.get_length() {
        fp.write_u8(write::WaveData::convert_from_float(pcm.s[n]))
            .unwrap(); /* 音データの書き出し */
    }
    if (pcm.length % 2) == 1 {
        /* 音データの長さが奇数のとき */

        fp.write_u8(0).unwrap(); /* 0パディング */
    }
}

#[allow(non_snake_case)]
pub fn wave_write_8bit_stereo_safer3(path: &str, pcm: &StereoPcm) {
    let mut fp = wave_write_header::<StereoPcm, u8>(path, pcm);
    for n in 0..pcm.length {
        fp.write_u8(write::WaveData::convert_from_float(pcm.s_l[n]))
            .unwrap(); /* 音データ（Lチャンネル）の書き出し */
        fp.write_u8(write::WaveData::convert_from_float(pcm.s_r[n]))
            .unwrap(); /* 音データ（Rチャンネル）の書き出し */
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
    let mut fp = wave_write_header::<MonoPcm, PCMU>(path, pcm);

    for n in 0..pcm.get_length() {
        let PCMU(dat) = write::WaveData::convert_from_float(pcm.s[n]);
        fp.write_u8(dat).unwrap(); /* 音データの書き出し */
    }
    if (pcm.length % 2) == 1 {
        /* 音データの長さが奇数のとき */

        fp.write_u8(0).unwrap(); /* 0パディング */
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

struct PCMA(u8);
struct PCMU(u8);

#[allow(non_snake_case)]
pub fn wave_write_PCMA_mono_safer3(path: &str, pcm: &MonoPcm) {
    let mut fp = wave_write_header::<MonoPcm, PCMA>(path, pcm);

    for n in 0..pcm.get_length() {
        let PCMA(dat) = write::WaveData::convert_from_float(pcm.s[n]);
        fp.write_u8(dat).unwrap(); /* 音データの書き出し */
    }
    if (pcm.length % 2) == 1 {
        /* 音データの長さが奇数のとき */

        fp.write_u8(0).unwrap(); /* 0パディング */
    }
}

#[allow(non_snake_case)]
pub fn wave_write_16bit_stereo_safer3(path: &str, pcm: &StereoPcm) {
    let mut fp = wave_write_header::<StereoPcm, i16>(path, pcm);
    for n in 0..pcm.length {
        fp.write_i16::<LittleEndian>(write::WaveData::convert_from_float(pcm.s_l[n]))
            .unwrap(); /* 音データ（Lチャンネル）の書き出し */
        fp.write_i16::<LittleEndian>(write::WaveData::convert_from_float(pcm.s_r[n]))
            .unwrap(); /* 音データ（Rチャンネル）の書き出し */
    }
}

pub fn wave_write_16bit_mono_safer3(path: &str, pcm: &MonoPcm) {
    let mut fp = wave_write_header::<MonoPcm, i16>(path, pcm);
    for n in 0..pcm.get_length() {
        fp.write_i16::<LittleEndian>(write::WaveData::convert_from_float(pcm.s[n]))
            .unwrap(); /* 音データの書き出し */
    }
}
