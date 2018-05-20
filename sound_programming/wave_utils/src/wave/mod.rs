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
        pcm_s[n] = data.convert_to_float(); /* 音データを-1以上1未満の範囲に正規化する */
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
        pcm_sL[n] = data.convert_to_float(); /* 音データを-1以上1未満の範囲に正規化する */
        let data = fp.read_u8().unwrap();
        pcm_sR[n] = data.convert_to_float(); /* 音データを-1以上1未満の範囲に正規化する */
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
        pcm_s[n] = data.convert_to_float(); /* 音データを-1以上1未満の範囲に正規化する */
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
        pcm_sL[n] = data.convert_to_float(); /* 音データを-1以上1未満の範囲に正規化する */
        let data = fp.read_i16::<LittleEndian>().unwrap();
        pcm_sR[n] = data.convert_to_float(); /* 音データを-1以上1未満の範囲に正規化する */
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
        fp.write_u8(WaveData::convert_from_float(pcm.s[n]))
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
        fp.write_u8(WaveData::convert_from_float(pcm.s_l[n]))
            .unwrap(); /* 音データ（Lチャンネル）の書き出し */
        fp.write_u8(WaveData::convert_from_float(pcm.s_r[n]))
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
        let PCMU(dat) = WaveData::convert_from_float(pcm.s[n]);
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
        let PCMA(dat) = WaveData::convert_from_float(pcm.s[n]);
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
        fp.write_i16::<LittleEndian>(WaveData::convert_from_float(pcm.s_l[n]))
            .unwrap(); /* 音データ（Lチャンネル）の書き出し */
        fp.write_i16::<LittleEndian>(WaveData::convert_from_float(pcm.s_r[n]))
            .unwrap(); /* 音データ（Rチャンネル）の書き出し */
    }
}

pub fn wave_write_16bit_mono_safer3(path: &str, pcm: &MonoPcm) {
    let mut fp = wave_write_header::<MonoPcm, i16>(path, pcm);
    for n in 0..pcm.get_length() {
        fp.write_i16::<LittleEndian>(WaveData::convert_from_float(pcm.s[n]))
            .unwrap(); /* 音データの書き出し */
    }
}


impl WaveData for PCMU {
    fn convert_to_float(&self) -> f64 {
        unimplemented!();
    }
    const MYSTERIOUS: i32 = 50;
    const BYTE_NUM: i32 = 1;
    const CHUNK_SIZE: i32 = 18;
    const WAVE_FORMAT_TYPE: i16 = 7;
    fn convert_from_float(d: f64) -> Self {
        let mut x: f64;
        let s: i16; /* 16bitの音データ */
        let c: u8; /* 8bitの圧縮データ */
        let sign: u8;
        let mut exponent: u8;
        let mantissa: u8;
        let mut magnitude: i32;
        let level: [i16; 8] = [
            0x00FF, 0x01FF, 0x03FF, 0x07FF, 0x0FFF, 0x1FFF, 0x3FFF, 0x7FFF
        ];

        x = (d + 1.0) / 2.0 * 65536.0;

        if x > 65535.0 {
            x = 65535.0; /* クリッピング */
        } else if x < 0.0 {
            x = 0.0; /* クリッピング */
        }

        s = ((x + 0.5) as i32 - 32768) as i16; /* 四捨五入とオフセットの調節 */

        if s < 0 {
            magnitude = -s as i32;
            sign = 0x80;
        } else {
            magnitude = s as i32;
            sign = 0x00;
        }

        magnitude += 0x84;
        if magnitude > 32767 {
            magnitude = 0x7FFF;
        }
        exponent = 0;
        while exponent < 8 {
            if magnitude <= level[exponent as usize] as i32 {
                break;
            }
            exponent += 1;
        }

        mantissa = ((magnitude >> (exponent + 3)) & 0x0F) as u8;

        c = !(sign | (exponent << 4) | mantissa);

        PCMU(c) /* 圧縮データの書き出し */
    }
}

pub trait WaveData {
    fn convert_from_float(d: f64) -> Self;
    fn convert_to_float(&self) -> f64;
    const BYTE_NUM: i32;
    const MYSTERIOUS: i32;
    const CHUNK_SIZE: i32;
    const WAVE_FORMAT_TYPE: i16;
}

impl WaveData for u8 {
    fn convert_to_float(&self) -> f64 {
        (*self as f64 - 128.0) / 128.0 /* 音データを-1以上1未満の範囲に正規化する */
    }
    fn convert_from_float(d: f64) -> u8 {
        let mut s = (d + 1.0) / 2.0 * 256.0;

        if s > 255.0 {
            s = 255.0; /* クリッピング */
        } else if s < 0.0 {
            s = 0.0; /* クリッピング */
        }

        ((s + 0.5) as i32) as u8 /* 四捨五入 */
    }
    const BYTE_NUM: i32 = 1;
    const MYSTERIOUS: i32 = 36;
    const CHUNK_SIZE: i32 = 16;
    const WAVE_FORMAT_TYPE: i16 = 1;
}

impl WaveData for i16 {
    fn convert_to_float(&self) -> f64{
        (*self as f64) / 32768.0
    }
    fn convert_from_float(d: f64) -> i16 {
        let mut s = (d + 1.0) / 2.0 * 65536.0;

        if s > 65535.0 {
            s = 65535.0; /* クリッピング */
        } else if s < 0.0 {
            s = 0.0; /* クリッピング */
        }

        ((s + 0.5) as i32 - 32768) as i16 /* 四捨五入とオフセットの調節 */
    }
    const BYTE_NUM: i32 = 2;
    const MYSTERIOUS: i32 = 36;
    const CHUNK_SIZE: i32 = 16;
    const WAVE_FORMAT_TYPE: i16 = 1;
}

impl WaveData for PCMA {
    fn convert_to_float(&self) -> f64 {
        unimplemented!();
    }
    const MYSTERIOUS: i32 = 50;
    const BYTE_NUM: i32 = 1;
    const CHUNK_SIZE: i32 = 18;
    const WAVE_FORMAT_TYPE: i16 = 6;
    fn convert_from_float(d: f64) -> PCMA {
        let mut x: f64 = (d + 1.0) / 2.0 * 65536.0;
        let level: [i16; 8] = [
            0x00FF, 0x01FF, 0x03FF, 0x07FF, 0x0FFF, 0x1FFF, 0x3FFF, 0x7FFF
        ];

        if x > 65535.0 {
            x = 65535.0; /* クリッピング */
        } else if x < 0.0 {
            x = 0.0; /* クリッピング */
        }

        let s = ((x + 0.5) as i32 - 32768) as i16; /* 四捨五入とオフセットの調節 */

        let (mut magnitude, sign): (i32, u8) = if s < 0 {
            (-s as i32, 0x80)
        } else {
            (s as i32, 0x00)
        };

        if magnitude > 32767 {
            magnitude = 0x7FFF;
        }

        let mut exponent = 0 as u8;
        while exponent < 8 {
            if magnitude <= level[exponent as usize] as i32 {
                break;
            }
            exponent += 1;
        }

        let mantissa: u8 = if exponent == 0 {
            (magnitude >> 4) & 0x0F
        } else {
            (magnitude >> (exponent + 3)) & 0x0F
        } as u8;

        PCMA((sign | (exponent << 4) | mantissa) ^ 0xD5)
    }
}

