extern crate byteorder;
use self::byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use MONO_PCM_CONST;
use MonoPcm;
use StereoPcm;
use libc::c_char;
use to_c_str;
use wave::read::read_i8x4;
use wave::read::read_partial_header;
use wave::read::read_u8x4;
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

    fn wave_write_IMA_ADPCM_mono(pcm: *const MONO_PCM_CONST, file_name: *const c_char);

}

#[allow(non_upper_case_globals)]
const step_size_table: [i32; 89] = [
    7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19, 21, 23, 25, 28, 31, 34, 37, 41, 45, 50, 55, 60, 66,
    73, 80, 88, 97, 107, 118, 130, 143, 157, 173, 190, 209, 230, 253, 279, 307, 337, 371, 408, 449,
    494, 544, 598, 658, 724, 796, 876, 963, 1060, 1166, 1282, 1411, 1552, 1707, 1878, 2066, 2272,
    2499, 2749, 3024, 3327, 3660, 4026, 4428, 4871, 5358, 5894, 6484, 7132, 7845, 8630, 9493,
    10442, 11487, 12635, 13899, 15289, 16818, 18500, 20350, 22385, 24623, 27086, 29794, 32767,
];

#[allow(non_upper_case_globals)]
const index_table: [i32; 16] = [-1, -1, -1, -1, 2, 4, 6, 8, -1, -1, -1, -1, 2, 4, 6, 8];

#[allow(non_snake_case)]
pub fn wave_read_IMA_ADPCM_mono_safer3(path: &str) -> MonoPcm {
    let (mut fp, samples_per_sec, block_size, _bits_per_sample) = read_partial_header(path);
    let _extra_size = fp.read_i16::<LittleEndian>().unwrap();
    let samples_per_block = fp.read_i16::<LittleEndian>().unwrap();
    let _fact_chunk_ID = read_i8x4(&mut fp);
    let _fact_chunk_size = fp.read_i32::<LittleEndian>().unwrap();
    let sample_length = fp.read_i32::<LittleEndian>().unwrap();
    let _data_chunk_ID = read_i8x4(&mut fp);
    let data_chunk_size = fp.read_i32::<LittleEndian>().unwrap();

    let number_of_block: i32 = data_chunk_size / block_size as i32;

    let pcm_fs = samples_per_sec; /* 標本化周波数 */
    let pcm_length = sample_length; /* 音データの長さ */
    let mut pcm_s = vec![0.0; pcm_length as usize]; /* メモリの確保 */

    for block in 0..number_of_block {
        let mut s: i16;
        let mut index: i32 = 0;
        let offset: i32 = samples_per_block as i32 * block as i32;
        let mut c: u8;
        let mut data: u8 = 0;
        let mut step_size: i32;
        let mut dp: i32;
        let mut sp: i32 = 0;
        let mut header: [u8; 4];
        for n in 0..samples_per_block {
            if n == 0 {
                header = read_u8x4(&mut fp);
                sp = (((header[1] as i8) as i16) << 8) as i32 + header[0] as i32;
                index = header[2] as i32;
                s = sp as i16;
            } else {
                /* 4bitの圧縮データ */
                c = if (n % 2) == 1 {
                    data = fp.read_u8().unwrap(); /* 圧縮データの読み取り */

                    (data & 0x0F) as u8 /* dataの下位4bit */
                } else {
                    ((data >> 4) & 0x0F) as u8 /* dataの上位4bit */
                };

                step_size = step_size_table[index as usize];

                /* 伸張 */
                dp = step_size >> 3;
                if (c & 0x1) == 0x1 {
                    dp += step_size >> 2;
                }
                if (c & 0x2) == 0x2 {
                    dp += step_size >> 1;
                }
                if (c & 0x4) == 0x4 {
                    dp += step_size;
                }
                if (c & 0x8) == 0x8 {
                    sp -= dp;
                } else {
                    sp += dp;
                }

                if sp > 32767 {
                    sp = 32767;
                } else if sp < -32768 {
                    sp = -32768;
                }

                index += index_table[c as usize];

                if index < 0 {
                    index = 0;
                } else if index > 88 {
                    index = 88;
                }

                s = sp as i16;
            }

            pcm_s[offset as usize + n as usize] = s as f64 / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
        }
    }

    MonoPcm {
        fs: pcm_fs as usize,
        bits: 16,
        length: pcm_length as usize,
        s: pcm_s,
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
        fp.write_u8(WaveData::convert_from_float(pcm.s[n])).unwrap(); /* 音データの書き出し */
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
    let (mut fp, pcm_fs, _, data_chunk_size) = read::read_header2(path, true);
    let pcm_length = data_chunk_size as usize;
    let mut pcm_s = vec![0.0; pcm_length];

    for n in 0..pcm_length {
        let data = PCMU(fp.read_u8().unwrap());
        pcm_s[n] = data.convert_to_float(); /* 音データを-1以上1未満の範囲に正規化する */
    }

    return MonoPcm {
        s: pcm_s,
        fs: pcm_fs as usize,
        bits: 16 as i32,
        length: pcm_length as usize,
    };
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
    let (mut fp, pcm_fs, _, data_chunk_size) = read::read_header2(path, true);
    let pcm_length = data_chunk_size as usize;
    let mut pcm_s = vec![0.0; pcm_length];

    for n in 0..pcm_length {
        let data = PCMA(fp.read_u8().unwrap());
        pcm_s[n] = data.convert_to_float(); /* 音データを-1以上1未満の範囲に正規化する */
    }

    return MonoPcm {
        s: pcm_s,
        fs: pcm_fs as usize,
        bits: 16 as i32,
        length: pcm_length as usize,
    };
}

struct PCMA(u8);
struct PCMU(u8);

#[allow(non_camel_case_types)]
struct IMA_ADPCM {}

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
        let s: i16; /* 16bitの音データ */
        let magnitude: i32;

        let PCMU(mut c) = self; /* 8bitの圧縮データ */
        c = !c;

        let sign: u8 = c & 0x80;
        let exponent: u8 = (c >> 4) & 0x07;
        let mantissa: u8 = c & 0x0F;

        magnitude = ((((mantissa as i32) << 3) + 0x84) << exponent) - 0x84;

        if sign == 0x80 {
            s = -(magnitude as i16);
        } else {
            s = magnitude as i16;
        }

        s as f64 / 32768.0 /* 音データを-1以上1未満の範囲に正規化する */
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
    fn convert_to_float(&self) -> f64 {
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
        let PCMA(mut c) = *self; /* 8bitの圧縮データ */
        c ^= 0xD5;
        let sign: u8 = c & 0x80;
        let exponent: u8 = (c >> 4) & 0x07;
        let mantissa: u8 = c & 0x0F;
        let magnitude: i32 = if exponent == 0 {
            ((mantissa as i32) << 4) + 0x0008
        } else {
            (((mantissa as i32) << 4) + 0x0108) << (exponent - 1)
        };

        /* 16bitの音データ */
        let s: i16 = if sign == 0x80 {
            -(magnitude as i16)
        } else {
            magnitude as i16
        };

        s as f64 / 32768.0 /* 音データを-1以上1未満の範囲に正規化する */
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
