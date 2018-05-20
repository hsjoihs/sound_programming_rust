extern crate byteorder;
use self::byteorder::{LittleEndian, WriteBytesExt};
use MonoPcm;
use StereoPcm;
use std::fs::File;
use wave::PCMA;
use wave::PCMU;

fn write_i8x4<T>(mut fp: T, arr: [i8; 4])
where
    T: byteorder::WriteBytesExt,
{
    for item in arr.iter() {
        fp.write_i8(*item).unwrap();
    }
}

#[allow(non_snake_case)]
pub fn wave_write_header<T, U>(path: &str, pcm: &T) -> File
where
    T: Pcm,
    U: WaveData,
{
    let channel_: i32 = T::CHANNEL;

    let riff_chunk_ID: [i8; 4] = ['R' as i8, 'I' as i8, 'F' as i8, 'F' as i8];
    let riff_chunk_size: i32 = U::MYSTERIOUS + pcm.get_length() as i32 * U::BYTE_NUM * channel_;
    let file_format_type: [i8; 4] = ['W' as i8, 'A' as i8, 'V' as i8, 'E' as i8];
    let fmt_chunk_ID: [i8; 4] = ['f' as i8, 'm' as i8, 't' as i8, ' ' as i8];
    let fmt_chunk_size: i32 = U::CHUNK_SIZE;
    let wave_format_type: i16 = U::WAVE_FORMAT_TYPE;
    let channel: i16 = channel_ as i16;
    let samples_per_sec: i32 = pcm.get_fs() as i32; /* 標本化周波数 */
    let bytes_per_sec: i32 = pcm.get_fs() as i32 * pcm.get_bits() / 8 * channel_;
    let block_size: i16 = (pcm.get_bits() / 8) as i16 * channel;
    let bits_per_sample: i16 = pcm.get_bits() as i16; /* 量子化精度 */
    let data_chunk_ID: [i8; 4] = ['d' as i8, 'a' as i8, 't' as i8, 'a' as i8];
    let data_chunk_size: i32 = pcm.get_length() as i32 * U::BYTE_NUM * channel_;

    let mut fp = File::create(path).expect("file cannot be created");
    write_i8x4(&mut fp, riff_chunk_ID);
    fp.write_i32::<LittleEndian>(riff_chunk_size).unwrap();
    write_i8x4(&mut fp, file_format_type);
    write_i8x4(&mut fp, fmt_chunk_ID);
    fp.write_i32::<LittleEndian>(fmt_chunk_size).unwrap();
    fp.write_i16::<LittleEndian>(wave_format_type).unwrap();
    fp.write_i16::<LittleEndian>(channel).unwrap();
    fp.write_i32::<LittleEndian>(samples_per_sec).unwrap();
    fp.write_i32::<LittleEndian>(bytes_per_sec).unwrap();
    fp.write_i16::<LittleEndian>(block_size).unwrap();
    fp.write_i16::<LittleEndian>(bits_per_sample).unwrap();
    if U::CHUNK_SIZE > 16 {
        let extra_size: i16 = 0;
        let fact_chunk_ID: [i8; 4] = ['f' as i8, 'a' as i8, 'c' as i8, 't' as i8];
        let fact_chunk_size: i32 = 4;
        let sample_length: i32 = pcm.get_length() as i32;

        fp.write_i16::<LittleEndian>(extra_size).unwrap();
        write_i8x4(&mut fp, fact_chunk_ID);
        fp.write_i32::<LittleEndian>(fact_chunk_size).unwrap();
        fp.write_i32::<LittleEndian>(sample_length).unwrap();
    }

    write_i8x4(&mut fp, data_chunk_ID);
    fp.write_i32::<LittleEndian>(data_chunk_size).unwrap();
    return fp;
}
/*


  int16_t extra_size;
  int8_t  fact_chunk_ID[4];
  int32_t fact_chunk_size;
  int32_t sample_length;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  
 
  
  static int16_t 
  
  

  for (n = 0; n < pcm->length; n++)
  {
 
  }
  
  if ((pcm->length % 2) == 1) /* 圧縮データの長さが奇数のとき */
  {
    c = 0;
    fwrite(&c, 1, 1, fp); /* 0パディング */
  }
  
  fclose(fp);
}

*/

impl WaveData for PCMU {
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
    const BYTE_NUM: i32;
    const MYSTERIOUS: i32;
    const CHUNK_SIZE: i32;
    const WAVE_FORMAT_TYPE: i16;
}

impl WaveData for u8 {
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

pub trait Pcm {
    fn get_fs(&self) -> usize;
    fn get_bits(&self) -> i32;
    fn get_length(&self) -> usize;
    const CHANNEL: i32;
}

impl Pcm for MonoPcm {
    fn get_fs(&self) -> usize {
        self.fs
    }
    fn get_bits(&self) -> i32 {
        self.bits
    }
    fn get_length(&self) -> usize {
        self.length
    }
    const CHANNEL: i32 = 1;
}

impl Pcm for StereoPcm {
    fn get_fs(&self) -> usize {
        self.fs
    }
    fn get_bits(&self) -> i32 {
        self.bits
    }
    fn get_length(&self) -> usize {
        self.length
    }
    const CHANNEL: i32 = 2;
}
