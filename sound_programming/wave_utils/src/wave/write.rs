extern crate byteorder;
use self::byteorder::{LittleEndian, WriteBytesExt};
use MonoPcm;
use StereoPcm;
use std::fs::File;
use wave::WaveData;

pub fn write_i8x4<T>(mut fp: T, arr: [i8; 4])
where
    T: byteorder::WriteBytesExt,
{
    for item in arr.iter() {
        fp.write_i8(*item).unwrap();
    }
}

pub fn write_u8x4<T>(mut fp: T, arr: [u8; 4])
where
    T: byteorder::WriteBytesExt,
{
    for item in arr.iter() {
        fp.write_u8(*item).unwrap();
    }
}

#[allow(non_snake_case)]
pub fn wave_write_header_partial(
    path: &str,
    umysterious: i32,
    bar: i32,
    uchunksize: i32,
    wave_format_type: i16,
    tchannel: i32,
    fs: usize,
    bytes_per_sec: i32,
    block_size: i16,
    bits_per_sample: i16,
) -> File {
    let mut fp = File::create(path).expect("file cannot be created");
    let riff_chunk_ID: [i8; 4] = ['R' as i8, 'I' as i8, 'F' as i8, 'F' as i8];
    let riff_chunk_size: i32 = umysterious + bar;
    let file_format_type: [i8; 4] = ['W' as i8, 'A' as i8, 'V' as i8, 'E' as i8];
    let fmt_chunk_ID: [i8; 4] = ['f' as i8, 'm' as i8, 't' as i8, ' ' as i8];
    let fmt_chunk_size: i32 = uchunksize;
    let wave_format_type: i16 = wave_format_type;
    let channel: i16 = tchannel as i16;
    let samples_per_sec: i32 = fs as i32; /* 標本化周波数 */

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

    fp
}

#[allow(non_snake_case)]
pub fn wave_write_header<T, U>(path: &str, pcm: &T) -> File
where
    T: Pcm,
    U: WaveData,
{
    let channel_: i32 = T::CHANNEL;
    let bits = U::BYTE_NUM * 8;

    let bits_per_sample: i16 = bits as i16; /* 量子化精度 */

    let mut fp = wave_write_header_partial(
        path,
        U::MYSTERIOUS,
        U::BYTE_NUM * channel_ * pcm.get_length() as i32,
        U::CHUNK_SIZE,
        U::WAVE_FORMAT_TYPE,
        T::CHANNEL,
        pcm.get_fs(),
        pcm.get_fs() as i32 * bits / 8 * channel_,
        (bits / 8) as i16 * channel_ as i16,
        bits_per_sample,
    );

    if U::CHUNK_SIZE == 18 {
        let extra_size: i16 = 0;
        let fact_chunk_ID: [i8; 4] = ['f' as i8, 'a' as i8, 'c' as i8, 't' as i8];
        let fact_chunk_size: i32 = 4;
        let sample_length: i32 = pcm.get_length() as i32;

        fp.write_i16::<LittleEndian>(extra_size).unwrap();
        write_i8x4(&mut fp, fact_chunk_ID);
        fp.write_i32::<LittleEndian>(fact_chunk_size).unwrap();
        fp.write_i32::<LittleEndian>(sample_length).unwrap();
    }
    let data_chunk_ID: [i8; 4] = ['d' as i8, 'a' as i8, 't' as i8, 'a' as i8];
    let data_chunk_size: i32 = pcm.get_length() as i32 * U::BYTE_NUM * channel_;
    write_i8x4(&mut fp, data_chunk_ID);
    fp.write_i32::<LittleEndian>(data_chunk_size).unwrap();
    return fp;
}

pub trait Pcm {
    fn get_fs(&self) -> usize;
    fn get_length(&self) -> usize;
    const CHANNEL: i32;
}

impl Pcm for MonoPcm {
    fn get_fs(&self) -> usize {
        self.fs
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
    fn get_length(&self) -> usize {
        self.length
    }
    const CHANNEL: i32 = 2;
}
