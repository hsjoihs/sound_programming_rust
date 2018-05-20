extern crate byteorder;
use self::byteorder::{LittleEndian, ReadBytesExt};
use std::fs::File;

fn read_i8x4<T>(mut fp: T) -> [i8; 4]
where
    T: byteorder::ReadBytesExt,
{
    let mut arr = [0; 4];
    for item in arr.iter_mut() {
        *item = fp.read_i8().unwrap();
    }
    return arr;
}

#[allow(non_snake_case)]
pub fn read_header(file_name: &str) -> (File, i32, i16, i32) {
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

    return (fp, samples_per_sec, bits_per_sample, data_chunk_size);
}
