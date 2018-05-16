/* code copied (and possibly modified) from http://floor13.sakura.ne.jp/book06/book06.html */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


typedef struct
{
  int fs; /* 標本化周波数 */
  int bits; /* 量子化精度 */
  int length; /* 音データの長さ */
  double *s; /* 音データ */
} MONO_PCM;

typedef struct
{
  int fs; /* 標本化周波数 */
  int bits; /* 量子化精度 */
  int length; /* 音データの長さ */
  double *sL; /* 音データ（Lチャンネル） */
  double *sR; /* 音データ（Rチャンネル） */
} STEREO_PCM;


#define READ_ARR(name, bit, length) int ## bit ## _t name[length]; fread((name), (bit)/8, (length), fp)
#define READ(name, bit) int ## bit ## _t name; fread(&(name), (bit)/8, 1, fp)

void wave_read_8bit_mono(MONO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  fp = fopen(file_name, "rb");

  int8_t riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  uint8_t data;
  int n;
  
  
  fread(riff_chunk_ID, 1, 4, fp);
  fread(&riff_chunk_size, 4, 1, fp);
  fread(file_format_type, 1, 4, fp);
  fread(fmt_chunk_ID, 1, 4, fp);
  fread(&fmt_chunk_size, 4, 1, fp);
  fread(&wave_format_type, 2, 1, fp);
  fread(&channel, 2, 1, fp);
  fread(&samples_per_sec, 4, 1, fp);
  fread(&bytes_per_sec, 4, 1, fp);
  fread(&block_size, 2, 1, fp);
  fread(&bits_per_sample, 2, 1, fp);
  fread(data_chunk_ID, 1, 4, fp);
  fread(&data_chunk_size, 4, 1, fp);
  
  pcm->fs = samples_per_sec; /* 標本化周波数 */
  pcm->bits = bits_per_sample; /* 量子化精度 */
  pcm->length = data_chunk_size; /* 音データの長さ */
  pcm->s = calloc(pcm->length, sizeof(double)); /* メモリの確保 */
  
  for (n = 0; n < pcm->length; n++)
  {
    fread(&data, 1, 1, fp); /* 音データの読み取り */
    pcm->s[n] = ((double)data - 128.0) / 128.0; /* 音データを-1以上1未満の範囲に正規化する */
  }
  
  fclose(fp);
}

void wave_write_8bit_mono(MONO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  double s;
  uint8_t data;
  int n;
  
  riff_chunk_ID[0] = 'R';
  riff_chunk_ID[1] = 'I';
  riff_chunk_ID[2] = 'F';
  riff_chunk_ID[3] = 'F';
  riff_chunk_size = 36 + pcm->length;
  file_format_type[0] = 'W';
  file_format_type[1] = 'A';
  file_format_type[2] = 'V';
  file_format_type[3] = 'E';
  
  fmt_chunk_ID[0] = 'f';
  fmt_chunk_ID[1] = 'm';
  fmt_chunk_ID[2] = 't';
  fmt_chunk_ID[3] = ' ';
  fmt_chunk_size = 16;
  wave_format_type = 1;
  channel = 1;
  samples_per_sec = pcm->fs; /* 標本化周波数 */
  bytes_per_sec = pcm->fs * pcm->bits / 8;
  block_size = pcm->bits / 8;
  bits_per_sample = pcm->bits; /* 量子化精度 */
  
  data_chunk_ID[0] = 'd';
  data_chunk_ID[1] = 'a';
  data_chunk_ID[2] = 't';
  data_chunk_ID[3] = 'a';
  data_chunk_size = pcm->length;
  
  fp = fopen(file_name, "wb");
  
  fwrite(riff_chunk_ID, 1, 4, fp);
  fwrite(&riff_chunk_size, 4, 1, fp);
  fwrite(file_format_type, 1, 4, fp);
  fwrite(fmt_chunk_ID, 1, 4, fp);
  fwrite(&fmt_chunk_size, 4, 1, fp);
  fwrite(&wave_format_type, 2, 1, fp);
  fwrite(&channel, 2, 1, fp);
  fwrite(&samples_per_sec, 4, 1, fp);
  fwrite(&bytes_per_sec, 4, 1, fp);
  fwrite(&block_size, 2, 1, fp);
  fwrite(&bits_per_sample, 2, 1, fp);
  fwrite(data_chunk_ID, 1, 4, fp);
  fwrite(&data_chunk_size, 4, 1, fp);
  
  for (n = 0; n < pcm->length; n++)
  {
    s = (pcm->s[n] + 1.0) / 2.0 * 256.0;
    
    if (s > 255.0)
    {
      s = 255.0; /* クリッピング */
    }
    else if (s < 0.0)
    {
      s = 0.0; /* クリッピング */
    }
    
    data = (unsigned char)((int)(s + 0.5)); /* 四捨五入 */
    fwrite(&data, 1, 1, fp); /* 音データの書き出し */
  }
  
  if ((pcm->length % 2) == 1) /* 音データの長さが奇数のとき */
  {
    data = 0;
    fwrite(&data, 1, 1, fp); /* 0パディング */
  }
  
  fclose(fp);
}

void wave_read_8bit_stereo(STEREO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  fp = fopen(file_name, "rb");
  READ_ARR(riff_chunk_ID, 8, 4);
  READ(riff_chunk_size, 32);
  READ_ARR(file_format_type, 8, 4);
  READ_ARR(fmt_chunk_ID, 8, 4);
  READ(fmt_chunk_size, 32);
  READ(wave_format_type,16);
  READ(channel,16);
  READ(samples_per_sec,32);
  READ(bytes_per_sec,32);
  READ(block_size,16);
  READ(bits_per_sample,16);
  READ_ARR(data_chunk_ID, 8, 4);
  READ(data_chunk_size, 32);

  
  
  
  pcm->fs = samples_per_sec; /* 標本化周波数 */
  pcm->bits = bits_per_sample; /* 量子化精度 */
  pcm->length = data_chunk_size / 2; /* 音データの長さ */
  pcm->sL = calloc(pcm->length, sizeof(double)); /* メモリの確保 */
  pcm->sR = calloc(pcm->length, sizeof(double)); /* メモリの確保 */
  
  for (int n = 0; n < pcm->length; n++)
  {
    uint8_t  data;
    fread(&data, 1, 1, fp); /* 音データ（Lチャンネル）の読み取り */
    pcm->sL[n] = ((double)data - 128.0) / 128.0; /* 音データを-1以上1未満の範囲に正規化する */
    
    fread(&data, 1, 1, fp); /* 音データ（Rチャンネル）の読み取り */
    pcm->sR[n] = ((double)data - 128.0) / 128.0; /* 音データを-1以上1未満の範囲に正規化する */
  }
  
  fclose(fp);
}

void wave_write_8bit_stereo(STEREO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  double sL;
  double sR;
  uint8_t data;
  int n;
  
  riff_chunk_ID[0] = 'R';
  riff_chunk_ID[1] = 'I';
  riff_chunk_ID[2] = 'F';
  riff_chunk_ID[3] = 'F';
  riff_chunk_size = 36 + pcm->length * 2;
  file_format_type[0] = 'W';
  file_format_type[1] = 'A';
  file_format_type[2] = 'V';
  file_format_type[3] = 'E';
  
  fmt_chunk_ID[0] = 'f';
  fmt_chunk_ID[1] = 'm';
  fmt_chunk_ID[2] = 't';
  fmt_chunk_ID[3] = ' ';
  fmt_chunk_size = 16;
  wave_format_type = 1;
  channel = 2;
  samples_per_sec = pcm->fs; /* 標本化周波数 */
  bytes_per_sec = pcm->fs * pcm->bits / 8 * 2;
  block_size = pcm->bits / 8 * 2;
  bits_per_sample = pcm->bits; /* 量子化精度 */
  
  data_chunk_ID[0] = 'd';
  data_chunk_ID[1] = 'a';
  data_chunk_ID[2] = 't';
  data_chunk_ID[3] = 'a';
  data_chunk_size = pcm->length * 2;
  
  fp = fopen(file_name, "wb");
  
  fwrite(riff_chunk_ID, 1, 4, fp);
  fwrite(&riff_chunk_size, 4, 1, fp);
  fwrite(file_format_type, 1, 4, fp);
  fwrite(fmt_chunk_ID, 1, 4, fp);
  fwrite(&fmt_chunk_size, 4, 1, fp);
  fwrite(&wave_format_type, 2, 1, fp);
  fwrite(&channel, 2, 1, fp);
  fwrite(&samples_per_sec, 4, 1, fp);
  fwrite(&bytes_per_sec, 4, 1, fp);
  fwrite(&block_size, 2, 1, fp);
  fwrite(&bits_per_sample, 2, 1, fp);
  fwrite(data_chunk_ID, 1, 4, fp);
  fwrite(&data_chunk_size, 4, 1, fp);
  
  for (n = 0; n < pcm->length; n++)
  {
    sL = (pcm->sL[n] + 1.0) / 2.0 * 256.0;
    
    if (sL > 255.0)
    {
      sL = 255.0; /* クリッピング */
    }
    else if (sL < 0.0)
    {
      sL = 0.0; /* クリッピング */
    }
    
    data = (uint8_t)((int)(sL + 0.5)); /* 四捨五入 */
    fwrite(&data, 1, 1, fp); /* 音データ（Lチャンネル）の書き出し */
    
    sR = (pcm->sR[n] + 1.0) / 2.0 * 256.0;
    
    if (sR > 255.0)
    {
      sR = 255.0; /* クリッピング */
    }
    else if (sR < 0.0)
    {
      sR = 0.0; /* クリッピング */
    }
    
    data = (uint8_t)((int)(sR + 0.5)); /* 四捨五入 */
    fwrite(&data, 1, 1, fp); /* 音データ（Rチャンネル）の書き出し */
  }
  
  fclose(fp);
}

void wave_read_16bit_mono(MONO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  fp = fopen(file_name, "rb");
  READ_ARR(riff_chunk_ID, 8, 4);
  READ(riff_chunk_size, 32);
  READ_ARR(file_format_type, 8, 4);
  READ_ARR(fmt_chunk_ID, 8, 4);
  READ(fmt_chunk_size, 32);
  READ(wave_format_type, 16);
  READ(channel, 16);
  READ(samples_per_sec, 32);
  READ(bytes_per_sec, 32);
  READ(block_size, 16);
  READ(bits_per_sample, 16);
  READ_ARR(data_chunk_ID, 8, 4);
  READ(data_chunk_size, 32);
  
  pcm->fs = samples_per_sec; /* 標本化周波数 */
  pcm->bits = bits_per_sample; /* 量子化精度 */
  pcm->length = data_chunk_size / 2; /* 音データの長さ */
  pcm->s = calloc(pcm->length, sizeof(double)); /* メモリの確保 */
  
  for (int n = 0; n < pcm->length; n++)
  {
    READ(data, 16); /* 音データの読み取り */
    pcm->s[n] = (double)data / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
  }
  
  fclose(fp);
}

void wave_write_16bit_mono(MONO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  double s;
  int16_t data;
  int n;
  
  riff_chunk_ID[0] = 'R';
  riff_chunk_ID[1] = 'I';
  riff_chunk_ID[2] = 'F';
  riff_chunk_ID[3] = 'F';
  riff_chunk_size = 36 + pcm->length * 2;
  file_format_type[0] = 'W';
  file_format_type[1] = 'A';
  file_format_type[2] = 'V';
  file_format_type[3] = 'E';
  
  fmt_chunk_ID[0] = 'f';
  fmt_chunk_ID[1] = 'm';
  fmt_chunk_ID[2] = 't';
  fmt_chunk_ID[3] = ' ';
  fmt_chunk_size = 16;
  wave_format_type = 1;
  channel = 1;
  samples_per_sec = pcm->fs; /* 標本化周波数 */
  bytes_per_sec = pcm->fs * pcm->bits / 8;
  block_size = pcm->bits / 8;
  bits_per_sample = pcm->bits; /* 量子化精度 */
  
  data_chunk_ID[0] = 'd';
  data_chunk_ID[1] = 'a';
  data_chunk_ID[2] = 't';
  data_chunk_ID[3] = 'a';
  data_chunk_size = pcm->length * 2;
  
  fp = fopen(file_name, "wb");
  
  fwrite(riff_chunk_ID, 1, 4, fp);
  fwrite(&riff_chunk_size, 4, 1, fp);
  fwrite(file_format_type, 1, 4, fp);
  fwrite(fmt_chunk_ID, 1, 4, fp);
  fwrite(&fmt_chunk_size, 4, 1, fp);
  fwrite(&wave_format_type, 2, 1, fp);
  fwrite(&channel, 2, 1, fp);
  fwrite(&samples_per_sec, 4, 1, fp);
  fwrite(&bytes_per_sec, 4, 1, fp);
  fwrite(&block_size, 2, 1, fp);
  fwrite(&bits_per_sample, 2, 1, fp);
  fwrite(data_chunk_ID, 1, 4, fp);
  fwrite(&data_chunk_size, 4, 1, fp);
  
  for (n = 0; n < pcm->length; n++)
  {
    s = (pcm->s[n] + 1.0) / 2.0 * 65536.0;
    
    if (s > 65535.0)
    {
      s = 65535.0; /* クリッピング */
    }
    else if (s < 0.0)
    {
      s = 0.0; /* クリッピング */
    }
    
    data = (int16_t)((int)(s + 0.5) - 32768); /* 四捨五入とオフセットの調節 */
    fwrite(&data, 2, 1, fp); /* 音データの書き出し */
  }
  
  fclose(fp);
}

void wave_read_16bit_stereo(STEREO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  int16_t data;
  int n;
  
  fp = fopen(file_name, "rb");
  
  fread(riff_chunk_ID, 1, 4, fp);
  fread(&riff_chunk_size, 4, 1, fp);
  fread(file_format_type, 1, 4, fp);
  fread(fmt_chunk_ID, 1, 4, fp);
  fread(&fmt_chunk_size, 4, 1, fp);
  fread(&wave_format_type, 2, 1, fp);
  fread(&channel, 2, 1, fp);
  fread(&samples_per_sec, 4, 1, fp);
  fread(&bytes_per_sec, 4, 1, fp);
  fread(&block_size, 2, 1, fp);
  fread(&bits_per_sample, 2, 1, fp);
  fread(data_chunk_ID, 1, 4, fp);
  fread(&data_chunk_size, 4, 1, fp);
  
  pcm->fs = samples_per_sec; /* 標本化周波数 */
  pcm->bits = bits_per_sample; /* 量子化精度 */
  pcm->length = data_chunk_size / 4; /* 音データの長さ */
  pcm->sL = calloc(pcm->length, sizeof(double)); /* メモリの確保 */
  pcm->sR = calloc(pcm->length, sizeof(double)); /* メモリの確保 */
  
  for (n = 0; n < pcm->length; n++)
  {
    fread(&data, 2, 1, fp); /* 音データ（Lチャンネル）の読み取り */
    pcm->sL[n] = (double)data / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
    
    fread(&data, 2, 1, fp); /* 音データ（Rチャンネル）の読み取り */
    pcm->sR[n] = (double)data / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
  }
  
  fclose(fp);
}

void wave_write_16bit_stereo(STEREO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  double sL;
  double sR;
  int16_t data;
  int n;
  
  riff_chunk_ID[0] = 'R';
  riff_chunk_ID[1] = 'I';
  riff_chunk_ID[2] = 'F';
  riff_chunk_ID[3] = 'F';
  riff_chunk_size = 36 + pcm->length * 4;
  file_format_type[0] = 'W';
  file_format_type[1] = 'A';
  file_format_type[2] = 'V';
  file_format_type[3] = 'E';
  
  fmt_chunk_ID[0] = 'f';
  fmt_chunk_ID[1] = 'm';
  fmt_chunk_ID[2] = 't';
  fmt_chunk_ID[3] = ' ';
  fmt_chunk_size = 16;
  wave_format_type = 1;
  channel = 2;
  samples_per_sec = pcm->fs; /* 標本化周波数 */
  bytes_per_sec = pcm->fs * pcm->bits / 8 * 2;
  block_size = pcm->bits / 8 * 2;
  bits_per_sample = pcm->bits; /* 量子化精度 */
  
  data_chunk_ID[0] = 'd';
  data_chunk_ID[1] = 'a';
  data_chunk_ID[2] = 't';
  data_chunk_ID[3] = 'a';
  data_chunk_size = pcm->length * 4;
  
  fp = fopen(file_name, "wb");
  
  fwrite(riff_chunk_ID, 1, 4, fp);
  fwrite(&riff_chunk_size, 4, 1, fp);
  fwrite(file_format_type, 1, 4, fp);
  fwrite(fmt_chunk_ID, 1, 4, fp);
  fwrite(&fmt_chunk_size, 4, 1, fp);
  fwrite(&wave_format_type, 2, 1, fp);
  fwrite(&channel, 2, 1, fp);
  fwrite(&samples_per_sec, 4, 1, fp);
  fwrite(&bytes_per_sec, 4, 1, fp);
  fwrite(&block_size, 2, 1, fp);
  fwrite(&bits_per_sample, 2, 1, fp);
  fwrite(data_chunk_ID, 1, 4, fp);
  fwrite(&data_chunk_size, 4, 1, fp);
  
  for (n = 0; n < pcm->length; n++)
  {
    sL = (pcm->sL[n] + 1.0) / 2.0 * 65536.0;
    
    if (sL > 65535.0)
    {
      sL = 65535.0; /* クリッピング */
    }
    else if (sL < 0.0)
    {
      sL = 0.0; /* クリッピング */
    }
    
    data = (int16_t)((int)(sL + 0.5) - 32768); /* 四捨五入とオフセットの調節 */
    fwrite(&data, 2, 1, fp); /* 音データ（Lチャンネル）の書き出し */
    
    sR = (pcm->sR[n] + 1.0) / 2.0 * 65536.0;
    
    if (sR > 65535.0)
    {
      sR = 65535.0; /* クリッピング */
    }
    else if (sR < 0.0)
    {
      sR = 0.0; /* クリッピング */
    }
    
    data = (int16_t)((int)(sR + 0.5) - 32768); /* 四捨五入とオフセットの調節 */
    fwrite(&data, 2, 1, fp); /* 音データ（Rチャンネル）の書き出し */
  }
  
  fclose(fp);
}

void wave_read_PCMU_mono(MONO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int16_t extra_size;
  int8_t  fact_chunk_ID[4];
  int32_t fact_chunk_size;
  int32_t sample_length;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  
  int16_t s; /* 16bitの音データ */
  unsigned char c; /* 8bitの圧縮データ */
  unsigned char sign, exponent, mantissa;
  int n, magnitude;
  
  fp = fopen(file_name, "rb");
  
  fread(riff_chunk_ID, 1, 4, fp);
  fread(&riff_chunk_size, 4, 1, fp);
  fread(file_format_type, 1, 4, fp);
  fread(fmt_chunk_ID, 1, 4, fp);
  fread(&fmt_chunk_size, 4, 1, fp);
  fread(&wave_format_type, 2, 1, fp);
  fread(&channel, 2, 1, fp);
  fread(&samples_per_sec, 4, 1, fp);
  fread(&bytes_per_sec, 4, 1, fp);
  fread(&block_size, 2, 1, fp);
  fread(&bits_per_sample, 2, 1, fp);
  fread(&extra_size, 2, 1, fp);
  fread(fact_chunk_ID, 1, 4, fp);
  fread(&fact_chunk_size, 4, 1, fp);
  fread(&sample_length, 4, 1, fp);
  fread(data_chunk_ID, 1, 4, fp);
  fread(&data_chunk_size, 4, 1, fp);
  
  pcm->fs = samples_per_sec; /* 標本化周波数 */
  pcm->bits = 16; /* 量子化精度 */
  pcm->length = sample_length; /* 音データの長さ */
  pcm->s = calloc(pcm->length, sizeof(double)); /* メモリの確保 */
  
  for (n = 0; n < pcm->length; n++)
  {
    fread(&c, 1, 1, fp); /* 圧縮データの読み取り */
    
    c = ~c;
    
    sign = c & 0x80;
    exponent = (c >> 4) & 0x07;
    mantissa = c & 0x0F;
    
    magnitude = ((((int)mantissa << 3) + 0x84) << exponent) - 0x84;
    
    if (sign == 0x80)
    {
      s = -(int16_t)magnitude;
    }
    else
    {
      s = (int16_t)magnitude;
    }
    
    pcm->s[n] = (double)s / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
  }
  
  fclose(fp);
}

void wave_write_PCMU_mono(MONO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int16_t extra_size;
  int8_t  fact_chunk_ID[4];
  int32_t fact_chunk_size;
  int32_t sample_length;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  
  double x;
  int16_t s; /* 16bitの音データ */
  unsigned char c; /* 8bitの圧縮データ */
  unsigned char sign, exponent, mantissa;
  int n, magnitude;
  
  static int16_t level[8] =
  {
    0x00FF, 0x01FF, 0x03FF, 0x07FF, 0x0FFF, 0x1FFF, 0x3FFF, 0x7FFF
  };
  
  riff_chunk_ID[0] = 'R';
  riff_chunk_ID[1] = 'I';
  riff_chunk_ID[2] = 'F';
  riff_chunk_ID[3] = 'F';
  riff_chunk_size = 50 + pcm->length;
  file_format_type[0] = 'W';
  file_format_type[1] = 'A';
  file_format_type[2] = 'V';
  file_format_type[3] = 'E';
  
  fmt_chunk_ID[0] = 'f';
  fmt_chunk_ID[1] = 'm';
  fmt_chunk_ID[2] = 't';
  fmt_chunk_ID[3] = ' ';
  fmt_chunk_size = 18;
  wave_format_type = 7;
  channel = 1;
  samples_per_sec = pcm->fs; /* 標本化周波数 */
  bytes_per_sec = samples_per_sec;
  block_size = 1;
  bits_per_sample = 8; /* 量子化精度 */
  extra_size = 0;
  
  fact_chunk_ID[0] = 'f';
  fact_chunk_ID[1] = 'a';
  fact_chunk_ID[2] = 'c';
  fact_chunk_ID[3] = 't';
  fact_chunk_size = 4;
  sample_length = pcm->length;
  
  data_chunk_ID[0] = 'd';
  data_chunk_ID[1] = 'a';
  data_chunk_ID[2] = 't';
  data_chunk_ID[3] = 'a';
  data_chunk_size = pcm->length;
  
  fp = fopen(file_name, "wb");
  
  fwrite(riff_chunk_ID, 1, 4, fp);
  fwrite(&riff_chunk_size, 4, 1, fp);
  fwrite(file_format_type, 1, 4, fp);
  fwrite(fmt_chunk_ID, 1, 4, fp);
  fwrite(&fmt_chunk_size, 4, 1, fp);
  fwrite(&wave_format_type, 2, 1, fp);
  fwrite(&channel, 2, 1, fp);
  fwrite(&samples_per_sec, 4, 1, fp);
  fwrite(&bytes_per_sec, 4, 1, fp);
  fwrite(&block_size, 2, 1, fp);
  fwrite(&bits_per_sample, 2, 1, fp);
  fwrite(&extra_size, 2, 1, fp);
  fwrite(fact_chunk_ID, 1, 4, fp);
  fwrite(&fact_chunk_size, 4, 1, fp);
  fwrite(&sample_length, 4, 1, fp);
  fwrite(data_chunk_ID, 1, 4, fp);
  fwrite(&data_chunk_size, 4, 1, fp);
  
  for (n = 0; n < pcm->length; n++)
  {
    x = (pcm->s[n] + 1.0) / 2.0 * 65536.0;
    
    if (x > 65535.0)
    {
      x = 65535.0; /* クリッピング */
    }
    else if (x < 0.0)
    {
      x = 0.0; /* クリッピング */
    }
    
    s = (int16_t)((int)(x + 0.5) - 32768); /* 四捨五入とオフセットの調節 */
    
    if (s < 0)
    {
      magnitude = -s;
      sign = 0x80;
    }
    else
    {
      magnitude = s;
      sign = 0x00;
    }
    
    magnitude += 0x84;
    if (magnitude > 0x7FFF)
    {
      magnitude = 0x7FFF;
    }
    
    for (exponent = 0; exponent < 8; exponent++)
    {
      if (magnitude <= level[exponent])
      {
        break;
      }
    }
    
    mantissa = (magnitude >> (exponent + 3)) & 0x0F;
    
    c = ~(sign | (exponent << 4) | mantissa);
    
    fwrite(&c, 1, 1, fp); /* 圧縮データの書き出し */
  }
  
  if ((pcm->length % 2) == 1) /* 圧縮データの長さが奇数のとき */
  {
    c = 0;
    fwrite(&c, 1, 1, fp); /* 0パディング */
  }
  
  fclose(fp);
}

void wave_read_PCMA_mono(MONO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int16_t extra_size;
  int8_t  fact_chunk_ID[4];
  int32_t fact_chunk_size;
  int32_t sample_length;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  
  int16_t s; /* 16bitの音データ */
  unsigned char c; /* 8bitの圧縮データ */
  unsigned char sign, exponent, mantissa;
  int n, magnitude;
  
  fp = fopen(file_name, "rb");
  
  fread(riff_chunk_ID, 1, 4, fp);
  fread(&riff_chunk_size, 4, 1, fp);
  fread(file_format_type, 1, 4, fp);
  fread(fmt_chunk_ID, 1, 4, fp);
  fread(&fmt_chunk_size, 4, 1, fp);
  fread(&wave_format_type, 2, 1, fp);
  fread(&channel, 2, 1, fp);
  fread(&samples_per_sec, 4, 1, fp);
  fread(&bytes_per_sec, 4, 1, fp);
  fread(&block_size, 2, 1, fp);
  fread(&bits_per_sample, 2, 1, fp);
  fread(&extra_size, 2, 1, fp);
  fread(fact_chunk_ID, 1, 4, fp);
  fread(&fact_chunk_size, 4, 1, fp);
  fread(&sample_length, 4, 1, fp);
  fread(data_chunk_ID, 1, 4, fp);
  fread(&data_chunk_size, 4, 1, fp);
  
  pcm->fs = samples_per_sec; /* 標本化周波数 */
  pcm->bits = 16; /* 量子化精度 */
  pcm->length = sample_length; /* 音データの長さ */
  pcm->s = calloc(pcm->length, sizeof(double)); /* メモリの確保 */
  
  for (n = 0; n < pcm->length; n++)
  {
    fread(&c, 1, 1, fp); /* 圧縮データの読み取り */
    
    c ^= 0xD5;
    
    sign = c & 0x80;
    exponent = (c >> 4) & 0x07;
    mantissa = c & 0x0F;
    
    if (exponent == 0)
    {
      magnitude = ((int)mantissa << 4) + 0x0008;
    }
    else
    {
      magnitude = (((int)mantissa << 4) + 0x0108) << (exponent - 1);
    }
    
    if (sign == 0x80)
    {
      s = -(int16_t)magnitude;
    }
    else
    {
      s = (int16_t)magnitude;
    }
    
    pcm->s[n] = (double)s / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
  }
  
  fclose(fp);
}

void wave_write_PCMA_mono(MONO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int16_t extra_size;
  int8_t  fact_chunk_ID[4];
  int32_t fact_chunk_size;
  int32_t sample_length;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  
  double x;
  int16_t s; /* 16bitの音データ */
  unsigned char c; /* 8bitの圧縮データ */
  unsigned char sign, exponent, mantissa;
  int n, magnitude;
  
  static int level[8] =
  {
    0x00FF, 0x01FF, 0x03FF, 0x07FF, 0x0FFF, 0x1FFF, 0x3FFF, 0x7FFF
  };
  
  riff_chunk_ID[0] = 'R';
  riff_chunk_ID[1] = 'I';
  riff_chunk_ID[2] = 'F';
  riff_chunk_ID[3] = 'F';
  riff_chunk_size = 50 + pcm->length;
  file_format_type[0] = 'W';
  file_format_type[1] = 'A';
  file_format_type[2] = 'V';
  file_format_type[3] = 'E';
  
  fmt_chunk_ID[0] = 'f';
  fmt_chunk_ID[1] = 'm';
  fmt_chunk_ID[2] = 't';
  fmt_chunk_ID[3] = ' ';
  fmt_chunk_size = 18;
  wave_format_type = 6;
  channel = 1;
  samples_per_sec = pcm->fs; /* 標本化周波数 */
  bytes_per_sec = samples_per_sec;
  block_size = 1;
  bits_per_sample = 8; /* 量子化精度 */
  extra_size = 0;
  
  fact_chunk_ID[0] = 'f';
  fact_chunk_ID[1] = 'a';
  fact_chunk_ID[2] = 'c';
  fact_chunk_ID[3] = 't';
  fact_chunk_size = 4;
  sample_length = pcm->length;
  
  data_chunk_ID[0] = 'd';
  data_chunk_ID[1] = 'a';
  data_chunk_ID[2] = 't';
  data_chunk_ID[3] = 'a';
  data_chunk_size = pcm->length;
  
  fp = fopen(file_name, "wb");
  
  fwrite(riff_chunk_ID, 1, 4, fp);
  fwrite(&riff_chunk_size, 4, 1, fp);
  fwrite(file_format_type, 1, 4, fp);
  fwrite(fmt_chunk_ID, 1, 4, fp);
  fwrite(&fmt_chunk_size, 4, 1, fp);
  fwrite(&wave_format_type, 2, 1, fp);
  fwrite(&channel, 2, 1, fp);
  fwrite(&samples_per_sec, 4, 1, fp);
  fwrite(&bytes_per_sec, 4, 1, fp);
  fwrite(&block_size, 2, 1, fp);
  fwrite(&bits_per_sample, 2, 1, fp);
  fwrite(&extra_size, 2, 1, fp);
  fwrite(fact_chunk_ID, 1, 4, fp);
  fwrite(&fact_chunk_size, 4, 1, fp);
  fwrite(&sample_length, 4, 1, fp);
  fwrite(data_chunk_ID, 1, 4, fp);
  fwrite(&data_chunk_size, 4, 1, fp);
  
  for (n = 0; n < pcm->length; n++)
  {
    x = (pcm->s[n] + 1.0) / 2.0 * 65536.0;
    
    if (x > 65535.0)
    {
      x = 65535.0; /* クリッピング */
    }
    else if (x < 0.0)
    {
      x = 0.0; /* クリッピング */
    }
    
    s = (int16_t)((int)(x + 0.5) - 32768); /* 四捨五入とオフセットの調節 */
    
    if (s < 0)
    {
      magnitude = -s;
      sign = 0x80;
    }
    else
    {
      magnitude = s;
      sign = 0x00;
    }
    
    if (magnitude > 0x7FFF)
    {
      magnitude = 0x7FFF;
    }
    
    for (exponent = 0; exponent < 8; exponent++)
    {
      if (magnitude <= level[exponent])
      {
        break;
      }
    }
    
    if (exponent == 0)
    {
      mantissa = (magnitude >> 4) & 0x0F;
    }
    else
    {
      mantissa = (magnitude >> (exponent + 3)) & 0x0F;
    }
    
    c = (sign | (exponent << 4) | mantissa) ^ 0xD5;
    
    fwrite(&c, 1, 1, fp); /* 圧縮データの書き出し */
  }
  
  if ((pcm->length % 2) == 1) /* 圧縮データの長さが奇数のとき */
  {
    c = 0;
    fwrite(&c, 1, 1, fp); /* 0パディング */
  }
  
  fclose(fp);
}

void wave_read_IMA_ADPCM_mono(MONO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int16_t extra_size;
  int16_t samples_per_block;
  char fact_chunk_ID[4];
  int32_t fact_chunk_size;
  int32_t sample_length;
  char data_chunk_ID[4];
  int32_t data_chunk_size;
  
  int16_t s; /* 16bitの音データ */
  unsigned char c; /* 4bitの圧縮データ */
  unsigned char header[4];
  unsigned char data;
  int n, sp, dp, offset, block, number_of_block, index, step_size;
  
  static int index_table[16] =
  {
    -1, -1, -1, -1, 2, 4, 6, 8,
    -1, -1, -1, -1, 2, 4, 6, 8
  };
  
  static int step_size_table[89] =
  {
    7, 8, 9, 10, 11, 12, 13, 14,
    16, 17, 19, 21, 23, 25, 28, 31,
    34, 37, 41, 45, 50, 55, 60, 66,
    73, 80, 88, 97, 107, 118, 130, 143,
    157, 173, 190, 209, 230, 253, 279, 307,
    337, 371, 408, 449, 494, 544, 598, 658,
    724, 796, 876, 963, 1060, 1166, 1282, 1411,
    1552, 1707, 1878, 2066, 2272, 2499, 2749, 3024,
    3327, 3660, 4026, 4428, 4871, 5358, 5894, 6484,
    7132, 7845, 8630, 9493, 10442, 11487, 12635, 13899,
    15289, 16818, 18500, 20350, 22385, 24623, 27086, 29794,
    32767
  };
  
  fp = fopen(file_name, "rb");
  
  fread(riff_chunk_ID, 1, 4, fp);
  fread(&riff_chunk_size, 4, 1, fp);
  fread(file_format_type, 1, 4, fp);
  fread(fmt_chunk_ID, 1, 4, fp);
  fread(&fmt_chunk_size, 4, 1, fp);
  fread(&wave_format_type, 2, 1, fp);
  fread(&channel, 2, 1, fp);
  fread(&samples_per_sec, 4, 1, fp);
  fread(&bytes_per_sec, 4, 1, fp);
  fread(&block_size, 2, 1, fp);
  fread(&bits_per_sample, 2, 1, fp);
  fread(&extra_size, 2, 1, fp);
  fread(&samples_per_block, 2, 1, fp);
  fread(fact_chunk_ID, 1, 4, fp);
  fread(&fact_chunk_size, 4, 1, fp);
  fread(&sample_length, 4, 1, fp);
  fread(data_chunk_ID, 1, 4, fp);
  fread(&data_chunk_size, 4, 1, fp);
  
  number_of_block = data_chunk_size / block_size;
  
  pcm->fs = samples_per_sec; /* 標本化周波数 */
  pcm->bits = 16; /* 量子化精度 */
  pcm->length = sample_length; /* 音データの長さ */
  pcm->s = calloc(pcm->length, sizeof(double)); /* メモリの確保 */
  
  for (block = 0; block < number_of_block; block++)
  {
    offset = samples_per_block * block;
    
    for (n = 0; n < samples_per_block; n++)
    {
      if (n == 0)
      {
        fread(header, 1, 4, fp); /* ヘッダの読み取り */
        
        sp = ((int16_t)(int8_t)header[1] << 8) + header[0];
        index = header[2];
        
        s = sp;
      }
      else
      {
        if ((n % 2) == 1)
        {
          fread(&data, 1, 1, fp); /* 圧縮データの読み取り */
          
          c = (uint8_t)(data & 0x0F); /* dataの下位4bit */
        }
        else
        {
          c = (uint8_t)((data >> 4) & 0x0F); /* dataの上位4bit */
        }
        
        step_size = step_size_table[index];
        
        /* 伸張 */
        dp = step_size >> 3;
        if ((c & 0x1) == 0x1)
        {
          dp += (step_size >> 2);
        }
        if ((c & 0x2) == 0x2)
        {
          dp += (step_size >> 1);
        }
        if ((c & 0x4) == 0x4)
        {
          dp += step_size;
        }
        
        if ((c & 0x8) == 0x8)
        {
          sp -= dp;
        }
        else
        {
          sp += dp;
        }
        
        if (sp > 32767)
        {
          sp = 32767;
        }
        else if (sp < -32768)
        {
          sp = -32768;
        }
        
        index += index_table[c];
        if (index < 0)
        {
          index = 0;
        }
        else if (index > 88)
        {
          index = 88;
        }
        
        s = sp;
      }
      
      pcm->s[offset + n] = (double)s / 32768.0; /* 音データを-1以上1未満の範囲に正規化する */
    }
  }
  
  fclose(fp);
}

void wave_write_IMA_ADPCM_mono(MONO_PCM *pcm, const char *file_name)
{
  FILE *fp;
  int8_t  riff_chunk_ID[4];
  int32_t riff_chunk_size;
  int8_t  file_format_type[4];
  int8_t  fmt_chunk_ID[4];
  int32_t fmt_chunk_size;
  int16_t wave_format_type;
  int16_t channel;
  int32_t samples_per_sec;
  int32_t bytes_per_sec;
  int16_t block_size;
  int16_t bits_per_sample;
  int16_t extra_size;
  int16_t samples_per_block;
  int8_t  fact_chunk_ID[4];
  int32_t fact_chunk_size;
  int32_t sample_length;
  int8_t  data_chunk_ID[4];
  int32_t data_chunk_size;
  
  double x;
  int16_t s; /* 16bitの音データ */
  unsigned char c; /* 4bitの圧縮データ */
  unsigned char header[4];
  unsigned char data;
  int n, sp, d, dp, offset, block, number_of_block, index, step_size;
  
  static int index_table[16] =
  {
    -1, -1, -1, -1, 2, 4, 6, 8,
    -1, -1, -1, -1, 2, 4, 6, 8
  };
  
  static int step_size_table[89] =
  {
    7, 8, 9, 10, 11, 12, 13, 14,
    16, 17, 19, 21, 23, 25, 28, 31,
    34, 37, 41, 45, 50, 55, 60, 66,
    73, 80, 88, 97, 107, 118, 130, 143,
    157, 173, 190, 209, 230, 253, 279, 307,
    337, 371, 408, 449, 494, 544, 598, 658,
    724, 796, 876, 963, 1060, 1166, 1282, 1411,
    1552, 1707, 1878, 2066, 2272, 2499, 2749, 3024,
    3327, 3660, 4026, 4428, 4871, 5358, 5894, 6484,
    7132, 7845, 8630, 9493, 10442, 11487, 12635, 13899,
    15289, 16818, 18500, 20350, 22385, 24623, 27086, 29794,
    32767
  };
  
  block_size = 256;
  samples_per_block = (block_size - 4) * 2 + 1;
  number_of_block = (int)(pcm->length / samples_per_block);
  
  riff_chunk_ID[0] = 'R';
  riff_chunk_ID[1] = 'I';
  riff_chunk_ID[2] = 'F';
  riff_chunk_ID[3] = 'F';
  riff_chunk_size = 52 + block_size * number_of_block;
  file_format_type[0] = 'W';
  file_format_type[1] = 'A';
  file_format_type[2] = 'V';
  file_format_type[3] = 'E';
  
  fmt_chunk_ID[0] = 'f';
  fmt_chunk_ID[1] = 'm';
  fmt_chunk_ID[2] = 't';
  fmt_chunk_ID[3] = ' ';
  fmt_chunk_size = 20;
  wave_format_type = 17;
  channel = 1;
  samples_per_sec = pcm->fs; /* 標本化周波数 */
  bytes_per_sec = (int32_t)(block_size * samples_per_sec / samples_per_block);
  bits_per_sample = 4; /* 量子化精度 */
  extra_size = 2;
  
  fact_chunk_ID[0] = 'f';
  fact_chunk_ID[1] = 'a';
  fact_chunk_ID[2] = 'c';
  fact_chunk_ID[3] = 't';
  fact_chunk_size = 4;
  sample_length = samples_per_block * number_of_block + 1;
  
  data_chunk_ID[0] = 'd';
  data_chunk_ID[1] = 'a';
  data_chunk_ID[2] = 't';
  data_chunk_ID[3] = 'a';
  data_chunk_size = block_size * number_of_block;
  
  fp = fopen(file_name, "wb");
  
  fwrite(riff_chunk_ID, 1, 4, fp);
  fwrite(&riff_chunk_size, 4, 1, fp);
  fwrite(file_format_type, 1, 4, fp);
  fwrite(fmt_chunk_ID, 1, 4, fp);
  fwrite(&fmt_chunk_size, 4, 1, fp);
  fwrite(&wave_format_type, 2, 1, fp);
  fwrite(&channel, 2, 1, fp);
  fwrite(&samples_per_sec, 4, 1, fp);
  fwrite(&bytes_per_sec, 4, 1, fp);
  fwrite(&block_size, 2, 1, fp);
  fwrite(&bits_per_sample, 2, 1, fp);
  fwrite(&extra_size, 2, 1, fp);
  fwrite(&samples_per_block, 2, 1, fp);
  fwrite(fact_chunk_ID, 1, 4, fp);
  fwrite(&fact_chunk_size, 4, 1, fp);
  fwrite(&sample_length, 4, 1, fp);
  fwrite(data_chunk_ID, 1, 4, fp);
  fwrite(&data_chunk_size, 4, 1, fp);
  
  for (block = 0; block < number_of_block; block++)
  {
    offset = samples_per_block * block;
    
    for (n = 0; n < samples_per_block; n++)
    {
      x = (pcm->s[offset + n] + 1.0) / 2.0 * 65536.0;
      
      if (x > 65535.0)
      {
        x = 65535.0; /* クリッピング */
      }
      else if (x < 0.0)
      {
        x = 0.0; /* クリッピング */
      }
      
      s = (int16_t)((int)(x + 0.5) - 32768); /* 四捨五入とオフセットの調節 */
      
      if (block == 0 && n == 0)
      {
        index = 0; /* 最初のブロックにおけるindexの初期値を0にする */
      }
      
      if (n == 0)
      {
        header[0] = (uint8_t)(s & 0x00FF); /* sの下位バイト */
        header[1] = (uint8_t)((s >> 8) & 0x00FF); /* sの上位バイト */
        header[2] = (uint8_t)index; /* インデックス */
        header[3] = 0;
        
        fwrite(header, 1, 4, fp); /* ヘッダの書き出し */
        
        sp = s; /* spの初期値をsにする */
      }
      else
      {
        d = s - sp;
        if (d < 0)
        {
          c = 0x8;
          d = -d;
        }
        else
        {
          c = 0x0;
        }
        
        step_size = step_size_table[index];
        
        /* 圧縮 */
        if (d >= step_size)
        {
          c |= 0x4;
          d -= step_size;
        }
        if (d >= (step_size >> 1))
        {
          c |= 0x2;
          d -= (step_size >> 1);
        }
        if (d >= (step_size >> 2))
        {
          c |= 0x1;
        }
        
        /* 伸張 */
        dp = step_size >> 3;
        if ((c & 0x1) == 0x1)
        {
          dp += (step_size >> 2);
        }
        if ((c & 0x2) == 0x2)
        {
          dp += (step_size >> 1);
        }
        if ((c & 0x4) == 0x4)
        {
          dp += step_size;
        }
        
        if ((c & 0x8) == 0x8)
        {
          sp -= dp;
        }
        else
        {
          sp += dp;
        }
        
        if (sp > 32767)
        {
          sp = 32767;
        }
        else if (sp < -32768)
        {
          sp = -32768;
        }
        
        index += index_table[c];
        if (index < 0)
        {
          index = 0;
        }
        else if (index > 88)
        {
          index = 88;
        }
        
        if ((n % 2) == 1)
        {
          data = c & 0xF; /* dataの下位4bit */
        }
        else
        {
          data |= (c & 0xF) << 4; /* dataの上位4bit */
          
          fwrite(&data, 1, 1, fp); /* 圧縮データの書き出し */
        }
      }
    }
  }
  
  fclose(fp);
}
