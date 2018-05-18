# sound_programming_rust
[「サウンドプログラミング入門――音響合成の基本とC言語による実装」](https://www.amazon.co.jp/dp/4774155225)のサンプルプログラムRustへの移植

## 内容
[著者様が配布されているサンプルプログラム](http://floor13.sakura.ne.jp/book06/book06.html)をRustに勝手に移植してみる。ライブラリとしての使用が意図されていると思われる `wave.h` (67回登場。テスト済、未移植), `sinc.h` (12回登場。テスト済、移植済), `iir_filter.h` (10回登場。テスト済、移植済), `adsr.h`(8回登場。テスト済、移植済), `window_function.h` (7回登場。テスト済、移植済), `fft.h` (4回登場。テスト済、移植済), `fir_filter.h` (4回登場。テスト済、移植済) を移植するとともに、サンプルプログラムも全て移植し、（乱数が登場するもの以外に関しては）元ファイルとmd5が一致することを目標とする。  
12章については、Windowsに依存し面倒であるため移植を行わない。  
なお、作者様に連絡が取れていないため、何かあったら公開を停止する。

元のCのソースは上記の[配布元](http://floor13.sakura.ne.jp/book06/book06.html)を参照のこと。

## 現状

* 2018/05/17 8章までのサンプルの移植が完了
* 2018/05/16 7章までのサンプルの移植が完了
* 2018/05/16 fir_filter.cの完全移植が完了
* 2018/05/16 iir_filter.cの完全移植が完了
* 2018/05/16 6章までのサンプルの移植が完了
* 2018/05/16 5章までのサンプルの移植が完了
* 2018/05/16 fft.cの完全移植が完了
* 2018/05/16 4章までのサンプルの移植が完了
* 2018/05/16 adsr.cの完全移植が完了
* 2018/05/15 3章までのサンプルの移植が完了
* 2018/05/15 window_function.cの完全移植が完了
* 2018/05/15 sinc.cの完全移植が完了

