extern crate cc;

fn main() {
    cc::Build::new()
        .file("src/adsr.c")
        .compile("libadsr.a");
    cc::Build::new()
        .file("src/fir_filter.c")
        .compile("libfir_filter.a");
    cc::Build::new()
        .file("src/iir_filter.c")
        .compile("libiir_filter.a");      
    cc::Build::new()
        .file("src/wave.c")
        .compile("libwave.a");    
}