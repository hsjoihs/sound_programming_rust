extern crate cc;

fn main() {
    cc::Build::new()
        .file("src/fir_filter.c")
        .compile("libfir_filter.a");   
    cc::Build::new()
        .file("src/wave.c")
        .compile("libwave.a");    
}