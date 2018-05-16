extern crate cc;

fn main() {
    cc::Build::new()
        .file("src/wave.c")
        .compile("libwave.a");    
}