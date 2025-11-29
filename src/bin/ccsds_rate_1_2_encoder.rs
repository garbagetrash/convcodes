use std::io::{BufRead, Write};

use convcodes::*;


fn main() {
    let polynomials = [[1, 1, 1, 1, 0, 0, 1], [1, 0, 1, 1, 0, 1, 1]];
    let mut cc = ConvolutionalEncoder::new(polynomials);
    let input = std::io::stdin();
    let mut output = std::io::stdout();

    // Read input bits from stdin, push through encoder, write to stdout.
    loop {
        let mut input_handle = input.lock();
        let buf = input_handle.fill_buf().expect("failed to fill buffer from stdin");
        let y = cc.push_byte_slice(&buf);
        output.write_all(&y).expect("failed to write_all to stdout");
    }
}
