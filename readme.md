# Convcodes

Small rust library to provide
[convolutional codes](https://en.wikipedia.org/wiki/Convolutional_code). Has
both an encoder and a decoder using the
[Viterbi algorithm](https://www2.isye.gatech.edu/~yxie77/ece587/viterbi_algorithm.pdf).

Currently missing support for puncturing, but the decoder should work fine if
you simply feed it 0.0 for any bit positions that are punctured.

To generate some BER curves for a small constraint K=3 or the
[CCSDS K=7](https://ccsds.org/Pubs/131x0b5.pdf) rate 1/2 codes, try the example
binaries:

```
$ cargo run --release --bin constraint_3
$ cargo run --release --bin ccsds_rate_1_2_ber_curve
```

![BER Curves](screenshots/BER_Curves.png?raw=true)
