use awgn::*;
use convcodes::*;

fn main() {
    let polynomials = [[1, 0, 1], [1, 1, 1]];
    let gamma = 5 * (polynomials[0].len() - 1);
    let ec = 1.0;
    let ecn0_dbs = vec![-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0];
    let mut bers = vec![0.0; ecn0_dbs.len()];

    let seed = 0;
    let mut rng = Generator::new(seed, 128);

    for (ii, ecn0_db) in ecn0_dbs.iter().enumerate() {
        let mut total_bit_errors: usize = 0;
        let mut total_bits: usize = 0;
        let ecn0 = 10_f32.powf(ecn0_db / 10.0);
        let sigma = (0.5 * ec / ecn0).sqrt();
        let mut cc = ConvolutionalEncoder::new(polynomials);
        let mut v = ViterbiDecoder::new(polynomials, gamma);
        let mut input_history = vec![0; gamma];
        let mut t_idx = 0;
        let mut cntr = 0;
        loop {
            /*
            if cntr > 0 && cntr % 1000 == 0 {
                print!(
                    "{}/{} = {}",
                    total_bit_errors,
                    total_bits,
                    total_bit_errors as f32 / total_bits as f32
                );
            }
            */
            let x = if rng.rand_f32() > 0.5 { 1 } else { 0 };
            input_history[t_idx] = x;
            let y = cc.push(x);

            // Received bit with noise
            let mut r = [0.0, 0.0];
            r[0] = 2.0 * y[0] as f32 - 1.0 + sigma * rng.randn_f64() as f32;
            r[1] = 2.0 * y[1] as f32 - 1.0 + sigma * rng.randn_f64() as f32;

            if let Some(xest) = v.push(r) {
                if input_history[(t_idx + 1) % gamma] != xest {
                    total_bit_errors += 1;
                }
                total_bits += 1;
            }

            if total_bit_errors > 10000 || cntr > 1_000_000 {
                break;
            }

            cntr += 1;
            t_idx += 1;
            t_idx %= gamma;
        }

        let ber = total_bit_errors as f32 / total_bits as f32;
        bers[ii] = ber;
        //println!("Ec/N0: {} dB, Pb: {}", ecn0_db, ber);
        println!("[{}, {}],", ecn0_db, ber);
    }
}
