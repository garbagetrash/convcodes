pub struct ConvolutionalEncoder<const N: usize, const M: usize> {
    polynomials: [[u8; M]; N],
    lfsr: Vec<u8>,
    idx: usize,
}

impl<const N: usize, const M: usize> ConvolutionalEncoder<N, M> {
    pub fn new(polynomials: [[u8; M]; N]) -> Self {
        Self {
            polynomials,
            lfsr: vec![0; 2 * M],
            idx: 0,
        }
    }

    pub fn push(&mut self, x: u8) -> [u8; N] {
        // Push new bit into LFSR
        self.idx = if self.idx == 0 { M - 1 } else { self.idx - 1 };
        self.lfsr[self.idx] = x;
        self.lfsr[self.idx + M] = x;

        // Calculate encoded bits
        let mut output = [0; N];
        for (i, g_row) in self.polynomials.iter().enumerate() {
            for (gg, mm) in g_row.iter().zip(self.lfsr[self.idx..self.idx + M].iter()) {
                output[i] ^= gg & mm;
            }
        }
        output
    }

    pub fn push_block(&mut self, x: &[u8]) -> Vec<u8> {
        x.iter().flat_map(|&xx| self.push(xx)).collect()
    }
}

pub fn int_to_bitvec<const N: usize>(x: usize) -> [u8; N] {
    let mut output = [0; N];
    for i in 0..N {
        output[N - i - 1] = ((x >> i) & 1) as u8;
    }
    output
}

pub fn bitvec_to_int(bitvec: &[u8]) -> usize {
    let n = bitvec.len();
    let mut output = 0;
    for i in 0..n {
        output += (1 << i) * bitvec[n - i - 1] as usize;
    }
    output
}

pub fn build_trellis<const N: usize, const M: usize>(
    polynomials: [[u8; M]; N],
) -> Vec<(usize, usize, u8, [u8; N])> {
    let s = 2_usize.pow(M as u32 - 1);
    let mut output = vec![];
    for ss in 0..s {
        let mut v = int_to_bitvec::<M>(ss);

        // Simulate pushing in `0` to state ss
        v[0] = 0;
        let mut outbits = [0; N];
        for (i, g_row) in polynomials.iter().enumerate() {
            for (gg, mm) in g_row.iter().zip(v.iter()) {
                outbits[i] ^= gg & mm;
            }
        }
        let mut next_state = v[1..].to_vec();
        next_state.rotate_right(1);
        next_state[0] = 0;
        output.push((ss, bitvec_to_int(&next_state), 0, outbits));

        // Simulate pushing in `1` to state ss
        v[0] = 1;
        let mut outbits = [0; N];
        for (i, g_row) in polynomials.iter().enumerate() {
            for (gg, mm) in g_row.iter().zip(v.iter()) {
                outbits[i] ^= gg & mm;
            }
        }
        next_state = v[1..].to_vec();
        next_state.rotate_right(1);
        next_state[0] = 1;
        output.push((ss, bitvec_to_int(&next_state), 1, outbits));
    }
    output
}

fn llr(x: &[f32], a: &[f32]) -> f32 {
    assert_eq!(x.len(), a.len());
    x.iter()
        .zip(a.iter())
        .map(|(xx, aa)| (xx - aa).powi(2))
        .sum::<f32>()
}

pub struct ViterbiDecoder<const N: usize, const M: usize> {
    trellis: Vec<(usize, usize, u8, [u8; N])>,
    gamma: usize,
    s: usize,
    prior_tups: Vec<Vec<(usize, usize, u8, [u8; N])>>,
    metrics: Vec<Vec<f32>>,
    prior_states: Vec<Vec<usize>>,
    t_idx: usize,
    valid: bool,
    valid_really: bool,
}

impl<const N: usize, const M: usize> ViterbiDecoder<N, M> {
    pub fn new(polynomials: [[u8; M]; N], gamma: usize) -> Self {
        let trellis = build_trellis(polynomials);
        let s = trellis.len() / 2;
        let mut prior_tups = vec![];
        for ss in 0..s {
            let mut pts = vec![];
            for ptup in &trellis {
                if ptup.1 == ss {
                    pts.push(*ptup);
                }
            }
            prior_tups.push(pts);
        }
        let mut metrics = vec![vec![f32::INFINITY; 2 * gamma]; s];
        metrics[0][0] = 0.0;
        metrics[0][gamma] = 0.0;
        Self {
            trellis,
            gamma,
            s,
            prior_tups,
            metrics,
            prior_states: vec![vec![0; 2 * gamma]; s],
            t_idx: 1,
            valid: false,
            valid_really: false,
        }
    }

    pub fn push(&mut self, r: [f32; N]) -> Option<u8> {
        let mut new_metrics = vec![0.0; self.s];
        for ss in 0..self.s {
            let pts: &Vec<_> = &self.prior_tups[ss];
            let mut costs = vec![];
            for (p, _, _, a) in pts {
                let last = self.metrics[*p][self.gamma + self.t_idx - 1];
                let aa: Vec<f32> = a.iter().map(|&aaa| 2.0 * aaa as f32 - 1.0).collect();
                costs.push(llr(&r, &aa) + last);
            }
            let idx = if costs[0] < costs[1] { 0 } else { 1 };
            new_metrics[ss] = costs[idx];
            self.prior_states[ss][self.t_idx] = pts[idx].0;
            self.prior_states[ss][self.t_idx + self.gamma] = pts[idx].0;
        }

        let minval = new_metrics.iter().fold(f32::INFINITY, |a, &b| a.min(b));
        for i in 0..self.s {
            self.metrics[i][self.t_idx] = new_metrics[i] - minval;
            self.metrics[i][self.gamma + self.t_idx] = new_metrics[i] - minval;
        }

        // Lowest cost path at time t
        let mut minidx = 0;
        let mut minval = self.metrics[0][self.t_idx];
        for i in 1..self.s {
            let tmp = self.metrics[i][self.t_idx];
            if tmp < minval {
                minval = tmp;
                minidx = i;
            }
        }

        // Reverse lookup the state transitions for the best path at time t
        let mut states = vec![0; self.gamma + 1];
        states[self.gamma] = minidx;
        for i in 0..self.gamma {
            states[self.gamma - i - 1] =
                self.prior_states[states[self.gamma - i]][self.gamma + self.t_idx - i];
        }

        // Advance our window pointer
        self.t_idx += 1;
        if self.t_idx >= self.gamma {
            self.t_idx = 0;
            self.valid = true;
        }

        // TODO: This is so stupid there has to be a better way.
        if self.valid && self.t_idx > 0 {
            self.valid_really = true;
        }

        // Look up the input bit that causes the best path state transition from
        // gamma time steps ago
        if self.valid_really {
            let sp = states[0];
            if states[1] == self.trellis[2 * sp].1 {
                Some(0)
            } else {
                Some(1)
            }
        } else {
            None
        }
    }

    pub fn push_block(&mut self, r: &[[f32; N]]) -> Vec<u8> {
        let mut output = vec![];
        for rr in r {
            if let Some(y) = self.push(*rr) {
                output.push(y);
            }
        }
        output
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_int_to_bitvec() {
        assert_eq!(int_to_bitvec::<5>(12), [0, 1, 1, 0, 0]);
    }

    #[test]
    fn test_bitvec_to_int() {
        assert_eq!(bitvec_to_int(&[0, 1, 1, 0, 0]), 12);
    }

    #[test]
    fn test_build_trellis() {
        let trellis = build_trellis([[1, 1, 1, 1, 0, 0, 1], [1, 0, 1, 1, 0, 1, 1]]);
        assert_eq!(trellis.len(), 128);
        for i in 0..128 {
            assert_eq!(trellis[i].0, i / 2);
            assert_eq!(trellis[i].2 as usize, i % 2);
        }
    }

    #[test]
    fn convolutional_encoder_constructor() {
        let mut cc =
            ConvolutionalEncoder::<2, 7>::new([[1, 1, 1, 1, 0, 0, 1], [1, 0, 1, 1, 0, 1, 1]]);

        let x = vec![0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1];
        let xlen = x.len();
        let y = cc.push_block(&x);
        assert_eq!(y.len(), 2 * xlen);
    }

    #[test]
    fn short_convolutional_encoder() {
        let mut cc = ConvolutionalEncoder::<2, 3>::new([[1, 0, 1], [1, 1, 1]]);

        let x = vec![1, 1, 0, 0, 1, 0, 1, 0, 0];
        let y = cc.push_block(&x);
        assert_eq!(
            y,
            vec![1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1]
        );
    }
}
