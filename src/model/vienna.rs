use crate::{DangleModel, sequence::Base};

use super::{
    ScoringModel,
    params::vienna::{
        BULGE37, DANGLE3_37, DANGLE5_37, HAIRPIN37, HEXALOOP37, INT11_37, INT21_37, INT22_37,
        INTERNAL_LOOP37, LXC37, MAX_NINIO, MISMATCH_EXT37, MISMATCH_H37, MISMATCH_I37,
        MISMATCH_M37, MISMATCH1N_I37, MISMATCH23_I37, ML_CLOSING37, ML_INTERN37, NINIO37, STACK37,
        TERMINAL_AU37, TETRALOOP37, TRILOOP37,
    },
};

const MAX_LOOP: usize = 30;
const TRI_LOOPS: [&[u8]; 2] = [b"CAACG", b"GUUAC"];
const TETRA_LOOPS: [&[u8]; 16] = [
    b"CAACGG", b"CCAAGG", b"CCACGG", b"CCCAGG", b"CCGAGG", b"CCGCGG", b"CCUAGG", b"CCUCGG",
    b"CUAAGG", b"CUACGG", b"CUCAGG", b"CUCCGG", b"CUGCGG", b"CUUAGG", b"CUUCGG", b"CUUUGG",
];
const HEXA_LOOPS: [&[u8]; 4] = [b"ACAGUACU", b"ACAGUGAU", b"ACAGUGCU", b"ACAGUGUU"];

pub(crate) struct ViennaModel {
    dangles: DangleModel,
    sequence: Vec<u8>,
}

impl ViennaModel {
    fn pair_type(left: u8, right: u8) -> usize {
        match (left, right) {
            (1, 4) => 5,
            (2, 3) => 1,
            (3, 2) => 2,
            (3, 4) => 3,
            (4, 3) => 4,
            (4, 1) => 6,
            _ => 0,
        }
    }

    fn table3(table: &[i32], a: usize, b: u8, c: u8) -> i32 {
        table[a * 25 + usize::from(b) * 5 + usize::from(c)]
    }

    fn stack(first: usize, second: usize) -> i32 {
        STACK37[first * 8 + second]
    }

    fn special_loop_from_motif(motif: &[u8], size: usize) -> Option<i32> {
        match size {
            3 => TRI_LOOPS
                .iter()
                .position(|candidate| motif == *candidate)
                .map(|index| TRILOOP37[index]),
            4 => TETRA_LOOPS
                .iter()
                .position(|candidate| motif == *candidate)
                .map(|index| TETRALOOP37[index]),
            6 => HEXA_LOOPS
                .iter()
                .position(|candidate| motif == *candidate)
                .map(|index| HEXALOOP37[index]),
            _ => None,
        }
    }

    #[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
    fn loop_energy(table: &[i32; 31], size: usize) -> i32 {
        if size <= MAX_LOOP {
            table[size]
        } else {
            table[MAX_LOOP] + (LXC37 * (size as f64 / MAX_LOOP as f64).ln()) as i32
        }
    }

    fn hairpin_core(size: usize, pair_type: usize, si1: u8, sj1: u8, motif: Option<&[u8]>) -> i32 {
        let energy = Self::loop_energy(&HAIRPIN37, size);
        if size < 3 {
            return energy;
        }
        if let Some(motif) = motif
            && let Some(special) = Self::special_loop_from_motif(motif, size)
        {
            return special;
        }
        if size == 3 {
            return energy + if pair_type > 2 { TERMINAL_AU37 } else { 0 };
        }
        energy + Self::table3(&MISMATCH_H37, pair_type, si1, sj1)
    }

    fn hairpin_energy(&self, i: usize, j: usize, nucs: &[u8]) -> i32 {
        let size = j - i - 1;
        let pair_type = Self::pair_type(nucs[i], nucs[j]);
        let motif = self.sequence.get(i..j + 1);
        Self::hairpin_core(size, pair_type, nucs[i + 1], nucs[j - 1], motif)
    }

    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn degree_two_energy(
        n1: usize,
        n2: usize,
        first_type: usize,
        second_type: usize,
        si1: u8,
        sj1: u8,
        sp1: u8,
        sq1: u8,
    ) -> i32 {
        let (long, short) = if n1 > n2 { (n1, n2) } else { (n2, n1) };

        if long == 0 {
            return Self::stack(first_type, second_type);
        }
        if short == 0 {
            let mut energy = Self::loop_energy(&BULGE37, long);
            if long == 1 {
                energy += Self::stack(first_type, second_type);
            } else {
                if first_type > 2 {
                    energy += TERMINAL_AU37;
                }
                if second_type > 2 {
                    energy += TERMINAL_AU37;
                }
            }
            return energy;
        }

        if short == 1 {
            if long == 1 {
                let index = (((first_type * 8 + second_type) * 5 + usize::from(si1)) * 5)
                    + usize::from(sj1);
                return INT11_37[index];
            }
            if long == 2 {
                let index = if n1 == 1 {
                    ((((first_type * 8 + second_type) * 5 + usize::from(si1)) * 5
                        + usize::from(sq1))
                        * 5)
                        + usize::from(sj1)
                } else {
                    ((((second_type * 8 + first_type) * 5 + usize::from(sq1)) * 5
                        + usize::from(si1))
                        * 5)
                        + usize::from(sp1)
                };
                return INT21_37[index];
            }
            let mut energy = Self::loop_energy(&INTERNAL_LOOP37, long + 1);
            energy += MAX_NINIO.min((long - short) as i32 * NINIO37);
            energy += Self::table3(&MISMATCH1N_I37, first_type, si1, sj1);
            energy += Self::table3(&MISMATCH1N_I37, second_type, sq1, sp1);
            return energy;
        }

        if short == 2 {
            if long == 2 {
                let index = (((((first_type * 8 + second_type) * 5 + usize::from(si1)) * 5
                    + usize::from(sp1))
                    * 5
                    + usize::from(sq1))
                    * 5)
                    + usize::from(sj1);
                return INT22_37[index];
            }
            if long == 3 {
                return INTERNAL_LOOP37[5]
                    + NINIO37
                    + Self::table3(&MISMATCH23_I37, first_type, si1, sj1)
                    + Self::table3(&MISMATCH23_I37, second_type, sq1, sp1);
            }
        }

        Self::loop_energy(&INTERNAL_LOOP37, long + short)
            + MAX_NINIO.min((long - short) as i32 * NINIO37)
            + Self::table3(&MISMATCH_I37, first_type, si1, sj1)
            + Self::table3(&MISMATCH_I37, second_type, sq1, sp1)
    }

    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn single_energy(i: usize, j: usize, p: usize, q: usize, nucs: &[u8]) -> i32 {
        let first_type = Self::pair_type(nucs[i], nucs[j]);
        let second_type = Self::pair_type(nucs[q], nucs[p]);
        let n1 = p - i - 1;
        let n2 = j - q - 1;
        Self::degree_two_energy(
            n1,
            n2,
            first_type,
            second_type,
            nucs[i + 1],
            nucs[j - 1],
            nucs[p - 1],
            nucs[q + 1],
        )
    }

    pub(crate) fn circular_unpaired_score(length: usize) -> i32 {
        let kt = (37.0 + 273.15) * 1.98717 / 1000.0;
        let penalty = (100.0 * kt * (4.385 + 1.5 * (length as f64).ln()) + 0.5).floor() as i32;
        -penalty
    }

    pub(crate) fn circular_hairpin_score(&self, i: usize, j: usize, nucs: &[u8]) -> i32 {
        let n = self.sequence.len();
        let u = n - j + i - 1;
        if u < 3 {
            return i32::MIN;
        }
        let pair_type = Self::pair_type(nucs[j], nucs[i]);
        let si1 = nucs[(j + 1) % n];
        let sj1 = nucs[(i + n - 1) % n];

        let mut motif_buf = [0u8; 8];
        let motif = if u == 3 || u == 4 || u == 6 {
            let total_len = u + 2;
            for (idx, byte) in motif_buf.iter_mut().enumerate().take(total_len) {
                let seq_idx = (j + idx) % n;
                *byte = self.sequence[seq_idx];
            }
            Some(&motif_buf[..total_len])
        } else {
            None
        };

        -Self::hairpin_core(u, pair_type, si1, sj1, motif)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn circular_internal_score(
        i: usize,
        j: usize,
        p: usize,
        q: usize,
        nucs: &[u8],
    ) -> i32 {
        let n = nucs.len();
        let n1 = p - j - 1;
        let n2 = i + n - q - 1;
        let first_type = Self::pair_type(nucs[j], nucs[i]);
        let second_type = Self::pair_type(nucs[q], nucs[p]);
        let si1 = nucs[j + 1];
        let sj1 = nucs[(i + n - 1) % n];
        let sp1 = nucs[p - 1];
        let sq1 = nucs[(q + 1) % n];

        -Self::degree_two_energy(n1, n2, first_type, second_type, si1, sj1, sp1, sq1)
    }

    pub(crate) fn circular_multiloop_stem_score(&self, i: usize, j: usize, nucs: &[u8]) -> i32 {
        let n = nucs.len();
        let pair_type = Self::pair_type(nucs[i], nucs[j]);
        let left = nucs[(i + n - 1) % n];
        let right = nucs[(j + 1) % n];
        -self.ml_stem(pair_type, Some(left), Some(right))
    }

    pub(crate) const fn circular_multiloop_closing_score() -> i32 {
        -ML_CLOSING37
    }

    fn ml_stem(&self, pair_type: usize, left: Option<u8>, right: Option<u8>) -> i32 {
        let mut energy = 0;
        if let (DangleModel::Both, Some(left), Some(right)) = (self.dangles, left, right) {
            energy += Self::table3(&MISMATCH_M37, pair_type, left, right);
        }
        if pair_type > 2 {
            energy += TERMINAL_AU37;
        }
        energy + ML_INTERN37
    }

    fn external_stem(&self, pair_type: usize, left: Option<u8>, right: Option<u8>) -> i32 {
        let mut energy = 0;
        if self.dangles == DangleModel::Both {
            energy += match (left, right) {
                (Some(left), Some(right)) => Self::table3(&MISMATCH_EXT37, pair_type, left, right),
                (Some(left), None) => DANGLE5_37[pair_type * 5 + usize::from(left)],
                (None, Some(right)) => DANGLE3_37[pair_type * 5 + usize::from(right)],
                (None, None) => 0,
            };
        }
        if pair_type > 2 {
            energy += TERMINAL_AU37;
        }
        energy
    }
}

impl ScoringModel for ViennaModel {
    type Score = i32;

    fn new(sequence: &[Base], dangles: DangleModel) -> Self {
        let sequence = sequence
            .iter()
            .map(|base| match base {
                Base::A => b'A',
                Base::C => b'C',
                Base::G => b'G',
                Base::U => b'U',
                Base::Unknown => b'N',
            })
            .collect();
        Self { dangles, sequence }
    }

    fn encode(base: Base) -> u8 {
        match base {
            Base::Unknown => 0,
            Base::A => 1,
            Base::C => 2,
            Base::G => 3,
            Base::U => 4,
        }
    }

    fn hairpin(&self, i: usize, j: usize, nucs: &[u8]) -> i32 {
        -self.hairpin_energy(i, j, nucs)
    }

    fn helix(&self, i: usize, j: usize, p: usize, q: usize, nucs: &[u8]) -> i32 {
        -Self::single_energy(i, j, p, q, nucs)
    }

    fn single(&self, i: usize, j: usize, p: usize, q: usize, nucs: &[u8]) -> i32 {
        -Self::single_energy(i, j, p, q, nucs)
    }

    fn multi(&self, i: usize, j: usize, nucs: &[u8]) -> i32 {
        let pair_type = Self::pair_type(nucs[j], nucs[i]);
        -(self.ml_stem(pair_type, Some(nucs[j - 1]), Some(nucs[i + 1])) + ML_CLOSING37)
    }

    fn multi_unpaired(&self, _count: usize) -> i32 {
        0
    }

    fn m1(&self, i: usize, j: usize, nucs: &[u8]) -> i32 {
        let pair_type = Self::pair_type(nucs[i], nucs[j]);
        -self.ml_stem(
            pair_type,
            i.checked_sub(1).map(|index| nucs[index]),
            nucs.get(j + 1).copied(),
        )
    }

    fn external_paired(&self, i: usize, j: usize, nucs: &[u8]) -> i32 {
        let pair_type = Self::pair_type(nucs[i], nucs[j]);
        -self.external_stem(
            pair_type,
            i.checked_sub(1).map(|index| nucs[index]),
            nucs.get(j + 1).copied(),
        )
    }

    fn external_unpaired(&self, _count: usize) -> i32 {
        0
    }
}
