mod contrafold;
mod params {
    pub(in crate::model) mod contrafold;
    pub(in crate::model) mod vienna;
}
mod vienna;

use std::{cmp::Ordering, fmt::Debug, ops::Add};

use crate::{DangleModel, sequence::Base};

pub(crate) use contrafold::ContraFoldModel;
pub(crate) use vienna::ViennaModel;

pub(crate) trait Score:
    Copy + Debug + PartialEq + PartialOrd + Add<Output = Self> + Send + Sync + 'static
{
    const ZERO: Self;
    const NEG_INFINITY: Self;

    fn total_cmp(self, other: Self) -> Ordering;
}

impl Score for i32 {
    const ZERO: Self = 0;
    const NEG_INFINITY: Self = i32::MIN;

    fn total_cmp(self, other: Self) -> Ordering {
        self.cmp(&other)
    }
}

impl Score for f64 {
    const ZERO: Self = 0.0;
    const NEG_INFINITY: Self = f64::NEG_INFINITY;

    fn total_cmp(self, other: Self) -> Ordering {
        f64::total_cmp(&self, &other)
    }
}

pub(crate) trait ScoringModel: Sized {
    type Score: Score;

    fn new(sequence: &[Base], dangles: DangleModel) -> Self;
    fn encode(base: Base) -> u8;
    fn hairpin(&self, i: usize, j: usize, nucs: &[u8]) -> Self::Score;
    fn helix(&self, i: usize, j: usize, p: usize, q: usize, nucs: &[u8]) -> Self::Score;
    fn single(&self, i: usize, j: usize, p: usize, q: usize, nucs: &[u8]) -> Self::Score;
    fn multi(&self, i: usize, j: usize, nucs: &[u8]) -> Self::Score;
    fn multi_unpaired(&self, count: usize) -> Self::Score;
    fn m1(&self, i: usize, j: usize, nucs: &[u8]) -> Self::Score;
    fn external_paired(&self, i: usize, j: usize, nucs: &[u8]) -> Self::Score;
    fn external_unpaired(&self, count: usize) -> Self::Score;
}

#[cfg(test)]
mod parameter_tests {
    use super::params::{contrafold as c, vienna as v};

    const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

    fn hash_bytes(mut hash: u64, bytes: &[u8]) -> u64 {
        for byte in bytes {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(FNV_PRIME);
        }
        hash
    }

    #[test]
    fn contrafold_parameter_checksum_matches_reference() {
        let tables: &[&[f64]] = &[
            &c::BASE_PAIR,
            &c::INTERNAL_1X1_NUCLEOTIDES,
            &c::HELIX_STACKING,
            &c::TERMINAL_MISMATCH,
            &c::BULGE_0X1_NUCLEOTIDES,
            &c::HELIX_CLOSING,
            &c::DANGLE_LEFT,
            &c::DANGLE_RIGHT,
            &c::INTERNAL_EXPLICIT,
            &c::HAIRPIN_LENGTH,
            &c::BULGE_LENGTH,
            &c::INTERNAL_LENGTH,
            &c::INTERNAL_SYMMETRIC_LENGTH,
            &c::INTERNAL_ASYMMETRY,
        ];
        let mut hash = FNV_OFFSET;
        for table in tables {
            for value in *table {
                hash = hash_bytes(hash, &value.to_le_bytes());
            }
        }
        assert_eq!(hash, 0xc4f6_86bd_8764_aaf5);
    }

    #[test]
    fn vienna_parameter_checksum_matches_reference() {
        let tables: &[&[i32]] = &[
            &v::TRILOOP37,
            &v::TETRALOOP37,
            &v::HEXALOOP37,
            &v::STACK37,
            &v::HAIRPIN37,
            &v::BULGE37,
            &v::INTERNAL_LOOP37,
            &v::MISMATCH_I37,
            &v::MISMATCH_H37,
            &v::MISMATCH_M37,
            &v::MISMATCH1N_I37,
            &v::MISMATCH23_I37,
            &v::MISMATCH_EXT37,
            &v::DANGLE5_37,
            &v::DANGLE3_37,
            &v::INT11_37,
            &v::INT21_37,
            &v::INT22_37,
        ];
        let mut hash = FNV_OFFSET;
        for table in tables {
            for value in *table {
                hash = hash_bytes(hash, &value.to_le_bytes());
            }
        }
        assert_eq!(hash, 0x2d15_0dd4_97ad_583a);
    }
}
