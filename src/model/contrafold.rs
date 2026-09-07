use crate::{DangleModel, sequence::Base};

use super::{
    ScoringModel,
    params::contrafold::{
        BASE_PAIR, BULGE_0X1_NUCLEOTIDES, BULGE_LENGTH, DANGLE_LEFT, DANGLE_RIGHT, EXTERNAL_PAIRED,
        EXTERNAL_UNPAIRED, HAIRPIN_LENGTH, HELIX_CLOSING, HELIX_STACKING, INTERNAL_1X1_NUCLEOTIDES,
        INTERNAL_ASYMMETRY, INTERNAL_EXPLICIT, INTERNAL_LENGTH, INTERNAL_SYMMETRIC_LENGTH,
        MULTI_BASE, MULTI_PAIRED, MULTI_UNPAIRED, TERMINAL_MISMATCH,
    },
};

const MAX_LOOP: usize = 30;
const EXPLICIT_MAX: usize = 4;
const SYMMETRIC_MAX: usize = 15;
const ASYMMETRY_MAX: usize = 28;

const fn build_single_cache() -> [[f64; MAX_LOOP + 1]; MAX_LOOP + 1] {
    let mut cache = [[0.0; MAX_LOOP + 1]; MAX_LOOP + 1];
    let mut left = 0;
    while left <= MAX_LOOP {
        let mut right = 0;
        while right <= MAX_LOOP {
            if left == 0 && right != 0 {
                cache[left][right] += BULGE_LENGTH[right];
            } else if right == 0 && left != 0 {
                cache[left][right] += BULGE_LENGTH[left];
            } else if left != 0 {
                let total = if left + right < MAX_LOOP {
                    left + right
                } else {
                    MAX_LOOP
                };
                cache[left][right] += INTERNAL_LENGTH[total];
                if left <= EXPLICIT_MAX && right <= EXPLICIT_MAX {
                    let index = if left <= right {
                        left * EXPLICIT_MAX + right
                    } else {
                        right * EXPLICIT_MAX + left
                    };
                    cache[left][right] += INTERNAL_EXPLICIT[index];
                }
                if left == right {
                    let index = if left < SYMMETRIC_MAX {
                        left
                    } else {
                        SYMMETRIC_MAX
                    };
                    cache[left][right] += INTERNAL_SYMMETRIC_LENGTH[index];
                } else {
                    let difference = left.abs_diff(right);
                    let index = if difference < ASYMMETRY_MAX {
                        difference
                    } else {
                        ASYMMETRY_MAX
                    };
                    cache[left][right] += INTERNAL_ASYMMETRY[index];
                }
            }
            right += 1;
        }
        left += 1;
    }
    cache
}

const SINGLE_CACHE: [[f64; MAX_LOOP + 1]; MAX_LOOP + 1] = build_single_cache();

pub(crate) struct ContraFoldModel;

impl ContraFoldModel {
    fn pair(left: u8, right: u8) -> f64 {
        BASE_PAIR[usize::from(right) * 5 + usize::from(left)]
    }

    fn stacking(i: u8, i1: u8, j_1: u8, j: u8) -> f64 {
        HELIX_STACKING
            [usize::from(i) * 125 + usize::from(j) * 25 + usize::from(i1) * 5 + usize::from(j_1)]
    }

    fn closing(left: u8, right: u8) -> f64 {
        HELIX_CLOSING[usize::from(left) * 5 + usize::from(right)]
    }

    fn mismatch(i: u8, i1: u8, j_1: u8, j: u8) -> f64 {
        TERMINAL_MISMATCH
            [usize::from(i) * 125 + usize::from(j) * 25 + usize::from(i1) * 5 + usize::from(j_1)]
    }

    fn dangle_left(i: u8, i1: u8, j: u8) -> f64 {
        DANGLE_LEFT[usize::from(i) * 25 + usize::from(j) * 5 + usize::from(i1)]
    }

    fn dangle_right(i: u8, j_1: u8, j: u8) -> f64 {
        DANGLE_RIGHT[usize::from(i) * 25 + usize::from(j) * 5 + usize::from(j_1)]
    }

    fn junction_a(i: usize, j: usize, nucs: &[u8]) -> f64 {
        let left = nucs[i];
        let right = nucs[j];
        Self::closing(left, right)
            + nucs
                .get(i + 1)
                .map_or(0.0, |&next| Self::dangle_left(left, next, right))
            + j.checked_sub(1)
                .map_or(0.0, |prev| Self::dangle_right(left, nucs[prev], right))
    }

    fn junction_b(i: usize, j: usize, nucs: &[u8]) -> f64 {
        Self::closing(nucs[i], nucs[j]) + Self::mismatch(nucs[i], nucs[i + 1], nucs[j - 1], nucs[j])
    }
}

impl ScoringModel for ContraFoldModel {
    type Score = f64;

    fn new(_sequence: &[Base], _dangles: DangleModel) -> Self {
        Self
    }

    fn encode(base: Base) -> u8 {
        base as u8
    }

    fn hairpin(&self, i: usize, j: usize, nucs: &[u8]) -> f64 {
        HAIRPIN_LENGTH[(j - i - 1).min(MAX_LOOP)] + Self::junction_b(i, j, nucs)
    }

    fn helix(&self, i: usize, j: usize, p: usize, q: usize, nucs: &[u8]) -> f64 {
        Self::stacking(nucs[i], nucs[i + 1], nucs[j - 1], nucs[j]) + Self::pair(nucs[p], nucs[q])
    }

    fn single(&self, i: usize, j: usize, p: usize, q: usize, nucs: &[u8]) -> f64 {
        let left = p - i - 1;
        let right = j - q - 1;
        let nucleotide = match (left, right) {
            (0, 1) => BULGE_0X1_NUCLEOTIDES[usize::from(nucs[q + 1])],
            (1, 0) => BULGE_0X1_NUCLEOTIDES[usize::from(nucs[p - 1])],
            (1, 1) => {
                INTERNAL_1X1_NUCLEOTIDES[usize::from(nucs[p - 1]) * 5 + usize::from(nucs[q + 1])]
            }
            _ => 0.0,
        };
        (Self::junction_b(i, j, nucs) + Self::junction_b(q, p, nucs))
            + ((SINGLE_CACHE[left][right] + Self::pair(nucs[p], nucs[q])) + nucleotide)
    }

    fn multi(&self, i: usize, j: usize, nucs: &[u8]) -> f64 {
        Self::junction_a(i, j, nucs) + MULTI_PAIRED + MULTI_BASE
    }

    #[allow(clippy::cast_precision_loss)]
    fn multi_unpaired(&self, count: usize) -> f64 {
        count as f64 * MULTI_UNPAIRED
    }

    fn m1(&self, i: usize, j: usize, nucs: &[u8]) -> f64 {
        let junction = Self::closing(nucs[j], nucs[i])
            + nucs
                .get(j + 1)
                .map_or(0.0, |&next| Self::dangle_left(nucs[j], next, nucs[i]))
            + i.checked_sub(1).map_or(0.0, |previous| {
                Self::dangle_right(nucs[j], nucs[previous], nucs[i])
            });
        junction + Self::pair(nucs[i], nucs[j]) + MULTI_PAIRED
    }

    fn external_paired(&self, i: usize, j: usize, nucs: &[u8]) -> f64 {
        let junction = Self::closing(nucs[j], nucs[i])
            + nucs
                .get(j + 1)
                .map_or(0.0, |&next| Self::dangle_left(nucs[j], next, nucs[i]))
            + i.checked_sub(1).map_or(0.0, |previous| {
                Self::dangle_right(nucs[j], nucs[previous], nucs[i])
            });
        junction + EXTERNAL_PAIRED + Self::pair(nucs[i], nucs[j])
    }

    #[allow(clippy::cast_precision_loss)]
    fn external_unpaired(&self, count: usize) -> f64 {
        count as f64 * EXTERNAL_UNPAIRED
    }
}
