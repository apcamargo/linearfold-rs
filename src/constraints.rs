use crate::{sequence::Base, FoldError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Constraint {
    Unknown,
    Unpaired,
    Paired(usize),
}

#[derive(Debug)]
pub(crate) struct Constraints {
    positions: Vec<Constraint>,
    next_forced_pair: Vec<usize>,
}

impl Constraints {
    pub(crate) fn parse(input: &[u8], sequence: &[Base]) -> Result<Self, FoldError> {
        if input.len() != sequence.len() {
            return Err(FoldError::ConstraintLengthMismatch {
                sequence_len: sequence.len(),
                constraint_len: input.len(),
            });
        }

        let mut positions = vec![Constraint::Unknown; input.len()];
        let mut openings = Vec::new();
        for (index, byte) in input.iter().copied().enumerate() {
            match byte {
                b'?' => positions[index] = Constraint::Unknown,
                b'.' => positions[index] = Constraint::Unpaired,
                b'(' => openings.push(index),
                b')' => {
                    let Some(left) = openings.pop() else {
                        return Err(FoldError::UnmatchedClosingParenthesis { index });
                    };
                    if !sequence[left].can_pair(sequence[index]) {
                        return Err(FoldError::NoncanonicalConstrainedPair { left, right: index });
                    }
                    positions[left] = Constraint::Paired(index);
                    positions[index] = Constraint::Paired(left);
                }
                _ => return Err(FoldError::InvalidConstraintByte { index, byte }),
            }
        }
        if let Some(index) = openings.pop() {
            return Err(FoldError::UnmatchedOpeningParenthesis { index });
        }

        let mut next_forced_pair = vec![input.len(); input.len()];
        let mut next = input.len();
        for index in (0..input.len()).rev() {
            next_forced_pair[index] = next;
            if matches!(positions[index], Constraint::Paired(_)) {
                next = index;
            }
        }

        Ok(Self {
            positions,
            next_forced_pair,
        })
    }

    pub(crate) fn can_unpair(&self, index: usize) -> bool {
        !matches!(self.positions[index], Constraint::Paired(_))
    }

    pub(crate) fn forced_partner(&self, index: usize) -> Option<usize> {
        match self.positions[index] {
            Constraint::Paired(partner) => Some(partner),
            Constraint::Unknown | Constraint::Unpaired => None,
        }
    }

    pub(crate) fn can_pair(&self, left: usize, right: usize) -> bool {
        let left_allowed = match self.positions[left] {
            Constraint::Unknown => true,
            Constraint::Paired(partner) => partner == right,
            Constraint::Unpaired => false,
        };
        let right_allowed = match self.positions[right] {
            Constraint::Unknown => true,
            Constraint::Paired(partner) => partner == left,
            Constraint::Unpaired => false,
        };
        left_allowed && right_allowed
    }

    pub(crate) fn pair_boundary_after(&self, index: usize) -> usize {
        self.next_forced_pair[index]
    }

    pub(crate) fn is_forced_unpaired(&self, index: usize) -> bool {
        matches!(self.positions[index], Constraint::Unpaired)
    }

    pub(crate) fn first_forced_position_at_or_after(&self, index: usize) -> usize {
        if matches!(self.positions[index], Constraint::Paired(_)) {
            index
        } else {
            self.next_forced_pair[index]
        }
    }

    pub(crate) fn range_can_unpair(&self, start: usize, end: usize) -> bool {
        if start >= end {
            return true;
        }
        self.first_forced_position_at_or_after(start) >= end
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::normalize;

    #[test]
    fn parses_pair_and_unpaired_constraints() {
        let bases = normalize(b"AUGU").unwrap();
        let constraints = Constraints::parse(b"(..)", &bases).unwrap();

        assert_eq!(constraints.forced_partner(0), Some(3));
        assert_eq!(constraints.forced_partner(3), Some(0));
        assert_eq!(constraints.forced_partner(1), None);
        assert!(!constraints.can_unpair(0));
        assert!(constraints.can_unpair(1));
        assert!(constraints.can_pair(0, 3));
        assert!(!constraints.can_pair(0, 1));
        assert!(!constraints.range_can_unpair(0, 4));
        assert!(constraints.range_can_unpair(1, 3));
    }

    #[test]
    fn rejects_invalid_constraint_forms() {
        let bases = normalize(b"ACGU").unwrap();
        assert!(matches!(
            Constraints::parse(b"???", &bases),
            Err(FoldError::ConstraintLengthMismatch {
                sequence_len: 4,
                constraint_len: 3
            })
        ));
        assert!(matches!(
            Constraints::parse(b"??x?", &bases),
            Err(FoldError::InvalidConstraintByte {
                index: 2,
                byte: b'x'
            })
        ));
        assert!(matches!(
            Constraints::parse(b"?)??", &bases),
            Err(FoldError::UnmatchedClosingParenthesis { index: 1 })
        ));
        assert!(matches!(
            Constraints::parse(b"?(??", &bases),
            Err(FoldError::UnmatchedOpeningParenthesis { index: 1 })
        ));
        assert!(matches!(
            Constraints::parse(b"(??)", &normalize(b"ACGC").unwrap()),
            Err(FoldError::NoncanonicalConstrainedPair { left: 0, right: 3 })
        ));
    }
}
