#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]

mod constraints;
mod engine;
mod model;
mod sequence;

use constraints::Constraints;
use model::{ContraFoldModel, ViennaModel};

/// The scoring model used by the folding engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Model {
    /// Vienna `RNAfold` thermodynamic parameters.
    ViennaRnafold,
    /// `CONTRAfold` log-linear parameters.
    ContraFold,
}

/// Treatment of dangling bases adjacent to helices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DangleModel {
    /// Ignore dangling-end energies.
    None,
    /// Include both adjacent dangling bases when available.
    Both,
}

/// Sequence topology for secondary structure prediction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceTopology {
    /// Standard linear sequence with open free ends.
    Linear,
    /// Covalently closed circular sequence.
    Circular,
}

/// Configuration for one optimal fold.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FoldOptions {
    /// Scoring model. The default is [`Model::ContraFold`].
    pub model: Model,
    /// Maximum states retained per beam. Zero disables pruning.
    pub beam_size: usize,
    /// Allow pairs separated by fewer than three unpaired bases.
    pub allow_sharp_turns: bool,
    /// Vienna dangling-end treatment.
    pub dangles: DangleModel,
    /// Topology of the RNA sequence.
    pub topology: SequenceTopology,
}

impl Default for FoldOptions {
    fn default() -> Self {
        Self {
            model: Model::ContraFold,
            beam_size: 100,
            allow_sharp_turns: false,
            dangles: DangleModel::Both,
            topology: SequenceTopology::Linear,
        }
    }
}

/// Model-specific score for a predicted structure.
#[derive(Debug, Clone, PartialEq)]
pub enum PredictionScore {
    /// Vienna free energy in kcal/mol.
    ViennaRnafold(f64),
    /// `CONTRAfold` log-linear score.
    ContraFold(f64),
}

/// Optimal secondary-structure prediction.
#[derive(Debug, Clone, PartialEq)]
pub struct Prediction {
    /// Dot-bracket structure with one byte per input base.
    pub structure: String,
    /// Score whose variant identifies its units.
    pub score: PredictionScore,
}

/// Failure returned by [`fold`].
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum FoldError {
    #[error("the sequence is empty")]
    EmptySequence,
    #[error("invalid sequence byte {byte:#04x} at index {index}")]
    InvalidSequenceByte { index: usize, byte: u8 },
    #[error("constraint length {constraint_len} does not match sequence length {sequence_len}")]
    ConstraintLengthMismatch {
        sequence_len: usize,
        constraint_len: usize,
    },
    #[error("invalid constraint byte {byte:#04x} at index {index}")]
    InvalidConstraintByte { index: usize, byte: u8 },
    #[error("closing parenthesis at index {index} has no opener")]
    UnmatchedClosingParenthesis { index: usize },
    #[error("opening parenthesis at index {index} has no closer")]
    UnmatchedOpeningParenthesis { index: usize },
    #[error("forced pair {left}..{right} is not AU, CG, or GU")]
    NoncanonicalConstrainedPair { left: usize, right: usize },
    #[error("no structure satisfies the supplied constraints")]
    NoValidStructure,
    #[error("circular folding is supported only by the Vienna RNAfold model")]
    CircularFoldingRequiresViennaRnafold,
    #[error("internal folding invariant failed: {0}")]
    InternalInvariant(&'static str),
}

/// Predicts the optimal secondary structure for one RNA or DNA sequence.
///
/// A beam size of zero disables pruning. Constraints use `?`, `.`, `(`, and `)`.
///
/// # Errors
///
/// Returns [`FoldError`] for invalid input, malformed constraints, or an
/// unsatisfiable constraint set.
pub fn fold(
    sequence: &[u8],
    constraints: Option<&[u8]>,
    options: &FoldOptions,
) -> Result<Prediction, FoldError> {
    let bases = sequence::normalize(sequence)?;
    let constraints = constraints
        .map(|input| Constraints::parse(input, &bases))
        .transpose()?;

    if options.topology == SequenceTopology::Circular && options.model != Model::ViennaRnafold {
        return Err(FoldError::CircularFoldingRequiresViennaRnafold);
    }

    match options.model {
        Model::ViennaRnafold => {
            let filled = engine::fill_charts::<ViennaModel>(&bases, constraints.as_ref(), options)?;
            let (structure, score) = match options.topology {
                SequenceTopology::Linear => engine::run_linear_from_filled(&filled, bases.len())?,
                SequenceTopology::Circular => {
                    engine::run_circular_vienna(filled, constraints.as_ref(), bases.len())?
                }
            };
            Ok(Prediction {
                structure,
                score: PredictionScore::ViennaRnafold(-(f64::from(score)) / 100.0),
            })
        }
        Model::ContraFold => {
            let (structure, score) =
                engine::run_linear::<ContraFoldModel>(&bases, constraints.as_ref(), options)?;
            Ok(Prediction {
                structure,
                score: PredictionScore::ContraFold(score),
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_the_public_contract() {
        let options = FoldOptions::default();
        assert_eq!(options.model, Model::ContraFold);
        assert_eq!(options.beam_size, 100);
        assert!(!options.allow_sharp_turns);
        assert_eq!(options.dangles, DangleModel::Both);
        assert_eq!(options.topology, SequenceTopology::Linear);
    }

    #[test]
    fn validation_precedes_topology_checks() {
        let options = FoldOptions {
            model: Model::ContraFold,
            topology: SequenceTopology::Circular,
            ..FoldOptions::default()
        };
        assert_eq!(fold(b"", None, &options), Err(FoldError::EmptySequence));
        assert_eq!(
            fold(b"AC GU", None, &options),
            Err(FoldError::InvalidSequenceByte {
                index: 2,
                byte: b' '
            })
        );
        assert_eq!(
            fold(b"ACGU", Some(b"???"), &options),
            Err(FoldError::ConstraintLengthMismatch {
                sequence_len: 4,
                constraint_len: 3
            })
        );
        assert_eq!(
            fold(b"ACGU", Some(b"??x?"), &options),
            Err(FoldError::InvalidConstraintByte {
                index: 2,
                byte: b'x'
            })
        );
        assert_eq!(
            fold(b"ACGU", None, &options),
            Err(FoldError::CircularFoldingRequiresViennaRnafold)
        );
    }

    #[test]
    fn rejects_malformed_and_unsatisfiable_constraints() {
        let options = FoldOptions::default();
        assert_eq!(
            fold(b"ACGU", Some(b"?)??"), &options),
            Err(FoldError::UnmatchedClosingParenthesis { index: 1 })
        );
        assert_eq!(
            fold(b"ACGU", Some(b"?(??"), &options),
            Err(FoldError::UnmatchedOpeningParenthesis { index: 1 })
        );
        assert_eq!(
            fold(b"ACGC", Some(b"(??)"), &options),
            Err(FoldError::NoncanonicalConstrainedPair { left: 0, right: 3 })
        );
        assert_eq!(
            fold(b"AU", Some(b"()"), &options)
                .expect("forced pairs allow sharp turns")
                .structure,
            "()"
        );

        let sequence = format!("A{}A{}UU", "C".repeat(31), "C".repeat(6));
        let constraints = format!("({}({}))", ".".repeat(31), ".".repeat(6));
        assert_eq!(
            fold(sequence.as_bytes(), Some(constraints.as_bytes()), &options),
            Err(FoldError::NoValidStructure)
        );
    }

    #[test]
    fn contrafold_known_examples_remain_stable() {
        let examples = [
            (
                b"UGAGUUCUCGAUCUCUAAAAUCG".as_slice(),
                ".......................",
                -0.22,
            ),
            (
                b"AAAACGGUCCUUAUCAGGACCAAACA".as_slice(),
                ".....((((((....)))))).....",
                4.91,
            ),
            (
                b"UCGGCCACAAACACACAAUCUACUGUUGGUCGA".as_slice(),
                "(((((((...................)))))))",
                0.99,
            ),
        ];

        for (sequence, structure, score) in examples {
            let prediction = fold(sequence, None, &FoldOptions::default()).unwrap();
            assert_eq!(prediction.structure, structure);
            let PredictionScore::ContraFold(actual_score) = prediction.score else {
                panic!("wrong score model");
            };
            assert!((actual_score - score).abs() < 0.005);
        }
    }

    #[test]
    fn vienna_known_examples_remain_stable() {
        let options = FoldOptions {
            model: Model::ViennaRnafold,
            ..FoldOptions::default()
        };
        let examples = [
            (
                b"UGAGUUCUCGAUCUCUAAAAUCG".as_slice(),
                ".(((........)))........",
                -1.80,
            ),
            (
                b"AAAACGGUCCUUAUCAGGACCAAACA".as_slice(),
                ".....((((((....)))))).....",
                -9.30,
            ),
            (
                b"AUUCUUGCUUCAACAGUGUUUGAACGGAAU".as_slice(),
                "(((((...(((((......))))).)))))",
                -6.80,
            ),
            (
                b"UCGGCCACAAACACACAAUCUACUGUUGGUCGA".as_slice(),
                "(((((((((..............))).))))))",
                -7.80,
            ),
            (
                b"GUUUUUAUCUUACACACGCUUGUGUAAGAUAGUUA".as_slice(),
                "....((((((((((((....))))))))))))...",
                -13.00,
            ),
        ];

        for (sequence, structure, score) in examples {
            let prediction = fold(sequence, None, &options).unwrap();
            assert_eq!(prediction.structure, structure);
            assert_eq!(prediction.score, PredictionScore::ViennaRnafold(score));
        }
    }

    #[test]
    fn constrained_predictions_remain_stable() {
        let sequence = b"AACUCCGCCAGGCCUGGAAGGGAGCAACGGUAGUGACACUCUCUGUGUGCGUAGGUUGCCUAGCUACCAUUU";
        let constraints =
            b"??(???(??????)?(????????)???(??????(???????)?)???????????)??.???????????";
        let prediction = fold(sequence, Some(constraints), &FoldOptions::default()).unwrap();
        assert_eq!(
            prediction.structure,
            "..(.(((......)((........))(((......(.......).))).....))..).............."
        );
        let PredictionScore::ContraFold(score) = prediction.score else {
            panic!("wrong score model");
        };
        assert!((score + 27.33).abs() < 0.005);

        let vienna = fold(
            sequence,
            Some(constraints),
            &FoldOptions {
                model: Model::ViennaRnafold,
                ..FoldOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            vienna.structure,
            "..(.(((......)((........))(((......(.......).))).....))..).............."
        );
        assert_eq!(vienna.score, PredictionScore::ViennaRnafold(13.40));
    }

    #[test]
    fn options_change_predictions_without_affecting_contrafold_dangles() {
        let sequence = b"GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA";
        let beam_twenty = FoldOptions {
            model: Model::ViennaRnafold,
            beam_size: 20,
            ..FoldOptions::default()
        };
        let prediction = fold(sequence, None, &beam_twenty).unwrap();
        assert_eq!(prediction.score, PredictionScore::ViennaRnafold(-31.50));

        let beam_one = FoldOptions {
            model: Model::ViennaRnafold,
            beam_size: 1,
            ..FoldOptions::default()
        };
        let beam_one_prediction = fold(sequence, None, &beam_one).unwrap();
        assert_eq!(
            beam_one_prediction.structure,
            ".(((((((......))))......)))..............((((((...........))).).))......"
        );
        assert_eq!(
            beam_one_prediction.score,
            PredictionScore::ViennaRnafold(-3.10)
        );

        let no_dangles = FoldOptions {
            model: Model::ViennaRnafold,
            dangles: DangleModel::None,
            ..FoldOptions::default()
        };
        let no_dangles_prediction = fold(sequence, None, &no_dangles).unwrap();
        assert_eq!(
            no_dangles_prediction.score,
            PredictionScore::ViennaRnafold(-25.50)
        );

        let unbounded = FoldOptions {
            model: Model::ViennaRnafold,
            beam_size: 0,
            ..FoldOptions::default()
        };
        assert_eq!(
            fold(sequence, None, &unbounded).unwrap().score,
            PredictionScore::ViennaRnafold(-31.50)
        );

        let contra_default = fold(sequence, None, &FoldOptions::default()).unwrap();
        let contra_no_dangles = fold(
            sequence,
            None,
            &FoldOptions {
                dangles: DangleModel::None,
                ..FoldOptions::default()
            },
        )
        .unwrap();
        assert_eq!(contra_no_dangles, contra_default);

        let sharp_sequence = b"GAUGUCAAACCCCGGGGGGA";
        let without_sharp_turns = fold(sharp_sequence, None, &FoldOptions::default()).unwrap();
        let with_sharp_turns = fold(
            sharp_sequence,
            None,
            &FoldOptions {
                allow_sharp_turns: true,
                ..FoldOptions::default()
            },
        )
        .unwrap();
        assert_eq!(without_sharp_turns.structure, ".........(((....))).");
        assert_eq!(with_sharp_turns.structure, ".........(((())))...");
    }

    #[test]
    fn long_unpairable_gaps_do_not_overflow_trace_padding() {
        let mut sequence =
            b"AACUCCGCCAGGCCUGGAAGGGAGCAACGGUAGUGACACUCUCUGUGUGCGUAGGUUGCCUAGCUACCAUUU".to_vec();
        sequence.extend(std::iter::repeat_n(b'N', 300));
        sequence.push(b'U');

        let prediction = fold(&sequence, None, &FoldOptions::default()).unwrap();
        assert_eq!(prediction.structure.len(), sequence.len());
        let PredictionScore::ContraFold(score) = prediction.score else {
            panic!("wrong score model");
        };
        assert!((score + 2.47).abs() < 0.005);
    }

    #[test]
    fn circular_vienna_cases_remain_stable() {
        let options = FoldOptions {
            model: Model::ViennaRnafold,
            topology: SequenceTopology::Circular,
            beam_size: 0,
            ..FoldOptions::default()
        };
        let cases = [
            (b"A".as_slice(), ".", 2.70),
            (b"GGGAAACCC".as_slice(), ".........", 4.73),
            (
                b"AUCGAUCGAUCGAUCGAUCG".as_slice(),
                ".((((((....))))))...",
                -1.10,
            ),
            (
                b"GGGAAACCCGGGAAACCC".as_slice(),
                "(((...)))(((...)))",
                -4.80,
            ),
            (
                b"AGGGGGAAAAAAAACCCCCA".as_slice(),
                "..((((........))))..",
                -2.20,
            ),
        ];
        for (sequence, structure, score) in cases {
            let prediction = fold(sequence, None, &options).unwrap();
            assert_eq!(prediction.structure, structure);
            assert_eq!(prediction.score, PredictionScore::ViennaRnafold(score));
        }

        let forced_unpaired = fold(b"GGGAAACCC", Some(b"........."), &options).unwrap();
        assert_eq!(forced_unpaired.structure, ".........");
        assert_eq!(forced_unpaired.score, PredictionScore::ViennaRnafold(4.73));

        let short_loop_options = FoldOptions {
            allow_sharp_turns: true,
            ..options
        };
        assert_eq!(
            fold(b"GGGAAACCC", Some(b"(((...)))"), &short_loop_options),
            Err(FoldError::NoValidStructure)
        );

        let degree_two_cases = [
            (
                b"AGGGGGAAAACCCCCAGGGGGAAAACCCCC".as_slice(),
                ".(((((....))))).(((((....)))))",
                -16.50,
            ),
            (
                b"AAGGGGGAAAACCCCCAAGGGGGAAAACCCCC".as_slice(),
                "..(((((....)))))..(((((....)))))",
                -16.20,
            ),
            (
                b"GGGGGAAAACCCCCAAAGGGGGAAAACCCCC".as_slice(),
                "(((((....)))))...(((((....)))))",
                -14.20,
            ),
        ];
        for (sequence, structure, score) in degree_two_cases {
            let prediction = fold(sequence, None, &options).unwrap();
            assert_eq!(prediction.structure, structure);
            assert_eq!(prediction.score, PredictionScore::ViennaRnafold(score));
        }
    }

    #[test]
    fn circular_multiloop_dangles_change_only_the_score() {
        let sequence = b"GCGCAAAAGCGCCGCGAAAACGCGGGCCAAAAGGCC";
        for (dangles, score) in [(DangleModel::Both, -12.30), (DangleModel::None, -8.10)] {
            let options = FoldOptions {
                model: Model::ViennaRnafold,
                topology: SequenceTopology::Circular,
                dangles,
                beam_size: 0,
                ..FoldOptions::default()
            };
            let prediction = fold(sequence, None, &options).unwrap();
            assert_eq!(prediction.structure, "((((....))))((((....))))((((....))))");
            assert_eq!(prediction.score, PredictionScore::ViennaRnafold(score));
        }
    }

    #[test]
    fn predictions_are_balanced_and_canonical_for_representative_sequences() {
        let sequences = [
            b"ACGU".as_slice(),
            b"GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCC".as_slice(),
            b"ACGRUN".as_slice(),
            b"A".as_slice(),
        ];
        for sequence in sequences {
            let prediction = fold(sequence, None, &FoldOptions::default()).unwrap();
            let mut openings = Vec::new();
            for (index, symbol) in prediction.structure.bytes().enumerate() {
                match symbol {
                    b'.' => {}
                    b'(' => openings.push(index),
                    b')' => {
                        let left = openings.pop().expect("balanced structure");
                        let left = sequence[left].to_ascii_uppercase();
                        let right = sequence[index].to_ascii_uppercase();
                        assert!(matches!(
                            (left, right),
                            (b'A' | b'G', b'U') | (b'U', b'A' | b'G') | (b'C', b'G') | (b'G', b'C')
                        ));
                    }
                    symbol => panic!("unexpected dot-bracket symbol {symbol:?}"),
                }
            }
            assert!(openings.is_empty());
        }
    }
}
