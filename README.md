# linearfold-rs

`linearfold-rs` is a Rust port of the [LinearFold](https://github.com/LinearFold/LinearFold) tool for RNA secondary structure prediction. It supports both the CONTRAfold and Vienna RNAfold models while adding new features and performance improvements over the original implementation.

## API

The crate exposes a single entry point, `fold`:

```rust,ignore
pub fn fold(
    sequence: &[u8],
    constraints: Option<&[u8]>,
    options: &FoldOptions,
) -> Result<Prediction, FoldError>
```

### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `sequence` | `&[u8]` | RNA or DNA sequence. Accepts ASCII letters (uppercased, `T` mapped to `U`, alphabetic ambiguity codes treated as unpairable unknown bases). Empty sequences, non-letter bytes, and whitespace are rejected. |
| `constraints` | `Option<&[u8]>` | Optional dot-bracket constraints of the same length as `sequence`. `?` permits either state, `.` forces an unpaired base, and matching `(` / `)` force a base pair. `None` disables constraints. |
| `options` | `&FoldOptions` | Folding configuration (scoring model, beam size, sharp turns, dangling ends, and sequence topology). |

`FoldOptions` is a struct that customizes the scoring model and prediction parameters used for folding:

| Field | Type | Default | Description |
| --- | --- | --- | --- |
| `model` | `Model` | `Model::ContraFold` | Scoring model (`Model::ContraFold` or `Model::ViennaRnafold`). |
| `beam_size` | `usize` | `100` | Maximum states retained per beam. Zero disables beam pruning. |
| `allow_sharp_turns` | `bool` | `false` | Allow base pairs separated by fewer than three unpaired bases. |
| `dangles` | `DangleModel` | `DangleModel::Both` | Vienna dangling-end treatment (`DangleModel::Both` or `DangleModel::None`). Ignored when using `Model::ContraFold`. |
| `topology` | `SequenceTopology` | `SequenceTopology::Linear` | Sequence topology (`SequenceTopology::Linear` or `SequenceTopology::Circular`). Supported only with `Model::ViennaRnafold`. |

### Example

```rs
use linearfold_rs::{FoldOptions, Model, SequenceTopology, fold};

let options = FoldOptions {
    model: Model::ViennaRnafold,
    topology: SequenceTopology::Linear,
    ..FoldOptions::default()
};
let prediction = fold(b"AAAACGGUCCUUAUCAGGACCAAACA", None, &options)?;
assert_eq!(prediction.structure, ".....((((((....)))))).....");
```

## Differences from the original implementation

- `linearfold-rs` adds support for folding circular RNAs via `SequenceTopology::Circular` in `FoldOptions`. This algorithm was ported from [ViennaRNA](https://github.com/ViennaRNA/ViennaRNA)'s `RNAfold` algorithm and is currently supported only when using `Model::ViennaRnafold`.

## Citation

If you use `linearfold-rs` in your work, please cite the original LinearFold paper:

> Huang, L., Zhang, H., Deng, D., Zhao, K., Liu, K., Hendrix, D. A., & Mathews, D. H. [**"LinearFold: linear-time approximate RNA folding by 5'-to-3'dynamic programming and beam search"**](https://doi.org/10.1093/bioinformatics/btz375). *Bioinformatics* 35.14 (2019): i295-i304.
