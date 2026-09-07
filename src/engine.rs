use std::{cmp::Ordering, collections::BinaryHeap};

use rustc_hash::FxHashMap;

use crate::{
    FoldError, FoldOptions,
    constraints::Constraints,
    model::{Score, ScoringModel},
    sequence::Base,
};

const MAX_SINGLE_LOOP: usize = 30;

type Beam<S> = FxHashMap<usize, State<S>>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Transition {
    None,
    HairpinCandidate,
    Hairpin,
    Single { left: u8, right: u8 },
    Helix,
    Multi { left: u8, right: usize },
    MultiExtended { left: u8, right: usize },
    PairFromMulti,
    M2FromMAndPair { split: usize },
    MFromM2,
    MFromMUnpaired,
    MFromPair,
    CompleteUnpaired,
    CompletePair { split: Option<usize> },
}

#[derive(Debug, Clone, Copy)]
struct State<S: Score> {
    score: S,
    transition: Transition,
}

impl<S: Score> State<S> {
    const fn unreachable() -> Self {
        Self {
            score: S::NEG_INFINITY,
            transition: Transition::None,
        }
    }

    const fn new(score: S, transition: Transition) -> Self {
        Self { score, transition }
    }
}

struct Charts<S: Score> {
    hairpin: Vec<Beam<S>>,
    paired: Vec<Beam<S>>,
    m2: Vec<Beam<S>>,
    multi: Vec<Beam<S>>,
    m: Vec<Beam<S>>,
    complete: Vec<State<S>>,
    sorted_m: Vec<Vec<(S, usize)>>,
}

#[derive(Debug, Clone, Copy)]
struct CubeCandidate<S: Score> {
    heuristic: S,
    paired_index: usize,
    m_index: usize,
}

impl<S: Score> PartialEq for CubeCandidate<S> {
    fn eq(&self, other: &Self) -> bool {
        self.heuristic.total_cmp(other.heuristic) == Ordering::Equal
            && self.paired_index == other.paired_index
            && self.m_index == other.m_index
    }
}

impl<S: Score> Eq for CubeCandidate<S> {}

impl<S: Score> PartialOrd for CubeCandidate<S> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<S: Score> Ord for CubeCandidate<S> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.heuristic
            .total_cmp(other.heuristic)
            .then_with(|| self.paired_index.cmp(&other.paired_index))
            .then_with(|| self.m_index.cmp(&other.m_index))
    }
}

impl<S: Score> Charts<S> {
    fn new(length: usize) -> Self {
        fn beams<S: Score>(length: usize) -> Vec<Beam<S>> {
            std::iter::repeat_with(Beam::default).take(length).collect()
        }
        Self {
            hairpin: beams(length),
            paired: beams(length),
            m2: beams(length),
            multi: beams(length),
            m: beams(length),
            complete: vec![State::unreachable(); length],
            sorted_m: vec![Vec::new(); length],
        }
    }
}

fn update<S: Score>(state: &mut State<S>, score: S, transition: Transition) {
    if state.transition == Transition::None || score > state.score {
        *state = State::new(score, transition);
    }
}

fn update_beam<S: Score>(beam: &mut Beam<S>, index: usize, score: S, transition: Transition) {
    let state = beam.entry(index).or_insert_with(State::unreachable);
    update(state, score, transition);
}

fn add_scores<S: Score>(left: S, right: S) -> S {
    if left == S::NEG_INFINITY || right == S::NEG_INFINITY {
        S::NEG_INFINITY
    } else {
        left + right
    }
}

fn prune<S: Score>(
    beam: &mut Beam<S>,
    complete: &[State<S>],
    beam_size: usize,
    scratch: &mut Vec<(S, usize)>,
) {
    if beam_size == 0 || beam.len() <= beam_size {
        return;
    }
    scratch.clear();
    scratch.extend(beam.iter().map(|(&left, state)| {
        let prefix = left
            .checked_sub(1)
            .map_or(S::ZERO, |index| complete[index].score);
        (add_scores(prefix, state.score), left)
    }));
    let threshold_index = scratch.len() - beam_size;
    let (_, threshold, _) =
        scratch.select_nth_unstable_by(threshold_index, |a, b| a.0.total_cmp(b.0));
    let threshold = threshold.0;
    for &(score, index) in scratch.iter() {
        if score < threshold {
            beam.remove(&index);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prune_removes_scores_below_the_cutoff() {
        let mut beam = Beam::default();
        beam.insert(0, State::new(1, Transition::Hairpin));
        beam.insert(1, State::new(2, Transition::Hairpin));
        let complete = vec![State::new(0, Transition::CompleteUnpaired)];
        let mut scratch = Vec::new();

        prune(&mut beam, &complete, 1, &mut scratch);

        assert_eq!(beam.len(), 1);
        assert!(beam.contains_key(&1));
    }

    #[test]
    fn prune_retains_cutoff_ties_after_prefix_scoring() {
        let mut beam = Beam::default();
        beam.insert(0, State::new(0, Transition::Hairpin));
        beam.insert(1, State::new(0, Transition::Hairpin));
        beam.insert(2, State::new(0, Transition::Hairpin));
        beam.insert(3, State::new(0, Transition::Hairpin));
        let complete = vec![
            State::new(1, Transition::CompleteUnpaired),
            State::new(1, Transition::CompleteUnpaired),
            State::new(2, Transition::CompleteUnpaired),
        ];
        let mut scratch = Vec::new();

        prune(&mut beam, &complete, 2, &mut scratch);

        assert_eq!(beam.len(), 3);
        assert!(!beam.contains_key(&0));
        assert!(beam.contains_key(&1));
        assert!(beam.contains_key(&2));
        assert!(beam.contains_key(&3));
    }

    #[test]
    fn zero_beam_size_disables_pruning() {
        let mut beam = Beam::default();
        beam.insert(0, State::new(1, Transition::Hairpin));
        beam.insert(1, State::new(2, Transition::Hairpin));
        let complete = vec![State::new(0, Transition::CompleteUnpaired)];
        let mut scratch = Vec::new();

        prune(&mut beam, &complete, 0, &mut scratch);

        assert_eq!(beam.len(), 2);
    }
}

fn expansion_entries<S: Score>(beam: &Beam<S>, scratch: &mut Vec<(usize, State<S>)>) {
    scratch.clear();
    scratch.extend(beam.iter().map(|(&index, &state)| (index, state)));
    scratch.sort_unstable_by_key(|(index, _)| std::cmp::Reverse(*index));
}

fn next_pairs(bases: &[Base], constraints: Option<&Constraints>) -> [Vec<Option<usize>>; 5] {
    std::array::from_fn(|kind| {
        let left = match kind {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            3 => Base::U,
            _ => Base::Unknown,
        };
        let mut result = vec![None; bases.len()];
        let mut next = None;
        for index in (0..bases.len()).rev() {
            result[index] = next;
            let can_use = constraints.is_none_or(|value| !value.is_forced_unpaired(index));
            if can_use && left.can_pair(bases[index]) {
                next = Some(index);
            }
        }
        result
    })
}

fn next_pair(table: &[Vec<Option<usize>>; 5], base: Base, after: usize) -> Option<usize> {
    table[base as usize][after]
}

fn can_pair(bases: &[Base], constraints: Option<&Constraints>, left: usize, right: usize) -> bool {
    bases[left].can_pair(bases[right])
        && constraints.is_none_or(|value| value.can_pair(left, right))
}

pub(crate) struct Filled<M: ScoringModel> {
    model: M,
    nucs: Vec<u8>,
    charts: Charts<M::Score>,
}

pub(crate) fn run_linear_from_filled<M: ScoringModel>(
    filled: &Filled<M>,
    length: usize,
) -> Result<(String, M::Score), FoldError> {
    run_linear_from_filled_impl(filled, length)
}

pub(crate) fn run_circular_vienna(
    filled: Filled<crate::model::ViennaModel>,
    constraints: Option<&Constraints>,
    length: usize,
) -> Result<(String, i32), FoldError> {
    run_circular_vienna_impl(filled, constraints, length)
}

#[allow(clippy::too_many_lines)]
pub(crate) fn fill_charts<M: ScoringModel>(
    bases: &[Base],
    constraints: Option<&Constraints>,
    options: &FoldOptions,
) -> Result<Filled<M>, FoldError> {
    let length = bases.len();
    let model = M::new(bases, options.dangles);
    let nucs: Vec<u8> = bases.iter().copied().map(M::encode).collect();
    let next = next_pairs(bases, constraints);
    let mut charts = Charts::new(length);
    let scratch_capacity = options.beam_size.min(length).max(1);
    let mut expansion_scratch = Vec::with_capacity(scratch_capacity);
    let mut prune_scratch = Vec::with_capacity(scratch_capacity);
    if constraints.is_none_or(|value| value.can_unpair(0)) {
        charts.complete[0] = State::new(model.external_unpaired(1), Transition::CompleteUnpaired);
    }
    if length > 1 && constraints.is_none_or(|value| value.can_unpair(0) && value.can_unpair(1)) {
        charts.complete[1] = State::new(model.external_unpaired(2), Transition::CompleteUnpaired);
    }

    for j in 0..length {
        prune(
            &mut charts.hairpin[j],
            &charts.complete,
            options.beam_size,
            &mut prune_scratch,
        );

        let mut first_right = next_pair(&next, bases[j], j);
        if !options.allow_sharp_turns {
            while first_right.is_some_and(|right| right - j < 4) {
                first_right = first_right.and_then(|right| next_pair(&next, bases[j], right));
            }
        }
        if let Some(value) = constraints.and_then(|constraint| constraint.forced_partner(j)) {
            first_right = (value > j).then_some(value);
        }
        if let Some(right) = first_right {
            let within_boundary =
                constraints.is_none_or(|value| right <= value.pair_boundary_after(j));
            if within_boundary && can_pair(bases, constraints, j, right) {
                let score = model.hairpin(j, right, &nucs);
                update_beam(
                    &mut charts.hairpin[right],
                    j,
                    score,
                    Transition::HairpinCandidate,
                );
            }
        }

        expansion_entries(&charts.hairpin[j], &mut expansion_scratch);
        for (i, state) in expansion_scratch.iter().copied() {
            update_beam(&mut charts.paired[j], i, state.score, Transition::Hairpin);
            let Some(right) = next_pair(&next, bases[i], j) else {
                continue;
            };
            let allowed = constraints.is_none_or(|value| {
                right <= value.pair_boundary_after(i) && value.can_pair(i, right)
            });
            if allowed {
                update_beam(
                    &mut charts.hairpin[right],
                    i,
                    model.hairpin(i, right, &nucs),
                    Transition::HairpinCandidate,
                );
            }
        }

        if j == 0 {
            continue;
        }

        prune(
            &mut charts.multi[j],
            &charts.complete,
            options.beam_size,
            &mut prune_scratch,
        );
        expansion_entries(&charts.multi[j], &mut expansion_scratch);
        for (i, state) in expansion_scratch.iter().copied() {
            update_beam(
                &mut charts.paired[j],
                i,
                state.score + model.multi(i, j, &nucs),
                Transition::PairFromMulti,
            );
            let Some(right) = next_pair(&next, bases[i], j) else {
                continue;
            };
            let allowed = constraints.is_none_or(|value| {
                right <= value.pair_boundary_after(j) && value.can_pair(i, right)
            });
            if allowed {
                let (Transition::Multi {
                    left: left_padding,
                    right: right_padding,
                }
                | Transition::MultiExtended {
                    left: left_padding,
                    right: right_padding,
                }) = state.transition
                else {
                    return Err(FoldError::InternalInvariant("invalid multi trace"));
                };
                let extra = right - j;
                let right_padding = right_padding + extra;
                update_beam(
                    &mut charts.multi[right],
                    i,
                    state.score + model.multi_unpaired(extra),
                    Transition::MultiExtended {
                        left: left_padding,
                        right: right_padding,
                    },
                );
            }
        }

        prune(
            &mut charts.paired[j],
            &charts.complete,
            options.beam_size,
            &mut prune_scratch,
        );
        expansion_entries(&charts.paired[j], &mut expansion_scratch);
        let use_cube_pruning = options.beam_size > 20 && expansion_scratch.len() > 20;
        for &(i, state) in &expansion_scratch {
            if i > 0 && j + 6 < length {
                update_beam(
                    &mut charts.m[j],
                    i,
                    state.score + model.m1(i, j, &nucs),
                    Transition::MFromPair,
                );
            }

            match i.checked_sub(1) {
                Some(k)
                    if !use_cube_pruning && j + 1 < length && k > 0 && !charts.m[k].is_empty() =>
                {
                    let branch_score = state.score + model.m1(i, j, &nucs);
                    for (&new_i, left_state) in &charts.m[k] {
                        update_beam(
                            &mut charts.m2[j],
                            new_i,
                            left_state.score + branch_score,
                            Transition::M2FromMAndPair { split: k },
                        );
                    }
                }
                Some(_) | None => {}
            }

            let split = i.checked_sub(1);
            let prefix = split.map_or(M::Score::ZERO, |index| charts.complete[index].score);
            if split.is_none() || prefix != M::Score::NEG_INFINITY {
                // `LinearFold.cpp:1136` sums the external-pair term first. `f64` addition is
                // not associative, so regrouping changes the last ULP and flips tie-breaking
                // between equally optimal structures.
                let score = model.external_paired(i, j, &nucs) + prefix + state.score;
                update(
                    &mut charts.complete[j],
                    score,
                    Transition::CompletePair { split },
                );
            }

            if i > 0 && j + 1 < length {
                let first_p = i.saturating_sub(MAX_SINGLE_LOOP);
                for p in (first_p..i).rev() {
                    if constraints.is_some_and(|value| p + 1 < i && !value.can_unpair(p + 1)) {
                        break;
                    }
                    let mut q = next_pair(&next, bases[p], j);
                    if let Some(forced) = constraints.and_then(|value| value.forced_partner(p)) {
                        if forced < p {
                            break;
                        }
                        q = Some(forced);
                    }
                    while let Some(right) = q {
                        let loop_size = i - p + right - j - 2;
                        if loop_size > MAX_SINGLE_LOOP {
                            break;
                        }
                        let allowed = constraints.is_none_or(|value| {
                            (right <= j + 1 || right <= value.pair_boundary_after(j))
                                && value.can_pair(p, right)
                        });
                        if !allowed {
                            break;
                        }
                        if p + 1 == i && right == j + 1 {
                            let score = state.score + model.helix(p, right, i, j, &nucs);
                            update_beam(&mut charts.paired[right], p, score, Transition::Helix);
                        } else {
                            let score = state.score + model.single(p, right, i, j, &nucs);
                            let left = u8::try_from(i - p).map_err(|_| {
                                FoldError::InternalInvariant("single left padding overflow")
                            })?;
                            let right_padding = u8::try_from(right - j).map_err(|_| {
                                FoldError::InternalInvariant("single right padding overflow")
                            })?;
                            update_beam(
                                &mut charts.paired[right],
                                p,
                                score,
                                Transition::Single {
                                    left,
                                    right: right_padding,
                                },
                            );
                        }
                        q = next_pair(&next, bases[p], right);
                    }
                }
            }
        }

        if use_cube_pruning && j + 1 < length {
            let mut valid_pairs = Vec::new();
            for &(i, state) in &expansion_scratch {
                let Some(k) = i.checked_sub(1) else {
                    continue;
                };
                if k == 0 || charts.sorted_m[k].is_empty() {
                    continue;
                }
                let branch_score = state.score + model.m1(i, j, &nucs);
                let is_candidate = charts.m2[j]
                    .get(&i)
                    .is_none_or(|existing| branch_score > existing.score);
                if is_candidate {
                    valid_pairs.push((i, k, branch_score));
                }
            }

            let mut heap = BinaryHeap::new();
            for (paired_index, &(_, k, branch_score)) in valid_pairs.iter().enumerate() {
                heap.push(CubeCandidate {
                    heuristic: branch_score + charts.sorted_m[k][0].0,
                    paired_index,
                    m_index: 0,
                });
            }

            let mut filled = 0;
            let mut previous = M::Score::NEG_INFINITY;
            let mut current = M::Score::NEG_INFINITY;
            while (filled < options.beam_size || current == previous) && !heap.is_empty() {
                let Some(candidate) = heap.pop() else {
                    break;
                };
                previous = current;
                current = candidate.heuristic;
                let (_, k, branch_score) = valid_pairs[candidate.paired_index];
                let new_i = charts.sorted_m[k][candidate.m_index].1;
                let new_score = branch_score + charts.m[k][&new_i].score;
                if !charts.m2[j].contains_key(&new_i) {
                    filled += 1;
                    update_beam(
                        &mut charts.m2[j],
                        new_i,
                        new_score,
                        Transition::M2FromMAndPair { split: k },
                    );
                }

                let mut next_index = candidate.m_index + 1;
                while next_index < charts.sorted_m[k].len() {
                    let candidate_i = charts.sorted_m[k][next_index].1;
                    if !charts.m2[j].contains_key(&candidate_i) {
                        heap.push(CubeCandidate {
                            heuristic: branch_score + charts.sorted_m[k][next_index].0,
                            paired_index: candidate.paired_index,
                            m_index: next_index,
                        });
                        break;
                    }
                    next_index += 1;
                }
            }
        }

        prune(
            &mut charts.m2[j],
            &charts.complete,
            options.beam_size,
            &mut prune_scratch,
        );
        expansion_entries(&charts.m2[j], &mut expansion_scratch);
        for (i, state) in expansion_scratch.iter().copied() {
            update_beam(&mut charts.m[j], i, state.score, Transition::MFromM2);
            let first_p = i.saturating_sub(MAX_SINGLE_LOOP);
            for p in (first_p..i).rev() {
                if constraints.is_some_and(|value| p + 1 < i && !value.can_unpair(p + 1)) {
                    break;
                }
                let mut q = next_pair(&next, bases[p], j);
                if let Some(forced) = constraints.and_then(|value| value.forced_partner(p)) {
                    if forced < p {
                        break;
                    }
                    q = Some(forced);
                }
                let Some(right) = q else {
                    continue;
                };
                let allowed = constraints.is_none_or(|value| {
                    (right <= j + 1 || right <= value.pair_boundary_after(j))
                        && value.can_pair(p, right)
                });
                if !allowed {
                    continue;
                }
                let left_padding = i - p;
                let right_padding = right - j;
                let score = (model.multi_unpaired(left_padding - 1)
                    + model.multi_unpaired(right_padding - 1))
                    + state.score;
                update_beam(
                    &mut charts.multi[right],
                    p,
                    score,
                    Transition::Multi {
                        left: u8::try_from(left_padding).map_err(|_| {
                            FoldError::InternalInvariant("multi left padding overflow")
                        })?,
                        right: right_padding,
                    },
                );
            }
        }

        prune(
            &mut charts.m[j],
            &charts.complete,
            options.beam_size,
            &mut prune_scratch,
        );
        charts.sorted_m[j] = charts.m[j]
            .iter()
            .map(|(&left, state)| {
                let prefix = left.checked_sub(1).map_or(M::Score::ZERO, |index| {
                    let score = charts.complete[index].score;
                    if constraints.is_some() && score == M::Score::NEG_INFINITY {
                        M::Score::ZERO
                    } else {
                        score
                    }
                });
                (add_scores(prefix, state.score), left)
            })
            .collect();
        charts.sorted_m[j].sort_unstable_by(|left, right| {
            right.0.total_cmp(left.0).then_with(|| right.1.cmp(&left.1))
        });
        expansion_entries(&charts.m[j], &mut expansion_scratch);
        if j + 1 < length {
            for (i, state) in expansion_scratch.iter().copied() {
                if constraints.is_some_and(|value| !value.can_unpair(j + 1)) {
                    continue;
                }
                update_beam(
                    &mut charts.m[j + 1],
                    i,
                    state.score + model.multi_unpaired(1),
                    Transition::MFromMUnpaired,
                );
            }

            if constraints.is_none_or(|value| value.can_unpair(j + 1))
                && charts.complete[j].transition != Transition::None
            {
                let score = charts.complete[j].score + model.external_unpaired(1);
                update(
                    &mut charts.complete[j + 1],
                    score,
                    Transition::CompleteUnpaired,
                );
            }
        }
    }

    Ok(Filled {
        model,
        nucs,
        charts,
    })
}

pub(crate) fn run_linear<M: ScoringModel>(
    bases: &[Base],
    constraints: Option<&Constraints>,
    options: &FoldOptions,
) -> Result<(String, M::Score), FoldError> {
    let filled = fill_charts::<M>(bases, constraints, options)?;
    run_linear_from_filled(&filled, bases.len())
}

fn run_linear_from_filled_impl<M: ScoringModel>(
    filled: &Filled<M>,
    length: usize,
) -> Result<(String, M::Score), FoldError> {
    let final_state = filled.charts.complete[length - 1];
    if final_state.transition == Transition::None {
        return Err(FoldError::NoValidStructure);
    }
    let structure = traceback(&filled.charts, length)?;
    Ok((structure, final_state.score))
}

#[derive(Debug, Clone, Copy)]
struct RootPair {
    left: usize,
    right: usize,
    score: i32,
}

struct PairInventory {
    pairs: Vec<RootPair>,
    offsets: Vec<usize>,
}

impl PairInventory {
    fn build(charts: &Charts<i32>, length: usize) -> Self {
        let mut counts = vec![0usize; length];
        for beam in &charts.paired {
            for &left in beam.keys() {
                counts[left] += 1;
            }
        }
        let mut offsets = vec![0usize; length + 1];
        for index in 0..length {
            offsets[index + 1] = offsets[index] + counts[index];
        }
        let mut cursors = offsets[..length].to_vec();
        let mut pairs = vec![
            RootPair {
                left: 0,
                right: 0,
                score: 0
            };
            offsets[length]
        ];
        for (right, beam) in charts.paired.iter().enumerate() {
            for (&left, state) in beam {
                let position = cursors[left];
                pairs[position] = RootPair {
                    left,
                    right,
                    score: state.score,
                };
                cursors[left] += 1;
            }
        }
        Self { pairs, offsets }
    }

    fn starting_at(&self, left: usize) -> &[RootPair] {
        &self.pairs[self.offsets[left]..self.offsets[left + 1]]
    }
}

fn run_circular_vienna_impl(
    filled: Filled<crate::model::ViennaModel>,
    constraints: Option<&Constraints>,
    length: usize,
) -> Result<(String, i32), FoldError> {
    let inventory = PairInventory::build(&filled.charts, length);
    let mut best_score = if constraints.is_none_or(|value| value.range_can_unpair(0, length)) {
        crate::model::ViennaModel::circular_unpaired_score(length)
    } else {
        i32::MIN
    };
    let mut best_roots = Vec::new();

    for pair in &inventory.pairs {
        let outside_left = constraints.is_none_or(|value| value.range_can_unpair(0, pair.left));
        let outside_right =
            constraints.is_none_or(|value| value.range_can_unpair(pair.right + 1, length));
        if !outside_left || !outside_right {
            continue;
        }
        let loop_score = filled
            .model
            .circular_hairpin_score(pair.left, pair.right, &filled.nucs);
        if loop_score != i32::MIN && pair.score + loop_score > best_score {
            best_score = pair.score + loop_score;
            best_roots.clear();
            best_roots.push((pair.left, pair.right));
        }
    }

    for first in &inventory.pairs {
        if first.left > MAX_SINGLE_LOOP {
            continue;
        }
        let remaining = MAX_SINGLE_LOOP - first.left;
        for a in 0..=remaining {
            for b in 0..=remaining - a {
                let p = first.right + 1 + a;
                let Some(q) = length.checked_sub(1 + b) else {
                    continue;
                };
                if p >= q {
                    continue;
                }
                let allowed = constraints.is_none_or(|value| {
                    value.range_can_unpair(first.right + 1, p)
                        && value.range_can_unpair(q + 1, length)
                        && value.range_can_unpair(0, first.left)
                });
                if !allowed {
                    continue;
                }
                let Some(second) = filled.charts.paired[q].get(&p) else {
                    continue;
                };
                let loop_score = crate::model::ViennaModel::circular_internal_score(
                    first.left,
                    first.right,
                    p,
                    q,
                    &filled.nucs,
                );
                let score = first.score + second.score + loop_score;
                if score > best_score {
                    best_score = score;
                    best_roots.clear();
                    best_roots.push((first.left, first.right));
                    best_roots.push((p, q));
                }
            }
        }
    }

    #[derive(Clone, Copy)]
    enum RootStep {
        Unreachable,
        Start,
        Unpaired {
            previous: usize,
            class: usize,
        },
        Branch {
            previous: usize,
            class: usize,
            left: usize,
            right: usize,
        },
    }

    let mut prefix = vec![i32::MIN; (length + 1) * 4];
    let mut steps = vec![RootStep::Unreachable; (length + 1) * 4];
    prefix[0] = 0;
    steps[0] = RootStep::Start;
    for boundary in 0..length {
        for class in 0..4 {
            let current = prefix[boundary * 4 + class];
            if current == i32::MIN {
                continue;
            }
            for pair in inventory.starting_at(boundary) {
                let next_class = (class + 1).min(3);
                let target = (pair.right + 1) * 4 + next_class;
                let score = current
                    + pair.score
                    + filled.model.circular_multiloop_stem_score(
                        pair.left,
                        pair.right,
                        &filled.nucs,
                    );
                if score > prefix[target] {
                    prefix[target] = score;
                    steps[target] = RootStep::Branch {
                        previous: boundary,
                        class,
                        left: pair.left,
                        right: pair.right,
                    };
                }
            }
            if constraints.is_none_or(|value| value.can_unpair(boundary)) {
                let target = (boundary + 1) * 4 + class;
                if current > prefix[target] {
                    prefix[target] = current;
                    steps[target] = RootStep::Unpaired {
                        previous: boundary,
                        class,
                    };
                }
            }
        }
    }

    let multi_index = length * 4 + 3;
    if prefix[multi_index] != i32::MIN {
        let score =
            prefix[multi_index] + crate::model::ViennaModel::circular_multiloop_closing_score();
        if score > best_score {
            best_score = score;
            best_roots.clear();
            let mut boundary = length;
            let mut class = 3;
            while boundary > 0 || class > 0 {
                let step = steps[boundary * 4 + class];
                match step {
                    RootStep::Branch {
                        previous,
                        class: previous_class,
                        left,
                        right,
                    } => {
                        best_roots.push((left, right));
                        boundary = previous;
                        class = previous_class;
                    }
                    RootStep::Unpaired {
                        previous,
                        class: previous_class,
                    } => {
                        boundary = previous;
                        class = previous_class;
                    }
                    RootStep::Start => break,
                    RootStep::Unreachable => {
                        return Err(FoldError::InternalInvariant("missing circular root trace"));
                    }
                }
            }
            best_roots.reverse();
        }
    }

    if best_score == i32::MIN {
        return Err(FoldError::NoValidStructure);
    }
    if best_roots.is_empty() {
        return Ok((".".repeat(length), best_score));
    }
    let structure = traceback_from_paired_roots(&filled.charts, length, &best_roots)?;
    Ok((structure, best_score))
}

fn state_from_beam<S: Score>(beam: &[Beam<S>], i: usize, j: usize) -> Result<State<S>, FoldError> {
    beam[j]
        .get(&i)
        .copied()
        .ok_or(FoldError::InternalInvariant("trace state is missing"))
}

#[derive(Clone, Copy)]
enum TraceKind {
    Complete,
    Paired,
    M,
    M2,
    Multi,
}

#[allow(clippy::too_many_lines)]
fn traceback<S: Score>(charts: &Charts<S>, length: usize) -> Result<String, FoldError> {
    traceback_with_stack(
        charts,
        length,
        vec![(
            TraceKind::Complete,
            0,
            length - 1,
            charts.complete[length - 1],
        )],
    )
}

fn traceback_from_paired_roots<S: Score>(
    charts: &Charts<S>,
    length: usize,
    roots: &[(usize, usize)],
) -> Result<String, FoldError> {
    let mut stack = Vec::with_capacity(roots.len());
    for &(left, right) in roots.iter().rev() {
        stack.push((
            TraceKind::Paired,
            left,
            right,
            state_from_beam(&charts.paired, left, right)?,
        ));
    }
    traceback_with_stack(charts, length, stack)
}

fn traceback_with_stack<S: Score>(
    charts: &Charts<S>,
    length: usize,
    mut stack: Vec<(TraceKind, usize, usize, State<S>)>,
) -> Result<String, FoldError> {
    let mut structure = vec![b'.'; length];

    while let Some((kind, i, j, state)) = stack.pop() {
        match state.transition {
            Transition::HairpinCandidate => {}
            Transition::Hairpin => {
                structure[i] = b'(';
                structure[j] = b')';
            }
            Transition::Single { left, right } => {
                structure[i] = b'(';
                structure[j] = b')';
                let p = i + usize::from(left);
                let q = j - usize::from(right);
                stack.push((
                    TraceKind::Paired,
                    p,
                    q,
                    state_from_beam(&charts.paired, p, q)?,
                ));
            }
            Transition::Helix => {
                structure[i] = b'(';
                structure[j] = b')';
                stack.push((
                    TraceKind::Paired,
                    i + 1,
                    j - 1,
                    state_from_beam(&charts.paired, i + 1, j - 1)?,
                ));
            }
            Transition::Multi { left, right } | Transition::MultiExtended { left, right } => {
                let p = i + usize::from(left);
                let q = j - right;
                stack.push((TraceKind::M2, p, q, state_from_beam(&charts.m2, p, q)?));
            }
            Transition::PairFromMulti => {
                structure[i] = b'(';
                structure[j] = b')';
                stack.push((
                    TraceKind::Multi,
                    i,
                    j,
                    state_from_beam(&charts.multi, i, j)?,
                ));
            }
            Transition::M2FromMAndPair { split } => {
                stack.push((
                    TraceKind::M,
                    i,
                    split,
                    state_from_beam(&charts.m, i, split)?,
                ));
                stack.push((
                    TraceKind::Paired,
                    split + 1,
                    j,
                    state_from_beam(&charts.paired, split + 1, j)?,
                ));
            }
            Transition::MFromM2 => {
                stack.push((TraceKind::M2, i, j, state_from_beam(&charts.m2, i, j)?));
            }
            Transition::MFromMUnpaired => stack.push((
                TraceKind::M,
                i,
                j - 1,
                state_from_beam(&charts.m, i, j - 1)?,
            )),
            Transition::MFromPair => stack.push((
                TraceKind::Paired,
                i,
                j,
                state_from_beam(&charts.paired, i, j)?,
            )),
            Transition::CompleteUnpaired => {
                if j > 0 {
                    stack.push((TraceKind::Complete, 0, j - 1, charts.complete[j - 1]));
                }
            }
            Transition::CompletePair { split } => {
                let left = split.map_or(0, |value| value + 1);
                stack.push((
                    TraceKind::Paired,
                    left,
                    j,
                    state_from_beam(&charts.paired, left, j)?,
                ));
                if let Some(split) = split {
                    stack.push((TraceKind::Complete, 0, split, charts.complete[split]));
                }
            }
            Transition::None => {
                return Err(FoldError::InternalInvariant(match kind {
                    TraceKind::Complete => "missing complete trace",
                    TraceKind::Paired => "missing paired trace",
                    TraceKind::M => "missing M trace",
                    TraceKind::M2 => "missing M2 trace",
                    TraceKind::Multi => "missing multi trace",
                }));
            }
        }
    }

    String::from_utf8(structure).map_err(|_| FoldError::InternalInvariant("invalid traceback"))
}
