use std::{
    io::{BufWriter, Write},
    path::Path,
};

use clap::{Parser, ValueEnum};
use itertools::Itertools;
use pa_types::{Base, Seq, Sequence};
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

/// Each `ErrorModel` creates a different type of sequence pair.
#[derive(ValueEnum, Default, Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum ErrorModel {
    /// Generate two independent sequences of length n.
    Independent,

    /// Make e*n random substitution/insertion/deletion, each with equal probability.
    /// TODO(ragnar): Support non-equal sub/ins/del probabilities?
    #[default]
    Uniform,

    /// Delete a region of size e*n.
    Delete,
    /// Insert a region of size e*n.
    Insert,
    /// Move a region of size e*n to a random position.
    Move,
    /// Duplicate a region of size e*n.
    Duplicate,

    /// Make an en/2 insertion and apply en/2 noise.
    NoisyInsert,
    /// Make an en/2 deletion and apply en/2 noise.
    NoisyDelete,
    /// Make an en/2 move and apply en/2 noise.
    NoisyMove,
    /// Duplicate a region of size en/2, and apply en/2 noise.
    NoisyDuplicate,

    /// Repeat a pattern of given `pattern_length`, and apply e*n mutations.
    Repeat,
    /// Repeat a pattern of given `pattern_length`, and apply e*n mutations to
    /// get A. Then apply e*n more mutations to get B.
    NestedRepeat,
    /// Repeat a pattern of given `pattern_length` to get sequences of length
    /// `n` and `n-en/2`. Then apply en/2 mutations to each.
    DifferentRepeat,
    /// Repeat a pattern of given `pattern_length`. Apply en/2 mutations to get each of A and B.
    SymmetricRepeat,
}

/// Options to generate a single pair of sequences.
#[derive(Parser, Clone, Serialize, Deserialize, Debug)]
pub struct GenerateOptions {
    /// Target length of each generated sequence.
    #[arg(short = 'n', long)]
    pub length: usize,

    /// Error rate between sequences.
    #[arg(short, long)]
    pub error_rate: f32,

    /// The type of error to generate.
    #[arg(long, value_enum, default_value_t, value_name = "MODEL")]
    pub error_model: ErrorModel,

    /// The length of a repeating pattern.
    ///
    /// Must be positive. Default: uniform in [1, sqrt(n)]
    #[arg(long, hide_short_help = true)]
    pub pattern_length: Option<usize>,
}

/// Options to generate multiple sequence pairs.
#[derive(Parser, Clone, Serialize, Deserialize, Debug)]
#[command(author, version, about)]
#[clap(group(
    clap::ArgGroup::new("total_size")
        .multiple(false)
        .args(&["cnt", "size"]),
))]
pub struct GenerateArgs {
    /// Options for generating a single pair.
    #[command(flatten)]
    pub options: GenerateOptions,

    /// The number of sequence pairs to generate.
    ///
    /// Conflicts with --size.
    #[arg(short = 'x', long, default_value = "1")]
    pub cnt: Option<usize>,

    /// The total input size to generate.
    ///
    /// Generate as many pairs as needed to get to a total number of bases as close
    /// as possible to 2*size.
    #[arg(short = 's', long)]
    pub size: Option<usize>,

    /// RNG seed. Randomized if not set.
    #[arg(long)]
    pub seed: Option<u64>,
}

/// List of characters that can be generated.
const ALPH: [Base; 4] = [b'A', b'C', b'G', b'T'];

fn rand_base(rng: &mut impl Rng) -> Base {
    ALPH[rng.gen_range(0..ALPH.len())] as Base
}

/// Random mutations that can be generated.
enum Mutation {
    // Replace char at pos.
    Substitution(usize, u8),
    // Insert char before pos.
    Insertion(usize, u8),
    // Delete char at pos.
    Deletion(usize),
}

fn random_mutation(len: usize, rng: &mut impl Rng) -> Mutation {
    // Substitution / insertion / deletion all with equal probability.
    // For length 0 sequences, only generate insertions.
    match if len == 0 { 1 } else { rng.gen_range(0..3) } {
        0 => Mutation::Substitution(rng.gen_range(0..len), rand_base(rng)),
        1 => Mutation::Insertion(rng.gen_range(0..len + 1), rand_base(rng)),
        2 => Mutation::Deletion(rng.gen_range(0..len)),
        _ => unreachable!(),
    }
}

fn mutate_once(seq: &mut ropey::Rope, rng: &mut impl Rng) {
    let m = random_mutation(seq.len_bytes(), rng);
    match m {
        Mutation::Substitution(i, c) => {
            seq.remove(i..=i);
            seq.insert(i, std::str::from_utf8(&[c]).unwrap());
        }
        Mutation::Insertion(i, c) => seq.insert(i, std::str::from_utf8(&[c]).unwrap()),
        Mutation::Deletion(i) => {
            seq.remove(i..=i);
        }
    }
}

fn mutate(seq: Seq, mutations: usize, rng: &mut impl Rng) -> Sequence {
    if mutations == 0 {
        return seq.to_vec();
    }
    let mut seq = ropey::Rope::from_str(std::str::from_utf8(seq).unwrap());
    for _ in 0..mutations {
        mutate_once(&mut seq, rng);
    }
    seq.to_string().into_bytes()
}

/// Generate a uniform random sequence.
pub fn random_sequence(len: usize, rng: &mut impl Rng) -> Sequence {
    (0..len).map(|_| rand_base(rng)).collect()
}

impl GenerateOptions {
    /// Generate a single sequence pair via the given error model.
    pub fn generate(&self, rng: &mut impl Rng) -> (Sequence, Sequence) {
        use ErrorModel::*;

        let a: Sequence = random_sequence(self.length, rng);
        let muts = (self.error_rate * self.length as f32).ceil() as usize;
        let (indel_len, noise) = match self.error_model {
            Independent => (0, 0),
            Uniform => (0, muts),
            Delete | Insert | Move | Duplicate => (muts, 0),
            NoisyInsert | NoisyDelete | NoisyMove | NoisyDuplicate => (muts / 2, muts / 2),
            Repeat | DifferentRepeat | NestedRepeat | SymmetricRepeat => (0, muts),
        };

        match self.error_model {
            Independent => (a, random_sequence(self.length, rng)),
            Uniform => (a.clone(), mutate(&a, noise, rng)),
            Delete | NoisyDelete => {
                let mut b = a.clone();
                let start = rng.gen_range(0..=b.len() - indel_len);
                b.drain(start..start + indel_len);
                (a, mutate(&b, noise, rng))
            }
            Insert | NoisyInsert => {
                let mut b = a.clone();
                let start = rng.gen_range(0..=b.len());
                b.splice(start..start, random_sequence(indel_len, rng));
                (a, mutate(&b, noise, rng))
            }
            Move | NoisyMove => {
                let mut b = a.clone();
                // deletion
                let start = rng.gen_range(0..=b.len() - indel_len);
                let piece = b.drain(start..start + indel_len).collect_vec();
                // insertion
                let start = rng.gen_range(0..=b.len());
                b.splice(start..start, piece);
                (a, mutate(&b, noise, rng))
            }
            Duplicate | NoisyDuplicate => {
                let mut b = a.clone();
                let start = rng.gen_range(0..=b.len() - indel_len);
                let piece = b.drain(start..start + indel_len).collect_vec();
                b.splice(start..start, piece);
                (a, mutate(&b, noise, rng))
            }
            Repeat | NestedRepeat | DifferentRepeat | SymmetricRepeat => {
                let (pre_noise, b_crop, a_noise, b_noise) = match self.error_model {
                    Repeat => (0, 0, 0, noise),
                    NestedRepeat => (noise / 2, 0, 0, noise / 2),
                    DifferentRepeat => (0, noise / 2, 0, noise / 2),
                    SymmetricRepeat => (0, 0, noise / 2, noise / 2),
                    _ => unreachable!(),
                };

                let len = if let Some(pattern_length) = self.pattern_length {
                    pattern_length
                } else {
                    rng.gen_range(1..=(self.length as f32).sqrt() as usize)
                };
                let pattern = random_sequence(len, rng);
                let mut base = pattern
                    .iter()
                    .cycle()
                    .take(self.length)
                    .copied()
                    .collect_vec();
                base = mutate(&base, pre_noise, rng);
                let a = mutate(&base, a_noise, rng);
                let b = mutate(&base[0..base.len() - b_crop], b_noise, rng);
                (a, b)
            }
        }
    }

    // Convenience functions.

    /// Generate a random pair.
    pub fn random(&self) -> (Sequence, Sequence) {
        self.generate(&mut get_rng(None))
    }

    /// Generate a seeded random pair.
    pub fn seeded(&self, seed: u64) -> (Sequence, Sequence) {
        self.generate(&mut get_rng(Some(seed)))
    }
}

fn get_rng(seed: Option<u64>) -> rand_chacha::ChaCha8Rng {
    match seed {
        Some(seed) => rand_chacha::ChaCha8Rng::seed_from_u64(seed as u64),
        None => rand_chacha::ChaCha8Rng::from_entropy(),
    }
}

impl GenerateArgs {
    pub fn generate(&self) -> Vec<(Sequence, Sequence)> {
        let mut pairs = vec![];
        let mut prev_len = 0;
        let mut total_len = 0;

        let rng = &mut get_rng(self.seed);

        loop {
            // Stop when either >= 2*self.size bp have been generated, or when
            // self.cnt pairs have been generated.
            if let Some(size) = self.size {
                if total_len >= 2 * size {
                    break;
                }
            } else {
                if pairs.len() == self.cnt.unwrap() {
                    break;
                }
            }

            let (a, b) = self.options.generate(rng);
            prev_len = total_len;
            total_len += a.len() + b.len();
            pairs.push((a, b));
        }

        // Remove the last pair if that gives gets closer to the target size.
        if let Some(size) = self.size {
            if prev_len > 0 && 2 * size - prev_len < total_len - 2 * size {
                pairs.pop();
            }
        }

        pairs
    }

    pub fn generate_file(&self, path: &Path) {
        assert_eq!(
            path.extension().unwrap_or_default(),
            "seq",
            "Output file must have .seq extension!"
        );
        let mut f = BufWriter::new(std::fs::File::create(path).unwrap());
        for (a, b) in self.generate() {
            write!(f, ">").unwrap();
            f.write_all(&a).unwrap();
            writeln!(f).unwrap();
            write!(f, "<").unwrap();
            f.write_all(&b).unwrap();
            writeln!(f).unwrap();
        }
    }
}

// Convenience functions.

/// Generate a random pair with length n and error rate e.
pub fn uniform_random(n: usize, e: f32) -> (Sequence, Sequence) {
    GenerateOptions {
        length: n,
        error_rate: e,
        error_model: ErrorModel::Uniform,
        pattern_length: None,
    }
    .random()
}

/// Generate a seeded random pair with length n and error rate e.
pub fn uniform_seeded(n: usize, e: f32, seed: u64) -> (Sequence, Sequence) {
    GenerateOptions {
        length: n,
        error_rate: e,
        error_model: ErrorModel::Uniform,
        pattern_length: None,
    }
    .seeded(seed)
}
