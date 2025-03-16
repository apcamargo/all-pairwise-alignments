use block_aligner::{cigar::*, scan_block::*, scores::*};
use itertools::Itertools;
use needletail::parse_fastx_file;
use rayon::iter::{ParallelBridge, ParallelIterator};
use std::env;
use std::process;

const MIN_BLOCK_SIZE: usize = 32;
const MAX_BLOCK_SIZE: usize = 256;
const GAP_OPEN: i8 = -11;
const GAP_EXTEND: i8 = -1;

fn align_sequences(seq1: &[u8], seq2: &[u8]) -> Cigar {
    let gaps = Gaps {
        open: GAP_OPEN,
        extend: GAP_EXTEND,
    };

    // Note that PaddedBytes, Block, and Cigar can be initialized with sequence length
    // and block size upper bounds and be reused later for shorter sequences, to avoid
    // repeated allocations.
    let r = PaddedBytes::from_bytes::<AAMatrix>(seq1, MAX_BLOCK_SIZE);
    let q = PaddedBytes::from_bytes::<AAMatrix>(seq2, MAX_BLOCK_SIZE);

    // Align with traceback, but no X-drop threshold (global alignment).
    let mut a = Block::<true, false>::new(q.len(), r.len(), MAX_BLOCK_SIZE);
    a.align(&q, &r, &BLOSUM62, gaps, MIN_BLOCK_SIZE..=MAX_BLOCK_SIZE, 0);
    let res = a.res();

    // println!("{:?}", res);

    let mut cigar = Cigar::new(res.query_idx, res.reference_idx);
    // Compute traceback and resolve =/X (matches/mismatches).
    let t = a.trace();
    t.cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);

    cigar
}

fn apply_alignment_to_sequences(seq1: &[u8], seq2: &[u8], cigar: &Cigar) -> (Vec<u8>, Vec<u8>) {
    let mut aligned_seq1 = Vec::new();
    let mut aligned_seq2 = Vec::new();

    let mut pos1 = 0;
    let mut pos2 = 0;

    for i in 0..cigar.len() {
        let operation = cigar.get(i);
        match operation.op {
            Operation::Eq | Operation::X | Operation::M => {
                // Match, mismatch, or generic alignment
                // They all have identical processing
                for _ in 0..operation.len {
                    if pos1 < seq1.len() && pos2 < seq2.len() {
                        aligned_seq1.push(seq1[pos1]);
                        aligned_seq2.push(seq2[pos2]);
                        pos1 += 1;
                        pos2 += 1;
                    }
                }
            }
            Operation::I => {
                // Insertion in the first sequence, gap in the second
                for _ in 0..operation.len {
                    if pos1 < seq1.len() {
                        aligned_seq1.push(seq1[pos1]);
                        aligned_seq2.push(b'-'); // Gap character
                        pos1 += 1;
                    }
                }
            }
            Operation::D => {
                // Deletion in the first sequence (insertion in the second)
                for _ in 0..operation.len {
                    if pos2 < seq2.len() {
                        aligned_seq1.push(b'-'); // Gap character
                        aligned_seq2.push(seq2[pos2]);
                        pos2 += 1;
                    }
                }
            }
            Operation::Sentinel => {
                // Sentinel is a marker operation, not actually part of the alignment
                // We can safely ignore it
            }
        }
    }
    (aligned_seq1, aligned_seq2)
}

fn calculate_identity(cigar: &Cigar) -> f64 {
    let mut matching_positions = 0;
    let mut total_alignment_length = 0;

    for i in 0..cigar.len() {
        let operation = cigar.get(i);
        match operation.op {
            // Exact matches contribute to both matching positions and total length
            Operation::Eq => {
                matching_positions += operation.len;
                total_alignment_length += operation.len;
            }
            // Mismatches, generic alignments, insertions and deletions
            // only contribute to the total alignment length
            Operation::X | Operation::M | Operation::I | Operation::D => {
                total_alignment_length += operation.len;
            }
            // Sentinel is ignored
            Operation::Sentinel => {}
        }
    }
    if total_alignment_length == 0 {
        return 0.0;
    }
    matching_positions as f64 / total_alignment_length as f64
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename>", args[0]);
        process::exit(1);
    }

    let filename = &args[1];
    let mut reader = parse_fastx_file(&filename).expect("valid path/file");

    let mut record_list: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    while let Some(Ok(record)) = reader.next() {
        record_list.push((record.id().to_vec(), record.seq().to_vec()));
    }

    record_list
        .iter()
        .combinations(2)
        .par_bridge()
        .for_each(|pair| {
            let cigar = align_sequences(&pair[0].1, &pair[1].1);
            let identity = calculate_identity(&cigar);
            let (aseq1, aseq2) = apply_alignment_to_sequences(&pair[0].1, &pair[1].1, &cigar);
            println!(
                "{}\t{}\t{}\t{:.6}\t{}\t{}",
                String::from_utf8_lossy(&pair[0].0),
                String::from_utf8_lossy(&pair[1].0),
                cigar,
                identity,
                String::from_utf8_lossy(&aseq1),
                String::from_utf8_lossy(&aseq2)
            );
        });
}
