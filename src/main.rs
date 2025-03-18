mod blosum30;
mod blosum45;

use crate::blosum30::blosum30;
use crate::blosum45::blosum45;

use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation::*;
use bio::alignment::pairwise::*;
use bio::scores::blosum62;

use itertools::Itertools;
use needletail::parse_fastx_file;
use rayon::iter::{ParallelBridge, ParallelIterator};
use std::env;
use std::process;

const GAP_OPEN: i32 = -14;
const GAP_EXTEND: i32 = -2;

enum AlignmentType {
    Global,
    Local,
}

fn align_sequences(seq1: &[u8], seq2: &[u8], alignment_type: &AlignmentType) -> Alignment {
    let mut aligner =
        Aligner::with_capacity(seq1.len(), seq2.len(), GAP_OPEN, GAP_EXTEND, &blosum45);
    match alignment_type {
        AlignmentType::Global => aligner.global(seq1, seq2),
        AlignmentType::Local => aligner.local(seq1, seq2),
    }
}

fn apply_alignment_to_sequences(
    seq1: &[u8],
    seq2: &[u8],
    alignment: &Alignment,
) -> (Vec<u8>, Vec<u8>) {
    let clipped_seq1 = &seq1[alignment.xstart..alignment.xend];
    let clipped_seq2 = &seq2[alignment.ystart..alignment.yend];

    let mut aligned_seq1 = Vec::new();
    let mut aligned_seq2 = Vec::new();

    let mut pos1 = 0;
    let mut pos2 = 0;

    for op in &alignment.operations {
        match op {
            Match | Subst => {
                aligned_seq1.push(clipped_seq1[pos1]);
                aligned_seq2.push(clipped_seq2[pos2]);
                pos1 += 1;
                pos2 += 1;
            }
            Ins => {
                // Insertion in the first sequence, gap in the second
                aligned_seq1.push(clipped_seq1[pos1]);
                aligned_seq2.push(b'-'); // Gap character
                pos1 += 1;
            }
            Del => {
                // Deletion in the first sequence (insertion in the second)
                aligned_seq1.push(b'-'); // Gap character
                aligned_seq2.push(clipped_seq2[pos2]);
                pos2 += 1;
            }
            // If Xclip or Yclip, do nothing (sequences were already clipped)
            _ => {}
        }
    }
    (aligned_seq1, aligned_seq2)
}

fn calculate_identity(alignment: &Alignment) -> f64 {
    let mut matching_positions = 0;
    let mut total_alignment_length = 0;

    for op in &alignment.operations {
        match op {
            // Matches contribute to both matching positions and total length
            Match => {
                matching_positions += 1;
                total_alignment_length += 1;
            }
            // Substitutions, insertions, and deletions only contribute to the
            // total alignment length
            Subst | Ins | Del => {
                total_alignment_length += 1;
            }
            // Xclip and Yclip are ignored
            _ => {}
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
    let alignment_type = AlignmentType::Local;

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
            let alignment = align_sequences(&pair[0].1, &pair[1].1, &alignment_type);
            let identity = calculate_identity(&alignment);
            let (aseq1, aseq2) = apply_alignment_to_sequences(&pair[0].1, &pair[1].1, &alignment);
            println!(
                "{}\t{}\t{:.6}\t{}\t{}",
                String::from_utf8_lossy(&pair[0].0),
                String::from_utf8_lossy(&pair[1].0),
                identity,
                String::from_utf8_lossy(&aseq1),
                String::from_utf8_lossy(&aseq2)
            );
        });
}
