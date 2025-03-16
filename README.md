# all-pairwise-alignments

Computes all pairwise alignments between the sequences in a FASTA file.

## Usage

```sh
all-pairwise-alignments <fasta_file>
```

The output has a tabular format containing the folowing columns:

1. First protein name
2. Second protein name
3. CIGAR string representing the alignment
4. Alignment identity
5. Aligned sequence of the first protein
6. Aligned sequence of the second protein

## Notes

- The input FASTA file should contain protein sequences
- Global alignments are calculated using the [`block-aligner`](https://github.com/Daniel-Liu-c0deb0t/block-aligner) library
