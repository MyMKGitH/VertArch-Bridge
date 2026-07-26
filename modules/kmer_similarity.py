"""
kmer_similarity.py
-------------------
Alignment-free, reference-free sequence similarity via k-mer MinHash sketching.

Rationale (why this exists in this project):
The job description highlights "reference-free whole-genome alignments"
(the Vertebrate Genomes Project / Cactus-style approach) as central to
answering whether genome architecture tracks mutation or selection across
lineages. Full reference-free whole-genome alignment (e.g. Cactus, minimap2
whole-genome mode) is out of scope for a lightweight prototype, but the same
underlying idea — comparing sequences without a coordinate-anchored
alignment — has a well-known, tractable approximation: k-mer MinHash
sketching + Jaccard similarity (the same alignment-free strategy behind
tools like Mash, and behind MK's own FishAMR-Link plasmid-similarity module).

This module is intentionally the same technique applied one level up: from
comparing bacterial plasmids (FishAMR-Link) to comparing vertebrate genomic
sequences/homologous regions here. It is not a substitute for Cactus-style
alignment, but it is a legitimate, fast, dependency-light way to get a first
similarity signal between two sequences (e.g. orthologous regions from two
species) before reaching for a full alignment pipeline.
"""

import hashlib
from dataclasses import dataclass
from typing import Iterable, List, Set


def kmer_set(seq: str, k: int = 21) -> Set[str]:
    seq = seq.upper()
    if len(seq) < k:
        return set()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def minhash_signature(kmers: Iterable[str], num_hashes: int = 128, seed: int = 0) -> List[int]:
    """
    Bottom-sketch MinHash: for each of `num_hashes` independent hash functions
    (simulated here by salting SHA-1 with an index), keep the minimum hash
    value observed across all k-mers. Two sequences that share more k-mers
    will agree on more of these minima.
    """
    kmers = list(kmers)
    if not kmers:
        return [0] * num_hashes
    sig = []
    for i in range(num_hashes):
        salt = str(i + seed).encode()
        min_h = min(int(hashlib.sha1(salt + km.encode()).hexdigest(), 16) for km in kmers)
        sig.append(min_h)
    return sig


def estimated_jaccard(sig_a: List[int], sig_b: List[int]) -> float:
    """Fraction of matching minima across two MinHash signatures — an unbiased
    estimator of true k-mer Jaccard similarity, cheap to compute even for
    long sequences since it never touches the full k-mer sets directly."""
    if len(sig_a) != len(sig_b) or not sig_a:
        return 0.0
    matches = sum(1 for a, b in zip(sig_a, sig_b) if a == b)
    return matches / len(sig_a)


def exact_jaccard(seq_a: str, seq_b: str, k: int = 21) -> float:
    """Ground-truth Jaccard on the full k-mer sets — useful for validating the
    MinHash estimate on shorter sequences where computing it directly is cheap."""
    ka, kb = kmer_set(seq_a, k), kmer_set(seq_b, k)
    if not ka and not kb:
        return 0.0
    return len(ka & kb) / len(ka | kb)


@dataclass
class SimilarityResult:
    estimated_jaccard: float
    num_hashes: int
    k: int


def compare_sequences(seq_a: str, seq_b: str, k: int = 21, num_hashes: int = 128) -> SimilarityResult:
    sig_a = minhash_signature(kmer_set(seq_a, k), num_hashes=num_hashes)
    sig_b = minhash_signature(kmer_set(seq_b, k), num_hashes=num_hashes)
    return SimilarityResult(
        estimated_jaccard=round(estimated_jaccard(sig_a, sig_b), 4),
        num_hashes=num_hashes,
        k=k,
    )
