"""
cpg_islands.py
--------------
Sliding-window CpG-island / methylation-candidate detector.

Rationale (why this exists in this project):
The Nord evolutionary-genomics position uses Oxford Nanopore sequencing to
measure germline DNA methylation directly. CpG islands are the classic,
well-established proxy for regions of regulatory / methylation interest —
they are typically hypomethylated in the germline and are the first thing
a genome-architecture study annotates before layering real Nanopore
methylation calls on top. This module is a lightweight, dependency-light
implementation of the same logic real tools (e.g. EMBOSS newcpgreport,
CpGcluster) use, so it is easy to read, audit, and extend.

Two published threshold sets are supported:
  - "classic"  : Gardiner-Garden & Frommer (1987)  — GC >= 50%, ObsExp >= 0.6
  - "strict"   : Takai & Jones (2002)               — GC >= 55%, ObsExp >= 0.65,
                 length >= 500 bp (reduces Alu/repeat false positives in
                 large, repeat-rich vertebrate genomes)

This is a detection heuristic, not a methylation caller — it flags CpG-dense
candidate regions from primary sequence alone. Pairing candidate regions
with real Nanopore methylation-calling output (e.g. from a tool such as
Nanopolish/f5c or dorado) is the natural next step once real Nanopore reads
are available, and is left as an extension point (see README).
"""

from dataclasses import dataclass
from typing import List, Tuple

PRESETS = {
    "classic": dict(min_gc=0.50, min_oe=0.60, min_len=200),
    "strict": dict(min_gc=0.55, min_oe=0.65, min_len=500),
}


@dataclass
class CpGRegion:
    start: int
    end: int
    length: int
    gc_fraction: float
    obs_exp_ratio: float


def _window_stats(window: str) -> Tuple[float, float]:
    c = window.count("C")
    g = window.count("G")
    cg_dinucs = sum(1 for i in range(len(window) - 1) if window[i] == "C" and window[i + 1] == "G")
    gc_fraction = (c + g) / len(window) if window else 0.0
    if c == 0 or g == 0:
        obs_exp = 0.0
    else:
        obs_exp = (cg_dinucs * len(window)) / (c * g)
    return gc_fraction, obs_exp


def scan_sequence(
    seq: str,
    window: int = 200,
    step: int = 50,
    preset: str = "strict",
) -> List[CpGRegion]:
    """
    Slide a window across `seq` and return merged candidate CpG-island regions.

    Args:
        seq: raw nucleotide sequence (any case; non-ACGT characters are ignored
             for composition counts but keep the window aligned to `seq` coordinates).
        window: window size in bp.
        step: step size in bp between successive windows.
        preset: "classic" or "strict" (see PRESETS above), or pass a dict with
                min_gc / min_oe / min_len to use custom thresholds.

    Returns:
        List of CpGRegion, merged where adjacent/overlapping windows both pass
        threshold, sorted by start position.
    """
    thresholds = PRESETS[preset] if isinstance(preset, str) else preset
    seq = seq.upper()
    n = len(seq)
    if n < window:
        return []

    raw_hits: List[Tuple[int, int, float, float]] = []
    for start in range(0, n - window + 1, step):
        w = seq[start:start + window]
        gc_frac, obs_exp = _window_stats(w)
        if gc_frac >= thresholds["min_gc"] and obs_exp >= thresholds["min_oe"]:
            raw_hits.append((start, start + window, gc_frac, obs_exp))

    merged: List[Tuple[int, int, float, float]] = []
    for s, e, gc, oe in raw_hits:
        if merged and s <= merged[-1][1]:
            ps, pe, pgc, poe = merged[-1]
            merged[-1] = (ps, max(pe, e), max(pgc, gc), max(poe, oe))
        else:
            merged.append((s, e, gc, oe))

    regions = [
        CpGRegion(start=s, end=e, length=e - s, gc_fraction=round(gc, 3), obs_exp_ratio=round(oe, 3))
        for s, e, gc, oe in merged
        if (e - s) >= thresholds["min_len"]
    ]
    return regions


def scan_fasta(path: str, window: int = 200, step: int = 50, preset: str = "strict") -> dict:
    """
    Run scan_sequence over every record in a FASTA file.
    Uses Biopython if available; falls back to a minimal built-in FASTA
    reader so the module has no hard dependency for this function alone.
    Returns {record_id: [CpGRegion, ...]}.
    """
    records = {}
    try:
        from Bio import SeqIO
        for rec in SeqIO.parse(path, "fasta"):
            records[rec.id] = str(rec.seq)
    except ImportError:
        records = _minimal_fasta_reader(path)

    return {rid: scan_sequence(seq, window=window, step=step, preset=preset) for rid, seq in records.items()}


def _minimal_fasta_reader(path: str) -> dict:
    records, current_id, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(buf)
                current_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if current_id is not None:
            records[current_id] = "".join(buf)
    return records
