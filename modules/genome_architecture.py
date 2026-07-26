"""
genome_architecture.py
-----------------------
Loads vertebrate genome-architecture metrics (genome size, GC content,
chromosome number, a generation-time proxy for mutation supply) and gives a
minimal scaffold for testing the mutation-vs-selection framing the job
description is built around: is genome architecture better explained by
the processes that generate mutations (approximated here, crudely, by
generation time as a proxy for germline replication/mutation-accumulation
rate) than by the selection that later acts on it?

IMPORTANT — data honesty note:
`data/example_species.csv` contains small, deliberately ROUNDED, illustrative
figures for a handful of well-known vertebrates, for demo/UI purposes only.
They are NOT sourced from a verified database call in this environment and
must not be cited or relied on for real analysis. Before using this for
actual research, replace them with assembly statistics pulled live from
NCBI Datasets, Ensembl, or the VGP assembly reports themselves — see
`fetch_ncbi_assembly_stats()` below for a starting point (untested against
a live endpoint in this sandbox — verify the request/response shape before
relying on it).
"""

import os
import pandas as pd

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "example_species.csv")


def load_example_dataset() -> pd.DataFrame:
    """Loads the small illustrative example dataset shipped with this prototype."""
    return pd.read_csv(DATA_PATH)


def load_dataset(path: str) -> pd.DataFrame:
    """Loads a user-supplied CSV with the same column schema as example_species.csv:
    species, common_name, vertebrate_class, genome_size_mb, gc_content_pct,
    chromosome_number, generation_time_years
    """
    return pd.read_csv(path)


def fetch_ncbi_assembly_stats(taxon: str) -> dict:
    """
    Best-effort live fetch of genome assembly statistics from the NCBI
    Datasets API v2 for a given species/taxon name.

    NOTE: this endpoint shape has not been exercised against a live server in
    the environment this prototype was built in (no network access at build
    time). Treat this as a starting scaffold — confirm the exact request URL,
    auth requirements, and response schema against the current NCBI Datasets
    API docs (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/) before
    relying on its output, and adjust the parsing below to match.
    """
    import requests

    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxon}/dataset_report"
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    payload = resp.json()

    reports = payload.get("reports", [])
    if not reports:
        raise ValueError(f"No assembly reports returned for taxon '{taxon}'")

    top = reports[0]
    stats = top.get("assembly_stats", {})
    return {
        "taxon": taxon,
        "organism_name": top.get("organism", {}).get("organism_name"),
        "genome_size_mb": stats.get("total_sequence_length"),
        "gc_content_pct": stats.get("gc_percent"),
        "chromosome_number": stats.get("total_number_of_chromosomes"),
        "assembly_accession": top.get("accession"),
        "source": "NCBI Datasets v2 (live, verify against current docs)",
    }
