"""
VertArch-Bridge
===============
A prototype bridging MK's existing bioinformatics toolchain (Python,
Biopython, Streamlit, Docker, alignment-free k-mer methods) toward the
Nord University evolutionary-genomics PhD project: whether vertebrate
genome architecture (genome size, GC content, chromosome number/size) is
better explained by mutational processes than by selection, using
chromosome-scale vertebrate assemblies and Nanopore-based germline
methylation.

This is a starting scaffold, not a finished research tool. Each tab maps to
one plank of that research question:

  1. Genome Architecture Explorer  -> genome size / GC% / chromosome number
                                       vs. a generation-time mutation-supply
                                       proxy, across species.
  2. CpG Island / Methylation Scanner -> sequence-level candidate regions for
                                       germline methylation, the natural
                                       precursor to real Nanopore methylation
                                       calls.
  3. Reference-free Similarity (MinHash) -> the same alignment-free k-mer
                                       technique already used in FishAMR-Link,
                                       applied here to compare genomic
                                       sequences without a full alignment —
                                       a lightweight analog of the
                                       reference-free whole-genome alignment
                                       approach central to the VGP.

Run locally:
    pip install -r requirements.txt
    streamlit run app.py
"""

import streamlit as st
import pandas as pd
import plotly.express as px

from modules import cpg_islands, kmer_similarity, genome_architecture

st.set_page_config(page_title="VertArch-Bridge", layout="wide")

st.title("VertArch-Bridge")
st.caption(
    "A prototype connecting existing bioinformatics tooling to the evolutionary-genomics "
    "question of whether vertebrate genome architecture tracks mutation supply or selection."
)

tab1, tab2, tab3 = st.tabs([
    "Genome Architecture Explorer",
    "CpG Island / Methylation Scanner",
    "Reference-free Similarity (MinHash)",
])

# ---------------------------------------------------------------------------
# TAB 1 — Genome Architecture Explorer
# ---------------------------------------------------------------------------
with tab1:
    st.subheader("Genome size, GC content, and chromosome number vs. a mutation-supply proxy")
    st.warning(
        "The bundled example dataset uses small, rounded, illustrative figures for a handful "
        "of well-known species — for demo purposes only. Replace with live-pulled NCBI/Ensembl/VGP "
        "assembly statistics before drawing any real conclusions (see genome_architecture.py).",
        icon="⚠️",
    )

    uploaded = st.file_uploader(
        "Optional: upload your own species CSV (same columns as the example dataset)",
        type=["csv"],
    )
    df = genome_architecture.load_dataset(uploaded) if uploaded else genome_architecture.load_example_dataset()
    st.dataframe(df, use_container_width=True)

    col1, col2 = st.columns(2)
    with col1:
        fig1 = px.scatter(
            df, x="generation_time_years", y="genome_size_mb",
            color="vertebrate_class", size="chromosome_number",
            hover_name="common_name", log_x=True,
            title="Genome size vs. generation time (mutation-supply proxy)",
        )
        st.plotly_chart(fig1, use_container_width=True)
    with col2:
        fig2 = px.scatter(
            df, x="gc_content_pct", y="chromosome_number",
            color="vertebrate_class", size="genome_size_mb",
            hover_name="common_name",
            title="GC content vs. chromosome number",
        )
        st.plotly_chart(fig2, use_container_width=True)

    st.markdown(
        "**How this maps to the project:** a mutation-driven hypothesis predicts genome "
        "architecture metrics correlate more strongly with proxies for mutation supply "
        "(here, generation time, as a very rough stand-in for germline replication rate) "
        "than with lineage-specific selective pressures. Swapping in real per-species "
        "germline mutation rate estimates and multiple sequentially aligned assemblies "
        "(rather than 7 illustrative rows) is the natural next step."
    )

# ---------------------------------------------------------------------------
# TAB 2 — CpG Island / Methylation Candidate Scanner
# ---------------------------------------------------------------------------
with tab2:
    st.subheader("Sliding-window CpG island / methylation-candidate detector")
    st.caption(
        "Flags CpG-dense candidate regions from raw sequence — the classic precursor step "
        "before layering real Oxford Nanopore methylation calls on top."
    )

    preset = st.radio(
        "Threshold preset",
        options=["strict", "classic"],
        index=0,
        help="strict = Takai & Jones (2002): GC≥55%, Obs/Exp≥0.65, len≥500bp — fewer false "
             "positives in large, repeat-rich genomes. classic = Gardiner-Garden & Frommer "
             "(1987): GC≥50%, Obs/Exp≥0.6, len≥200bp.",
        horizontal=True,
    )

    input_mode = st.radio("Input", ["Paste sequence", "Upload FASTA"], horizontal=True)
    seq_text = None
    if input_mode == "Paste sequence":
        seq_text = st.text_area("Paste a nucleotide sequence (FASTA header optional)", height=150)
        if seq_text and seq_text.strip().startswith(">"):
            seq_text = "".join(line for line in seq_text.splitlines() if not line.startswith(">"))
    else:
        fasta_file = st.file_uploader("Upload a FASTA file", type=["fa", "fasta", "txt"])
        if fasta_file:
            content = fasta_file.read().decode("utf-8", errors="ignore")
            seq_text = "".join(line for line in content.splitlines() if not line.startswith(">"))

    if seq_text and st.button("Scan for CpG islands"):
        regions = cpg_islands.scan_sequence(seq_text, preset=preset)
        if not regions:
            st.info("No candidate CpG islands found at this threshold.")
        else:
            result_df = pd.DataFrame([r.__dict__ for r in regions])
            st.dataframe(result_df, use_container_width=True)
            st.success(f"Found {len(regions)} candidate region(s).")

# ---------------------------------------------------------------------------
# TAB 3 — Reference-free Similarity (MinHash)
# ---------------------------------------------------------------------------
with tab3:
    st.subheader("Reference-free sequence similarity (k-mer MinHash / Jaccard)")
    st.caption(
        "The same alignment-free k-mer technique used for plasmid comparison in FishAMR-Link, "
        "applied here to compare genomic sequences directly — a lightweight analog of "
        "reference-free whole-genome alignment (e.g. the VGP's Cactus-based approach)."
    )

    col_a, col_b = st.columns(2)
    with col_a:
        seq_a = st.text_area("Sequence A", height=150, key="seq_a")
    with col_b:
        seq_b = st.text_area("Sequence B", height=150, key="seq_b")

    k = st.slider("k-mer size", min_value=11, max_value=31, value=21, step=2)
    num_hashes = st.slider("Number of MinHash functions", min_value=32, max_value=256, value=128, step=32)

    if seq_a and seq_b and st.button("Compare sequences"):
        result = kmer_similarity.compare_sequences(seq_a, seq_b, k=k, num_hashes=num_hashes)
        st.metric("Estimated Jaccard similarity", f"{result.estimated_jaccard:.3f}")
        st.caption(f"k={result.k}, {result.num_hashes} MinHash functions")
