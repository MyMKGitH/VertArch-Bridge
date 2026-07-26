# VertArch-Bridge v1.0.0
> **Comparative Vertebrate Genome Architecture & Alignment-Free Sequence Analysis Suite**
> 

## 📌 Executive Summary
**VertArch-Bridge** is an interactive computational framework engineered to explore fundamental questions in evolutionary genomics: **Is vertebrate genome architecture (genome size, GC content, chromosome topology) primarily shaped by mutational processes or by natural selection?**
The application bridges macro-scale genomic metadata analysis with micro-scale sequence profiling. It provides a three-tiered analytical pipeline combining genome-wide evolutionary trait charting, sliding-window CpG island scanning for germline methylation analysis, and reference-free k-mer MinHash sequence similarity estimation.
## 🛠️ Core Analytical Modules
```
VertArch-Bridge Architecture
 ├── 1. Genome Architecture Explorer  ---> [Macro-Scale] Assembly Stats vs. Mutation Supply
 ├── 2. CpG Island & Methylation Scanner --> [Micro-Scale] Windowed GC% & Obs/Exp CpG Profiling
 └── 3. MinHash Sequence Similarity   ---> [Alignment-Free] Bottom-k SHA-1 k-mer Sketching

```
### 1. Genome Architecture Explorer (genome_architecture.py)
 * **Objective:** Analyzes whether macro-genomic architecture (genome size in Mb, GC%, chromosome count) correlates with generation-time proxies for germline mutation supply across vertebrate classes (*Mammalia*, *Aves*, *Actinopterygii*).
 * **Live Ingestion:** Integrated with the **NCBI Datasets API v2** to retrieve live assembly metrics across user-specified taxa.
### 2. CpG Island & Methylation Scanner (cpg_islands.py)
 * **Objective:** Identifies hypomethylated germline candidate regions—the foundational precursor for downstream Nanopore methylation profiling.
 * **Supported Presets:**
   * **Classic Preset** (*Gardiner-Garden & Frommer, 1987*): Window \ge 200\text{ bp}, \text{GC\%} \ge 50\%, \text{Obs/Exp CpG} \ge 0.60.
   * **Strict Preset** (*Takai & Jones, 2002*): Min length \ge 500\text{ bp}, \text{GC\%} \ge 55\%, \text{Obs/Exp CpG} \ge 0.65 (mitigates repeat/Alu false positives in large vertebrate genomes).
 * **Mathematical Formula:**
   
### 3. Reference-Free MinHash Similarity Engine (kmer_similarity.py)
 * **Objective:** Computes fast, coordinate-free sequence containment and Jaccard similarity between genomic regions without requiring computationally expensive whole-genome alignments (e.g., Cactus/minimap2).
 * **Algorithm:** Bottom-k MinHash sketching using SHA-1 salted hash functions over sliding k-mers (11 \le k \le 31).
 * **Jaccard Estimator:**
   
## 📁 Tabs

1. **Genome Architecture Explorer** — plots genome size / GC% / chromosome
   number against a generation-time proxy for mutation supply, across a
   small illustrative set of vertebrates. Swap in the bundled example CSV
   for live-pulled NCBI Datasets / Ensembl / VGP assembly statistics
   (`modules/genome_architecture.py::fetch_ncbi_assembly_stats` is an
   untested starting point — verify the endpoint against the current NCBI
   Datasets API docs before relying on it).

2. **CpG Island / Methylation Candidate Scanner** — a sliding-window CpG-
   island detector (Gardiner-Garden & Frommer 1987 / Takai & Jones 2002
   thresholds) over pasted or uploaded FASTA sequence. This is the classic
   sequence-level precursor to real Oxford Nanopore methylation calls —
   pairing candidate regions here with actual Nanopore methylation output
   (e.g. via Nanopolish/f5c/dorado) is the natural next step once real
   long-read data is available.

3. **Reference-free Similarity (MinHash)** — the same alignment-free k-mer
   Jaccard technique already used for plasmid comparison in FishAMR-Link,
   applied one level up to compare genomic sequences directly, as a
   lightweight analog of the reference-free whole-genome alignment approach
   (Cactus / VGP) central to the position.

## 🚀 Quick Start & Installation
### Local Setup
```bash
# Clone repository
git clone https://github.com/MyMKGitH/VertArch-Bridge.git
cd VertArch-Bridge

# Install dependencies
pip install -r requirements.txt

# Launch Streamlit interface
streamlit run app.py

```
### Docker Deployment
```bash
# Build Docker image
docker build -t vertarch-bridge:v1.0.0 .

# Run containerized application
docker run -d -p 8501:8501 vertarch-bridge:v1.0.0

```
## 📁 Repository Structure
```
VertArch-Bridge/
│
├── app.py                      # Main Streamlit multi-tab application UI
├── cpg_islands.py              # Sliding-window CpG island candidate detector
├── genome_architecture.py      # NCBI API fetcher & macro-genomic trait loader
├── kmer_similarity.py          # Alignment-free MinHash k-mer similarity estimator
├── example_species.csv         # Illustrative vertebrate genomic dataset
├── Dockerfile                  # Production containerization specification
├── requirements.txt            # Python dependencies
└── README.md                   # Technical documentation

```
## 🔮 Future Development Roadmap
 * [ ] **VGP Assembly Pipeline Integration:** Direct ingestion of Vertebrate Genomes Project (VGP) chromosome-level assemblies.
 * [ ] **Nanopore Methylation Validation:** Pairing predicted CpG candidate regions with direct Oxford Nanopore modified base calls (e.g., dorado / modkit).
 * [ ] **Phylogenetic Comparative Controls:** Incorporating Phylogenetic Generalized Least Squares (PGLS) to adjust for evolutionary non-independence.
## 📜 License
Distributed under the **MIT License**.
