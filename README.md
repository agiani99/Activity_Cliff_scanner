# Activity Cliff Scanner

A three-script pipeline for detecting, enriching, and visualising **activity cliffs** in ChEMBL — structurally similar compound pairs with large potency differences — using ECFP4 fingerprints, RDKit, and the ChEMBL 37 SQLite database.
Final output is a dashboard ![dashboard]('https://github.com/agiani99/Activity_Cliff_scanner/blob/main/Screenshot 2026-08-30 173354.png') can be helpful

---

## Table of contents

1. [What is an activity cliff?](#what-is-an-activity-cliff)
2. [Scripts overview](#scripts-overview)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Script 1 — `activity_cliff_scanner.py`](#script-1--activity_cliff_scannerpy)
6. [Script 2 — `cliff_enrich.py`](#script-2--cliff_enrichpy)
7. [Script 3 — `cliff_dashboard.py`](#script-3--cliff_dashboardpy)
8. [Recommended workflow](#recommended-workflow)
9. [Output column reference](#output-column-reference)
10. [Scientific notes](#scientific-notes)

---

## What is an activity cliff?

An activity cliff is a pair of molecules that are **structurally very similar** but show a **large difference in biological potency**. They are valuable for medicinal chemistry because they identify structural features that are critical for activity — a single atom or bond change that switches a nanomolar compound into a micromolar one.

### Definition used in this pipeline

| Criterion | Threshold | Note |
|---|---|---|
| Tanimoto similarity (ECFP4, r=2, 2048 bits) | ≥ 0.95 | Single-atom changes typically give 0.95–0.99 |
| \|ΔpXC50\| | ≥ 2.0 log units | 100-fold potency difference |
| No structural–activity intermediate | optional | No third compound bridges the gap |
| Assay confidence score | ≥ 8 | Direct assay evidence on single protein |

### Cliff types

| Type | Active | Inactive | Comparability |
|---|---|---|---|
| `XC50` | pXC50 from dose-response assay | pXC50 from **same assay** | Highest — identical conditions |
| `pct_fallback` | pXC50 ≥ 5.0 (≤ 10 µM) | % inhibition < 50 % at ~10 µM | Lower — cross-assay |

---

## Scripts overview

```
activity_cliff_scanner.py   ──►  all_cliffs.csv
                                      │
                   cliff_enrich.py ◄──┘──► all_cliffs_enriched.csv
                                                      │
                          cliff_dashboard.py ◄────────┘──► dashboard.html
```

---

## Requirements

```
Python >= 3.9
rdkit >= 2023.03
pandas
numpy
requests
tqdm
```

**No `chembl-webresource-client` required.** The scanner uses direct REST API calls and/or a local SQLite database.

---

## Installation

```bash
pip install rdkit pandas numpy requests tqdm
```

**Strongly recommended**: download the ChEMBL 37 SQLite database (~5 GB uncompressed) for fastest, network-free operation:

```
https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/
```

**Optional but highly recommended**: download SIFTS flat files for offline protein family classification and PDB structure annotation (see [cliff_enrich.py](#script-2--cliff_enrichpy)):

```bash
# Download to a local .\MAPS\ folder
cd MAPS
curl -O https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_pdb.csv.gz
curl -O https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_enzyme.csv.gz
curl -O https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_pfam.csv.gz
curl -O https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz
gzip -d *.gz
```

---

## Script 1 — `activity_cliff_scanner.py`

Scans ChEMBL assays and writes cliff pairs to CSV with incremental checkpointing.

### Data sources (in order of preference)

| Source | Flag | Speed | Network |
|---|---|---|---|
| ChEMBL 37 SQLite | `--sqlite PATH` | Fastest | None |
| Pre-exported CSV | `--input-csv FILE` | Fast | None |
| ChEMBL REST API | *(default)* | Slow | Required |
| Built-in demo data | `--demo` | Instant | None |

### Usage

```bash
# Full ChEMBL 37, EGFR only — recommended first run
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --target CHEMBL279 \
    --output EGFR_cliffs.csv

# All targets — full proteome scan (hours)
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --output all_cliffs.csv

# Strictest filter: same-assay XC50 pairs only
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --xc50-only \
    --output all_cliffs_strict.csv

# Middle ground: same assay-type (B↔B or F↔F)
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --same-assay-type \
    --output all_cliffs_typed.csv

# Resume an interrupted run
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --resume \
    --output all_cliffs.csv

# Demo (validates pipeline, no data needed)
python activity_cliff_scanner.py --demo
```

### CLI flags

| Flag | Default | Description |
|---|---|---|
| `--sqlite PATH` | — | Path to ChEMBL `.db` file or directory |
| `--input-csv FILE` | — | Pre-exported activity CSV |
| `--demo` | — | Run on built-in example data |
| `--target CHEMBL_ID` | all targets | Restrict to one ChEMBL target (e.g. `CHEMBL279`) |
| `--output FILE` | `activity_cliffs.csv` | Output CSV path |
| `--tanimoto FLOAT` | `0.95` | ECFP4 Tanimoto cutoff |
| `--delta-pxic50 FLOAT` | `2.0` | Minimum \|ΔpXC50\| in log units |
| `--no-intermediate-check` | off | Skip Criterion 3 (faster, more pairs) |
| `--xc50-only` | off | Only XC50–XC50 same-assay pairs; disables pct_fallback entirely |
| `--same-assay-type` | off | pct_fallback only between matching assay types (B↔B, F↔F) |
| `--max-records INT` | 0 (unlimited) | Cap rows fetched from SQLite or API |
| `--resume` | off | Restart from checkpoint; appends to existing CSV |

### Assay-type strictness hierarchy

| Flag | Filter applied | Expected pairs |
|---|---|---|
| *(none)* | None — any cross-assay pair | Most (268 K+ in full proteome run) |
| `--same-assay-type` | B↔B or F↔F only | Subset — removes cross-biology pairs |
| `--xc50-only` | Same assay, XC50 only | Fewest — highest confidence |

### Measurement types included

IC50, EC50, DC50, AC50, CC50, GI50, Ki, Kd, Kb, ED50, Potency

*MIC excluded* — not convertible to nM without molecular weight.

### Identity filters (duplicate prevention)

Every pair is rejected if any of the following holds:

1. Same `molecule_chembl_id`
2. Tanimoto = 1.0 (identical fingerprint)
3. Same raw InChIKey (same molecular graph)
4. Same **parent** InChIKey (salt / free-base / zwitterion forms of the same compound)

### Checkpointing

Two files are created alongside the output:

| File | Contents |
|---|---|
| `all_cliffs.ckpt` | Processed assay ChEMBL IDs (one per line) |
| `all_cliffs_raw.csv` | XC50–XC50 pairs before enrichment |

Use `--resume` to skip already-processed assays on restart. Nothing is lost if the run is interrupted.

### Assay confidence scores

| Score | Meaning |
|---|---|
| 9 | Direct assay, single protein, exact compound tested |
| 8 | Direct assay, single protein |
| **≤ 7** | **Excluded** — homologous protein, multi-protein complex, or inferred |

---

## Script 2 — `cliff_enrich.py`

Post-hoc enrichment of the scanner output. Adds seven annotation layers without re-running the scanner.

### What it adds

| Column(s) | Source | Network |
|---|---|---|
| `tanimoto_ecfp4` | RDKit ECFP4 from SMILES | None |
| `active_assay_chembl_id` | Renamed from `assay_chembl_id` | None |
| `inactive_assay_chembl_id` | Same as active (XC50–XC50); SQLite lookup for pct_fallback | SQLite only |
| `protein_class_l1/l2/l3`, `protein_family` | **SIFTS EC + Pfam** (primary, offline) → SQLite → ChEMBL REST | SIFTS preferred |
| `n_target_pdb_structures` | SIFTS UniProt→PDB count | None |
| `target_pdb_sample` | First 5 PDB IDs for target | None |
| `active_compound_pdb` | RCSB combined chemical + UniProt search | RCSB REST |
| `active_pdb_resolution_A` | Resolution of co-crystal structure | RCSB REST |
| `active_pdb_method` | X-RAY / ELECTRON MICROSCOPY / NMR | RCSB REST |

### Usage

```bash
# Full enrichment — SIFTS classification + SQLite + PDB annotation (recommended)
python cliff_enrich.py \
    --input all_cliffs.csv \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --maps-dir ".\MAPS" \
    --output all_cliffs_enriched.csv

# SIFTS offline only — no network, no SQLite (Tanimoto + protein class + PDB counts)
python cliff_enrich.py \
    --input all_cliffs.csv \
    --maps-dir ".\MAPS" \
    --no-compound-pdb \
    --output all_cliffs_enriched.csv

# Tanimoto + assay IDs only (fastest, no classification)
python cliff_enrich.py \
    --input all_cliffs.csv \
    --skip-protein-class \
    --output all_cliffs_enriched.csv
```

### CLI flags

| Flag | Default | Description |
|---|---|---|
| `--input FILE` | required | Input cliff CSV from scanner |
| `--output FILE` | `<input>_enriched.csv` | Output CSV |
| `--sqlite PATH` | — | ChEMBL SQLite for inactive assay IDs and protein class fallback |
| `--maps-dir DIR` | — | Directory with SIFTS flat files (see below) |
| `--sifts FILE` | — | Single SIFTS file for PDB structure counts only (superseded by `--maps-dir`) |
| `--skip-protein-class` | off | Skip all protein class annotation |
| `--skip-uniprot-family` | off | Skip UniProt SIMILARITY REST fallback |
| `--no-compound-pdb` | off | Skip per-compound RCSB search (faster; still annotates target PDB counts) |

### Protein class annotation — priority order

```
1. SIFTS offline (--maps-dir)          pdb_chain_enzyme.csv  → EC number  → Kinase / Protease / PDE
                                        pdb_chain_pfam.csv    → Pfam domain → GPCR / NHR / Ion channel
   Coverage: >90% for drug targets. Fully offline. Seconds.

2. ChEMBL SQLite  protein_class table  Two-step: target → protein_class_id → L1/L2/L3
   Coverage: depends on SQLite build; may be absent.

3. ChEMBL REST    /protein_class/{id}  Two API calls per unique class ID (cached).
   Coverage: good but slow (~1 s per target).

4. UniProt REST   SIMILARITY comment   "Belongs to the protein kinase superfamily."
   Only for targets still unresolved after steps 1–3.
```

### SIFTS files used by `--maps-dir`

Auto-discovered inside the directory by filename. Plain `.csv` and `.gz` both accepted.

| File | Used for |
|---|---|
| `uniprot_pdb.csv` or `pdb_chain_uniprot.csv` | UniProt → PDB ID list (target PDB counts, RCSB search filtering) |
| `pdb_chain_enzyme.csv` | EC numbers → Kinase / Protease / PDE / Phosphatase / ... |
| `pdb_chain_pfam.csv` | Pfam domains → GPCR / Nuclear hormone receptor / Ion channel / ... |

Download all from: `https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/`

### EC number → protein class mapping (selected)

| EC prefix | L2 label |
|---|---|
| 2.7.10.1 | Receptor Tyrosine Kinase |
| 2.7.11.1 | Serine/Threonine Kinase |
| 2.7.12 | Dual-specificity Kinase |
| 3.4.21 | Serine Protease |
| 3.4.24 | Metalloprotease |
| 3.1.4.17 | Cyclic Nucleotide PDE |
| 3.1.3.48 | Tyrosine Phosphatase |
| 3.5.1.98 | HDAC |
| 1.14.13/14 | Cytochrome P450 |
| 3.6.5 | GTPase |

### Pfam → protein class mapping (selected)

| Pfam ID | L2 label |
|---|---|
| PF00069, PF07714 | Kinase |
| PF00001/2/3 | GPCR (family A/B/C) |
| PF00104, PF00105 | Nuclear Hormone Receptor |
| PF00520 | Ion Channel |
| PF00089, PF00026 | Protease |
| PF00233 | Phosphodiesterase |
| PF02132 | HDAC |
| PF00439, PF00628 | Bromodomain |

### Performance on 268 K rows

| Step | Source | Time |
|---|---|---|
| Tanimoto (vectorised unique-SMILES cache) | RDKit | ~60–120 s |
| Protein class via SIFTS EC + Pfam | Local files | ~30–60 s |
| Target PDB counts from SIFTS | Local files | ~5 s |
| Inactive assay IDs | SQLite | ~30 s |
| Active compound in PDB (RCSB search) | Network | ~1 s per unique pair |

---

## Script 3 — `cliff_dashboard.py`

Generates a **self-contained HTML dashboard** from the cliff CSV. No server required — open in any browser.

### Features

- Side-by-side 2D depictions with **MCS scaffold** (pale blue) and **diff atoms** highlighted (green = active-unique, red = inactive-unique)
- 2D layouts **aligned on MCS scaffold** — both images drawn in the same orientation
- **RMSD** from lowest-energy MMFF conformer alignment on scaffold atoms (ETKDGv3)
- MCS found via two-pass strategy: strict ring matching first, then ring-size-relaxed fallback (catches cyclohexyl→cyclopentyl cliffs); timed-out partial results accepted
- Sortable / filterable DataTable (Bootstrap 5 + DataTables)
- Click any row → full-size modal with detail view
- Self-contained HTML (no local server needed)

### Usage

```bash
# Top 200 pairs with 3D RMSD (~2 min)
python cliff_dashboard.py \
    --input all_cliffs_enriched.csv \
    --output dashboard.html

# Top 500, skip 3D (faster)
python cliff_dashboard.py \
    --input all_cliffs_enriched.csv \
    --top-n 500 --skip-3d \
    --output dashboard.html

# Single target, all pairs
python cliff_dashboard.py \
    --input all_cliffs_enriched.csv \
    --target CHEMBL279 --top-n 0 \
    --output EGFR_dashboard.html

# XC50-only pairs with 3D RMSD
python cliff_dashboard.py \
    --input all_cliffs_enriched.csv \
    --cliff-type XC50 --top-n 200 \
    --output dashboard_xc50.html
```

### CLI flags

| Flag | Default | Description |
|---|---|---|
| `--input FILE` | required | Cliff CSV (scanner or enricher output) |
| `--output FILE` | `<input>_dashboard.html` | Output HTML |
| `--top-n INT` | 200 | Pairs to render (0 = all; >500 makes large files) |
| `--min-delta FLOAT` | 2.0 | Minimum \|ΔpXC50\|; NaN-delta pct_fallback pairs are always kept |
| `--target CHEMBL_ID` | all | Filter to one target |
| `--cliff-type` | `all` | `XC50`, `pct_fallback`, or `all` |
| `--skip-3d` | off | Skip ETKDGv3 + MMFF conformer generation and RMSD |

### RMSD interpretation

| RMSD (Å) | Interpretation |
|---|---|
| < 0.5 | Scaffold rigid; cliff is likely **electronic** (changed atom alters charge, H-bond, polarisability) |
| 0.5–1.5 | Mixed — electronic + minor conformational perturbation |
| > 1.5 | Structural change **perturbs the preferred conformation** — possible binding-mode change |

---

## Recommended workflow

```bash
# Step 1: scan — same assay-type filter (B↔B, F↔F)
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --same-assay-type \
    --output all_cliffs.csv

# Step 2: enrich — SIFTS offline protein class + PDB annotation
python cliff_enrich.py \
    --input all_cliffs.csv \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --maps-dir ".\MAPS" \
    --no-compound-pdb \
    --output all_cliffs_enriched.csv

# Step 3: dashboard — XC50 pairs, top 200 with 3D RMSD
python cliff_dashboard.py \
    --input all_cliffs_enriched.csv \
    --cliff-type XC50 \
    --top-n 200 \
    --output dashboard_xc50.html
```

---

## Output column reference

### Scanner output (`all_cliffs.csv`)

| Column | Type | Notes |
|---|---|---|
| `active_chembl_id` | str | ChEMBL compound ID of the more potent molecule |
| `active_smiles` | str | Canonical SMILES (registered form, salt intact) |
| `active_pXC50` | float | −log₁₀(IC50/EC50/Ki/Kd in M); ChEMBL `pchembl_value` preferred |
| `active_value_nM` | float | Raw potency in nM |
| `active_std_type` | str | IC50, EC50, Ki, Kd, etc. |
| `inactive_chembl_id` | str | ChEMBL compound ID of the less potent molecule |
| `inactive_smiles` | str | Canonical SMILES |
| `inactive_pXC50` | float | Populated for XC50 pairs; NaN for pct_fallback |
| `inactive_value_nM` | float | Raw potency in nM; NaN for pct_fallback |
| `inactive_std_type` | str | Measurement type |
| `inactive_pct_activity` | float | % inhibition at ~10 µM; pct_fallback only |
| `delta_pXC50` | float | \|active_pXC50 − inactive_pXC50\|; NaN for pct_fallback |
| `assay_chembl_id` | str | ChEMBL assay ID (renamed to `active_assay_chembl_id` by enricher) |
| `target_chembl_id` | str | ChEMBL target ID |
| `target_name` | str | Target preferred name |
| `target_organism` | str | e.g. Homo sapiens |
| `uniprot_id` | str | UniProt accession(s); pipe-separated for multi-component targets |
| `assay_confidence_score` | int | 8 or 9 |
| `active_n_tautomers` | int | Number of enumerated tautomers (RDKit, cap 64) |
| `active_tautomer_important` | bool | True if ≥ 2 tautomers |
| `inactive_n_tautomers` | int | — |
| `inactive_tautomer_important` | bool | — |
| `active_charge_ph7` | str | Ionisable groups at pH 7 with pKa and expected charge |
| `inactive_charge_ph7` | str | — |
| `cliff_type` | str | `XC50` or `pct_fallback` |

### Additional columns after `cliff_enrich.py`

| Column | Source | Notes |
|---|---|---|
| `tanimoto_ecfp4` | RDKit | ECFP4 Tanimoto similarity (4 d.p.) |
| `active_assay_chembl_id` | Renamed | Assay of the active measurement |
| `inactive_assay_chembl_id` | SQLite / same | Assay of the inactive; differs for pct_fallback |
| `protein_class_l1` | SIFTS / SQLite / REST | Enzyme, Membrane receptor, Nuclear receptor, Ion channel, Transporter |
| `protein_class_l2` | SIFTS / SQLite / REST | **Kinase, Protease, GPCR, NHR, Phosphodiesterase, ...** |
| `protein_class_l3` | SIFTS / SQLite / REST | Serine Protease, Family A GPCR, Cyclic Nucleotide PDE, ... |
| `protein_family` | SIFTS / REST | Protein Kinase, Metalloprotease, NHR ligand-binding domain, ... |
| `n_target_pdb_structures` | SIFTS | Count of PDB entries for target UniProt accession |
| `target_pdb_sample` | SIFTS | Comma-separated list of first 5 PDB IDs |
| `active_compound_pdb` | RCSB REST | PDB ID where exact active compound is co-crystallised with target |
| `active_pdb_resolution_A` | RCSB REST | Resolution (Å) of co-crystal structure |
| `active_pdb_method` | RCSB REST | X-RAY DIFFRACTION / ELECTRON MICROSCOPY / SOLUTION NMR |

---

## Scientific notes

### Why Tanimoto ≥ 0.95?

At this threshold, pairs typically differ by a **single atom or small functional group** (F→Cl, H→CH₃, CH₂ insertion, ring-size change). This makes structural differences unambiguous and directly interpretable. Lowering to 0.85 admits scaffold hops where the activity difference may not be attributable to one specific structural change.

### Why confidence score ≥ 8?

Scores 8 and 9 indicate **direct binding/functional measurement on a single, defined protein target** (not inferred from cellular readouts or homology). Lower scores introduce uncertainty about which molecular target is responsible for the observed activity.

### Assay-type filter rationale

ChEMBL `assay_type` codes relevant to potency cliffs:

| Code | Biology | Typical measurements |
|---|---|---|
| B | Biochemical binding | Ki, Kd, SPR, ITC, TR-FRET |
| F | Functional | IC50, EC50, cell-based, enzyme activity |

Comparing a Ki (B) with a cellular IC50 (F) conflates target-binding affinity with cell penetration, efflux, and off-target effects. Using `--same-assay-type` keeps B↔B and F↔F comparisons only.

### pct_fallback pairs

These pairs have a quantified pXC50 for the active compound but only a single-point % inhibition measurement for the inactive. They are useful for identifying active compounds with clear structural neighbours that failed to show any activity at 10 µM. They are scientifically weaker than XC50–XC50 pairs because the "inactivity" may reflect insufficient dose rather than true structural failure. Use `--xc50-only` for high-confidence SAR analysis.

### Salt and duplicate handling

Before fingerprint computation, molecules are deduplicated per assay using three passes:

1. Same `molecule_chembl_id` + assay → keep most potent
2. Same raw InChIKey + assay → keep most potent
3. Same **parent** InChIKey (largest fragment + neutralised) + assay → keep most potent

This collapses free-base / hydrochloride / sodium salt forms of the same compound to a single entry, preventing spurious "cliffs" from salt differences.

### SIFTS-based protein classification

EC numbers from `pdb_chain_enzyme.csv` are the primary protein class source because they are assigned experimentally and encode enzyme function precisely. Pfam domain IDs from `pdb_chain_pfam.csv` supplement for non-enzymes (GPCRs, nuclear receptors, ion channels) which have no EC number by definition. The combination gives >90% coverage for drug targets in ChEMBL.

The `active_compound_pdb` RCSB search uses a combined filter: chemical (InChIKey exact match) AND protein entity (UniProt accession). This prevents false positives where the same ligand is found in a different protein target.
