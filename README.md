# Activity Cliff Scanner

A three-script pipeline for detecting, enriching, and visualising **activity cliffs** in ChEMBL — structurally similar compound pairs with large potency differences — using ECFP4 fingerprints, RDKit, and the ChEMBL 37 SQLite database.

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
| No structural–activity intermediate | required | No third compound bridges the gap (optional) |
| Assay confidence score | ≥ 8 | Direct assay evidence on single protein |

### Cliff types

| Type | Active | Inactive | Comparability |
|---|---|---|---|
| `XC50` | pXC50 from dose-response assay | pXC50 from same assay | Highest — identical conditions |
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

**Recommended**: download the ChEMBL 37 SQLite database (~5 GB uncompressed) for fastest, network-free operation:

```
https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/
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
| `--xc50-only` | off | Only XC50–XC50 same-assay pairs; disables pct_fallback |
| `--same-assay-type` | off | pct_fallback only between matching assay types (B↔B, F↔F) |
| `--max-records INT` | 0 (unlimited) | Cap rows fetched from SQLite or API |
| `--resume` | off | Restart from checkpoint; appends to existing CSV |

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

ChEMBL confidence scores for target assignment:

| Score | Meaning |
|---|---|
| 9 | Direct assay, single protein, exact compound tested |
| 8 | Direct assay, single protein |
| **≤ 7** | **Excluded** — homologous protein, multi-protein complex, or inferred |

---

## Script 2 — `cliff_enrich.py`

Post-hoc enrichment of the scanner output. Adds four annotation layers without re-running the scanner.

### What it adds

| Column | Source | Network needed |
|---|---|---|
| `tanimoto_ecfp4` | RDKit ECFP4 from SMILES | No |
| `active_assay_chembl_id` | Renamed from `assay_chembl_id` | No |
| `inactive_assay_chembl_id` | Same as active (XC50–XC50); SQLite lookup for pct_fallback | SQLite only |
| `protein_class_l1/l2/l3` | ChEMBL `protein_class` table | SQLite preferred |
| `protein_family` | ChEMBL hierarchy + UniProt SIMILARITY comment | UniProt REST fallback |

### Usage

```bash
# Full enrichment with SQLite (recommended)
python cliff_enrich.py \
    --input all_cliffs.csv \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --output all_cliffs_enriched.csv

# Tanimoto only, no network
python cliff_enrich.py \
    --input all_cliffs.csv \
    --skip-protein-class \
    --output all_cliffs_enriched.csv
```

### CLI flags

| Flag | Default | Description |
|---|---|---|
| `--input FILE` | required | Input cliff CSV |
| `--output FILE` | `<input>_enriched.csv` | Output CSV |
| `--sqlite PATH` | — | ChEMBL SQLite for protein class + inactive assay IDs |
| `--skip-protein-class` | off | Skip protein class lookup (Steps 3–4) |
| `--skip-uniprot-family` | off | Skip UniProt SIMILARITY comment lookup |

### Protein class fallback strategy

1. `protein_class` table (full L1→L3 hierarchy) — present in most ChEMBL 37 SQLite builds  
2. `target_class` table (denormalized view) — some builds  
3. `target_dictionary.target_type` (broad: SINGLE PROTEIN, GPCR…) + UniProt REST for `protein_family`

### Performance on 268 K rows

| Step | Time |
|---|---|
| Tanimoto (vectorised, unique-SMILES cache) | ~60–120 s |
| Protein class from SQLite | < 5 s (one batch query) |
| UniProt family (per unique accession, cached) | ~1 s/ID |
| Inactive assay IDs from SQLite | < 30 s |

---

## Script 3 — `cliff_dashboard.py`

Generates a **self-contained HTML dashboard** from the cliff CSV. No server required — open in any browser.

### Features

- Side-by-side 2D depictions with **MCS scaffold** (pale blue) and **diff atoms** highlighted (green = active-unique, red = inactive-unique)
- 2D layouts are **aligned on the MCS scaffold** — both images are drawn in the same orientation for direct visual comparison
- **RMSD** from lowest-energy MMFF conformer alignment on scaffold atoms
- Sortable / filterable DataTable (Bootstrap 5 + DataTables)
- Click any row → full-size modal with detail view
- Self-contained HTML (no local server needed)

### Usage

```bash
# Top 200 pairs, with 3D RMSD (~2 min)
python cliff_dashboard.py \
    --input all_cliffs.csv \
    --output dashboard.html

# Top 500, skip 3D (faster)
python cliff_dashboard.py \
    --input all_cliffs.csv \
    --top-n 500 --skip-3d \
    --output dashboard.html

# Single target, all pairs, no 3D
python cliff_dashboard.py \
    --input all_cliffs.csv \
    --target CHEMBL279 --top-n 0 --skip-3d \
    --output EGFR_dashboard.html

# XC50-only pairs with 3D RMSD
python cliff_dashboard.py \
    --input all_cliffs.csv \
    --cliff-type XC50 --top-n 200 \
    --output dashboard_xc50.html
```

### CLI flags

| Flag | Default | Description |
|---|---|---|
| `--input FILE` | required | Cliff CSV (scanner or enricher output) |
| `--output FILE` | `<input>_dashboard.html` | Output HTML |
| `--top-n INT` | 200 | Pairs to render (0 = all; >500 makes large files) |
| `--min-delta FLOAT` | 2.0 | Minimum \|ΔpXC50\| (NaN-delta pct_fallback pairs are kept) |
| `--target CHEMBL_ID` | all | Filter to one target |
| `--cliff-type` | `all` | `XC50`, `pct_fallback`, or `all` |
| `--skip-3d` | off | Skip ETKDGv3 + MMFF conformer generation and RMSD |

### RMSD interpretation

| RMSD (Å) | Interpretation |
|---|---|
| < 0.5 | Scaffold is rigid; cliff is likely **electronic** (changed atom alters charge, H-bond, polarisability) |
| 0.5–1.5 | Mixed — electronic + minor conformational perturbation |
| > 1.5 | Structural change **perturbs the preferred conformation** — possible binding-mode change |

---

## Recommended workflow

```bash
# Step 1: scan with assay-type filter
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --same-assay-type \
    --output all_cliffs.csv

# Step 2: enrich
python cliff_enrich.py \
    --input all_cliffs.csv \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --output all_cliffs_enriched.csv

# Step 3: dashboard — XC50 pairs, top 200 with RMSD
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
| `inactive_pct_activity` | float | % inhibition at ~10 µM; populated only for pct_fallback |
| `delta_pXC50` | float | \|active_pXC50 − inactive_pXC50\|; NaN for pct_fallback |
| `assay_chembl_id` | str | ChEMBL assay ID |
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

| Column | Notes |
|---|---|
| `tanimoto_ecfp4` | ECFP4 Tanimoto similarity (4 d.p.) |
| `active_assay_chembl_id` | Assay of the active measurement |
| `inactive_assay_chembl_id` | Assay of the inactive measurement; may differ for pct_fallback |
| `protein_class_l1` | Broad class: Enzyme, Ion channel, GPCR, Transcription factor… |
| `protein_class_l2` | e.g. Kinase, Phosphatase, Protease |
| `protein_class_l3` | e.g. Protein Kinase, Serine Protease |
| `protein_family` | Full ChEMBL pref_name or UniProt SIMILARITY comment |

---

## Scientific notes

### Why Tanimoto ≥ 0.95?

At this threshold, pairs typically differ by a **single atom or small functional group** (F→Cl, H→CH₃, CH₂ insertion, ring-size change). This makes structural differences unambiguous and directly interpretable. Lowering to 0.85 admits scaffold hops where the activity difference may not be attributable to one specific change.

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

These pairs have a quantified pXC50 for the active compound but only a single-point % inhibition measurement for the inactive. They are useful for:
- Identifying active compounds with clear structural neighbours that failed to show any activity at 10 µM
- High-throughput enumeration of potential cliffs from HTS data

They are scientifically weaker than XC50–XC50 pairs because the "inactivity" may reflect insufficient dose rather than true structural failure. Use `--xc50-only` for high-confidence SAR analysis.

### Salt and duplicate handling

Before fingerprint computation, molecules are deduplicated per assay using three passes:

1. Same `molecule_chembl_id` + assay → keep most potent
2. Same raw InChIKey + assay → keep most potent  
3. Same **parent** InChIKey (largest fragment + neutralised) + assay → keep most potent

This collapses free-base / hydrochloride / sodium salt forms of the same compound to a single entry, preventing spurious "cliffs" from salt differences.
