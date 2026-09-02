# Activity Cliff Scanner

A three-script pipeline for detecting, enriching, and visualising **activity cliffs** in ChEMBL — structurally similar compound pairs with large potency differences — using ECFP4 fingerprints, RDKit, and the ChEMBL 37 SQLite database.

Final output is an helpful dashboard ![dashboard](https://github.com/agiani99/Activity_Cliff_scanner/blob/main/Screenshot.png).

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
| `--cross-assay` | off | Find XC50–XC50 pairs **across assays** on the same target. Uses each molecule's best pXC50 from any assay. Essential for well-characterised targets (e.g. EGFR, CA-II) where IC50s are spread across hundreds of assays from different papers. Without this flag, both compounds must be in the exact same assay — rarely satisfied. |
| `--min-confidence` | `8` | Minimum ChEMBL assay confidence score. Default 8 = direct assay, single protein. Use `--min-confidence 7` when the target is a protein complex or heterodimer (most assays receive score 7 rather than 8). |
| `--relax-quality` | off | Skip the `potential_duplicate = 0` filter. ChEMBL sometimes over-flags replicated multi-lab measurements as duplicates. The `data_validity_comment` outlier filter is still applied. |
| `--chain-only` | off | **Biomimetics mode**: only consecutive pXC50-sorted pairs per assay. Eliminates star topology. Use with `--delta-pxic50 < 1.5`. |
| `--best-partner-only` | off | Keep only the highest-SALI inactive partner per active compound per assay. Lighter than `--chain-only`. |
| `--min-pxc50` | none | Exclude compounds weaker than this pXC50 (e.g. `--min-pxc50 5.0` drops IC50 > 10 µM). |

### Assay-type strictness hierarchy

| Flag | Filter applied | Expected pairs |
|---|---|---|
| *(none)* | None — any cross-assay pair | Most (268 K+ in full proteome run) |
| `--same-assay-type` | B↔B or F↔F only | Subset — removes cross-biology pairs |
| `--xc50-only` | Same assay, XC50 only | Fewest — highest confidence |

### Biomimetics vs activity cliffs — pair-reduction modes

At low `--delta-pxic50` (e.g. 0.5) a focused series of N compounds generates up to
N(N−1)/2 pairs — each "active" compound paired with every less-active structural
neighbor (star topology). Three mechanisms control this:

| Mode | Flag | Mechanism | Pairs produced |
|---|---|---|---|
| Unrestricted | *(none)* | All pairwise combinations | O(N²) |
| Best-partner | `--best-partner-only` | Keep highest-SALI partner per active | ≤ N |
| Chain | `--chain-only` | Consecutive pXC50 neighbours only | ≤ N−1 |

**SALI** (Structure-Activity Landscape Index) = `|ΔpXC50| / (1 − Tanimoto)` is added as
a column to every output. Use it to distinguish:
- SALI > 10 → classic activity cliff (large potency jump, very similar structure)
- SALI 1–10 → biomimetic / gradual SAR region (small change gives modest gain)

**Recommended biomimetics run:**
```bash
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --target CHEMBL2608 \
    --tanimoto 0.85 \
    --delta-pxic50 0.5 \
    --chain-only \
    --min-pxc50 5.0 \
    --output CHEMBL2608_biomimetics.csv
```

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

| Score | Meaning | Default included |
|---|---|---|
| 9 | Direct assay, single protein, exact compound tested | ✓ |
| 8 | Direct assay, single protein | ✓ |
| **7** | **Direct assay, protein complex / heterodimer** | use `--min-confidence 7` |
| ≤ 6 | Inferred / homologous / family-level target | ✗ |

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

New columns visible in the table and detail modal when present in the input CSV:
- **SALI** — sortable column showing Structure-Activity Landscape Index per pair
- **Raw XC50** — original concentration (nM or µM) shown under pXC50 in each molecule cell
- **Measurement type** — IC50 / EC50 / Ki / Kd shown as small label

### RMSD interpretation

| RMSD (Å) | Interpretation |
|---|---|
| < 0.5 | Scaffold rigid; cliff is likely **electronic** (changed atom alters charge, H-bond, polarisability) |
| 0.5–1.5 | Mixed — electronic + minor conformational perturbation |
| > 1.5 | Structural change **perturbs the preferred conformation** — possible binding-mode change |

---


## Cross-assay XC50 pairing — `--cross-assay`

### The problem: pct_fallback dominates without this flag

The scanner's default XC50–XC50 mode requires **both compounds to be in the exact same ChEMBL assay**.  
For well-characterised targets, ChEMBL organises IC50 data like this:

```
Paper 1 → Assay CHEMBL_A:  Cpd1 (5 nM),  Cpd2 (50 nM)
Paper 2 → Assay CHEMBL_B:  Cpd3 (8 nM),  Cpd4 (500 nM)
Paper 3 → Assay CHEMBL_C:  Cpd5 (12 nM), Cpd1 (6 nM)
HTS     → Assay CHEMBL_D:  All compounds → % inhibition
```

| Mode | Grouping key | Pairs found |
|---|---|---|
| XC50–XC50 (default) | `assay_chembl_id` | Only pairs where BOTH compounds are in the **same assay** (often 0–2 per target) |
| pct_fallback | `target_chembl_id` | Any XC50 compound vs any %inh compound from **any assay on the target** |
| XC50–XC50 `--cross-assay` | `target_chembl_id` | Best pXC50 per molecule per target, compared across all assays |

Result without `--cross-assay`: the dashboard shows **only pct_fallback pairs** because same-assay XC50 grouping finds almost nothing for targets with data spread across many publications.

### How `--cross-assay` works

1. For each *(molecule, target)* pair, keep only the single **best pXC50** value across all assays.
2. Group all molecules by target (instead of by assay).
3. Run the normal Tanimoto + ΔpXC50 cliff detection across the full compound set.
4. The original `assay_chembl_id` is preserved in the output row so the measurement provenance is always known.

This mirrors how a medicinal chemist reads an SAR table — they compare the best reported IC50 for each compound, regardless of which paper or lab produced it.

### When to use it

```bash
# Targets with IC50s spread across many assays (most well-characterised targets)
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --target CHEMBL1977 \
    --cross-assay \
    --min-confidence 7 \
    --output CHEMBL1977_cliffs.csv

# Combine with biomimetics mode
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --target CHEMBL1977 \
    --cross-assay \
    --tanimoto 0.85 --delta-pxic50 0.5 \
    --chain-only --min-pxc50 5.0 \
    --output CHEMBL1977_biomimetics.csv
```

### Diagnostic SQL — why are thousands of IC50s missing?

Save and run against your ChEMBL SQLite to pinpoint which filter is removing the most records:

```sql
-- 1. Total raw IC50s (no scanner filters)
SELECT COUNT(*) AS total_raw
FROM activities act
JOIN assays ass ON act.assay_id = ass.assay_id
JOIN target_dictionary td ON ass.tid = td.tid
WHERE td.chembl_id = 'CHEMBL1977' AND act.standard_type = 'IC50';

-- 2. Confidence score distribution (reveals if score < 8 holds most data)
SELECT ass.confidence_score, COUNT(*) AS n_ic50
FROM activities act
JOIN assays ass ON act.assay_id = ass.assay_id
JOIN target_dictionary td ON ass.tid = td.tid
WHERE td.chembl_id = 'CHEMBL1977' AND act.standard_type = 'IC50'
GROUP BY ass.confidence_score ORDER BY ass.confidence_score DESC;

-- 3. Relation distribution (reveals how many are >, <, =)
SELECT act.standard_relation, COUNT(*) AS n_ic50
FROM activities act
JOIN assays ass ON act.assay_id = ass.assay_id
JOIN target_dictionary td ON ass.tid = td.tid
WHERE td.chembl_id = 'CHEMBL1977' AND act.standard_type = 'IC50'
GROUP BY act.standard_relation ORDER BY n_ic50 DESC;

-- 4. After all scanner filters (what the scanner actually retrieves)
SELECT COUNT(*) AS scanner_retrieves
FROM activities act
JOIN assays ass ON act.assay_id = ass.assay_id
JOIN target_dictionary td ON ass.tid = td.tid
JOIN molecule_dictionary md ON act.molregno = md.molregno
JOIN compound_structures cs ON act.molregno = cs.molregno
WHERE td.chembl_id = 'CHEMBL1977'
  AND act.standard_type = 'IC50'
  AND ass.confidence_score >= 8
  AND act.standard_relation = '='
  AND act.standard_value IS NOT NULL
  AND cs.canonical_smiles IS NOT NULL
  AND act.standard_units IN ('nM','uM','µM','μM','pM','nmol/l','umol/l','mM','M')
  AND (act.potential_duplicate = 0 OR act.potential_duplicate IS NULL)
  AND (act.data_validity_comment IS NULL OR act.data_validity_comment = 'Manually validated');
```

The difference between query 1 and query 4 shows exactly which filter is responsible.  
Most common culprits:

| Filter | Impact |
|---|---|
| `confidence_score >= 8` | Protein complexes / heterodimers receive score 7 → use `--min-confidence 7` |
| `standard_relation = '='` | Removes `>10000 nM` (inactives) and `<1 nM` (very potent) — legitimate data |
| `potential_duplicate = 0` | ChEMBL over-flags multi-lab replicates → use `--relax-quality` |

## Recommended workflow

```bash
# Step 1a: scan — classic activity cliffs (cross-assay, default delta=2.0)
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --cross-assay \
    --same-assay-type \
    --output all_cliffs.csv

# Step 1b: biomimetics / low-delta scan on a specific target
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --target CHEMBL1977 \
    --cross-assay \
    --tanimoto 0.85 --delta-pxic50 0.5 \
    --chain-only --min-pxc50 5.0 \
    --min-confidence 7 \
    --output CHEMBL1977_biomimetics.csv

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
| `sali` | float | Structure-Activity Landscape Index = \|ΔpXC50\| / (1 − Tanimoto). High = true cliff; Low = biomimetic. |
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
| `sali` | Scanner | Structure-Activity Landscape Index = \|ΔpXC50\| / (1−Tanimoto). In output from scanner directly. |
| `tanimoto_ecfp4` | RDKit (enricher) | ECFP4 Tanimoto similarity (4 d.p.) — added by cliff_enrich.py |
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

---

## External data — `--external-csv`

Analyse proprietary or non-ChEMBL data without any database access.

### Requirements

| Column | Content | Notes |
|---|---|---|
| SMILES | Compound structure | Required |
| pChEMBL | −log₁₀(IC50/EC50/Ki/Kd in M) | Required — pChEMBL units, NOT raw nM |
| Compound ID | Molecule name / ID | Optional — auto-generated (`EXT_00001`…) if absent |
| Target / Assay | Grouping key | Optional but strongly recommended (see below) |

**The activity column must contain pChEMBL values (e.g. 8.3), not raw concentrations in nM.** Raw concentrations should be pre-converted: `pChEMBL = −log₁₀(IC50_M)`.

### Usage

```bash
# Auto-detect columns (tries common synonyms)
python activity_cliff_scanner.py \
    --external-csv my_data.csv \
    --output my_cliffs.csv

# Explicit column names (override auto-detection)
python activity_cliff_scanner.py \
    --external-csv my_data.csv \
    --smiles-col   "SMILES" \
    --activity-col "pChEMBL" \
    --id-col       "MOLECULE NAME" \
    --target-col   "TARGET NAME" \
    --output my_cliffs.csv

# With assay-level grouping (finer-grained than target)
python activity_cliff_scanner.py \
    --external-csv my_data.csv \
    --smiles-col   "SMILES" \
    --activity-col "pIC50" \
    --id-col       "Cpd_ID" \
    --assay-col    "Assay_ID" \
    --target-col   "Gene" \
    --output my_cliffs.csv

# Combine with any scanner flag
python activity_cliff_scanner.py \
    --external-csv my_data.csv \
    --smiles-col   "SMILES" \
    --activity-col "pChEMBL" \
    --id-col       "MOLECULE NAME" \
    --target-col   "TARGET NAME" \
    --tanimoto 0.85 --delta-pxic50 1.5 \
    --xc50-only \
    --output my_cliffs.csv
```

### Compound grouping — critical for multi-target datasets

Cliff detection only compares compounds **within the same group**. The grouping key determines which compounds are eligible to form pairs:

| Flags provided | Grouping key | Behaviour |
|---|---|---|
| `--assay-col` | Assay column | Pairs within each assay; multiple assays per target allowed |
| `--target-col` only | **Target column** | Pairs within each target; cross-target pairs excluded |
| Neither | Single `EXTERNAL_ASSAY` | All compounds compared against all — only correct for single-target datasets |

> **If your CSV contains multiple targets, always provide `--target-col`.**  
> Without it, all compounds land in one virtual assay and cross-target pairs are incorrectly generated.

### Auto-detected column name synonyms

| Field | Recognised names (case-insensitive) |
|---|---|
| SMILES | `smiles`, `canonical_smiles`, `structure`, `mol`, `molecule`, `smi` |
| pChEMBL | `pchembl_value`, `pxc50`, `pic50`, `pki`, `pkd`, `pec50`, `activity`, `potency`, `value` |
| Compound ID | `id`, `compound_id`, `molecule_id`, `name`, `compound_name`, `mol_id`, `cmpd_id` |
| Assay group | `assay_id`, `assay`, `group`, `series`, `project`, `batch`, `dataset` |
| Target | `target_id`, `target`, `protein`, `gene`, `receptor` |

### What happens internally

- pChEMBL is back-converted to nM (`standard_value = 10^(9 − pChEMBL)`) so the standard preprocessing pipeline works unchanged.
- All identity filters apply: InChIKey deduplication, parent InChIKey (salt stripping), Tanimoto = 1.0 guard.
- Extra columns in your CSV are preserved in the output with an `ext_` prefix.
- Compound IDs are auto-generated as `EXT_00001`, `EXT_00002`… if no ID column is found.
- Assay confidence score is set to 9 (trusted external data).
- No ChEMBL lookup is performed — `uniprot_id`, `target_name` etc. are populated from the target column or left as `EXTERNAL_TARGET`.


---

## Biomimetics and low-delta mode

When `--delta-pxic50` is lowered to 0.5–1.0, the scanner switches from finding
**activity cliffs** (large potency gaps) to finding **biomimetics** — structurally
very similar compounds with small but reproducible activity differences.  These
pairs reveal:

- Subtle pharmacophoric contributions of individual atoms or functional groups
- Bioisosteric replacements with measurable potency changes
- SAR gradients within a focused optimisation series

The tradeoff: at delta=0.5 and Tanimoto=0.85, one highly-active compound can
legitimately pair against 5–10 structural neighbours, producing a redundant
**star topology** where the same active compound is repeated at the centre of
many pairs.  Three new flags address this.

### New scanner flags

| Flag | Default | Description |
|---|---|---|
| `--chain-only` | off | **Biomimetics mode.** Within each assay, sort compounds by pXC50 descending and only report *consecutive* pairs. Compound A is compared only to its immediate lower-activity neighbour. Eliminates star-topology redundancy. Recommended with `--delta-pxic50 < 1.5`. |
| `--best-partner-only` | off | For each active compound keep only its single highest-SALI inactive partner. Less strict than `--chain-only`. Useful when multiple distinct structural series are present in one assay. |
| `--min-pxc50` | none | Exclude compounds below this pXC50 from the comparison. With delta=0.5 you would otherwise pair two very weak binders (e.g. pXC50=4.5 vs 4.0). `--min-pxc50 5.0` keeps only compounds with IC50 ≤ 10 µM. |

### SALI — Structure-Activity Landscape Index

Every pair now includes a `sali` column:

```
SALI = |ΔpXC50| / (1 − Tanimoto)
```

| SALI range | Interpretation |
|---|---|
| < 3 | **Biomimetic** — small structural change, small activity gain; valuable for SAR optimisation |
| 3–8 | Intermediate — notable activity change relative to structural similarity |
| > 8 | **Activity cliff** — large potency gap from a minor structural change; reveals critical pharmacophore |

SALI is displayed as a colour-coded badge in the dashboard:
- 🔵 Blue `< 3` — biomimetic
- 🟠 Orange `3–8` — intermediate
- 🔴 Red `> 8` — true cliff

### Recommended biomimetics run

```bash
python activity_cliff_scanner.py \
    --sqlite "C:/path/to/chembl_37_sqlite" \
    --target CHEMBL2608 \
    --tanimoto 0.85 \
    --delta-pxic50 0.5 \
    --chain-only \
    --min-pxc50 5.0 \
    --output CHEMBL2608_biomimetics.csv
```

### Topology comparison

| Mode | Pairs for 6 analogues | SAR picture |
|---|---|---|
| Default (delta=2.0) | 0–2 | Identifies hard cliffs only |
| Low delta, no flags (delta=0.5) | Up to 15 (star) | Redundant — same active repeated |
| `--chain-only` (delta=0.5) | ≤ 5 (ladder) | Clean step-by-step SAR |
| `--best-partner-only` (delta=0.5) | ≤ 6 (one per active) | Collapsed star — best partner retained |

