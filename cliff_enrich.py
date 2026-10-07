#!/usr/bin/env python3
"""
cliff_enrich.py  --  Post-hoc enrichment of activity_cliff_scanner output
==========================================================================
Adds four annotation layers to an existing all_cliffs.csv:

  a) tanimoto_ecfp4          ECFP4 Tanimoto similarity (4 d.p.)
                             Computed from SMILES with RDKit -- no network.
                             Vectorised: unique-SMILES fingerprint cache makes
                             this fast even for 268 K rows.

  b) uniprot_id              Already in the CSV from the scanner.  This script
                             fills any "N/A" gaps using the SQLite or ChEMBL
                             REST API.

  c) protein_class_l1        ChEMBL protein class hierarchy (Enzyme, Ion
     protein_class_l2        channel, GPCR, Kinase ...) from the SQLite
     protein_class_l3        protein_class table -- one batch query on the
     protein_family          unique set of target_chembl_ids.
                             Falls back to UniProt "SIMILARITY" comment when
                             SQLite is absent.

  d) active_assay_chembl_id  Renamed from the existing assay_chembl_id.
     inactive_assay_chembl_id
                             For XC50-XC50 cliff pairs: identical to
                             active_assay_chembl_id (both compounds come from
                             the same assay by definition of the scanner).
                             For pct_fallback pairs: the inactive was measured
                             in a different assay (% inhibition screen).  This
                             script recovers it with a SQLite lookup on the
                             inactive_chembl_id + target_chembl_id.  Without
                             SQLite, the column is set to "N/A (needs SQLite)".

DATA SOURCES (in order of preference)
--------------------------------------
  --sqlite PATH    Local ChEMBL .db file or directory -- recommended.
                   Required for: protein class, inactive assay IDs.
  REST API         ChEMBL + UniProt REST.  Used automatically when SQLite
                   is absent or a particular ID has no SQLite entry.
                   Needs www.ebi.ac.uk and rest.uniprot.org reachable.

USAGE
-----
  # Recommended: SQLite supplies everything except UniProt family text
  python cliff_enrich.py \\
      --input  all_cliffs.csv \\
      --sqlite "C:/path/to/chembl_37_sqlite" \\
      --output all_cliffs_enriched.csv

  # API-only (slower, needs network):
  python cliff_enrich.py \\
      --input  all_cliffs.csv \\
      --output all_cliffs_enriched.csv

  # Just Tanimoto + assay rename, skip network entirely:
  python cliff_enrich.py \\
      --input  all_cliffs.csv \\
      --skip-protein-class \\
      --output all_cliffs_enriched.csv

PERFORMANCE (268 K rows)
------------------------
  Tanimoto          :  ~60-120 s (vectorised, unique-SMILES cache)
  Protein class     :  <5 s    (one SQL batch on unique target IDs)
  UniProt family    :  ~1 s/ID (cached; typically <2000 unique IDs)
  Inactive assay ID :  <30 s   (one SQL batch on unique inactive IDs)
"""

import argparse
import glob
import logging
import os
import sqlite3
import sys
import time
from typing import Optional

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, DataStructs

# ---------------------------------------------------------------------------
# SALT / COUNTER-ION REMOVAL -- applied to every molecule parsed by RDKit
# ---------------------------------------------------------------------------
from rdkit.Chem.MolStandardize import rdMolStandardize as _rdMS
_DESALT_LFC = _rdMS.LargestFragmentChooser(preferOrganic=True)
_DESALT_UN  = _rdMS.Uncharger()


def desalt_mol(mol):
    """Parent molecule: largest organic fragment, neutralised.
    Strips counter-ions / solvates (HCl, Na+, TFA, water ...). None -> None."""
    if mol is None:
        return None
    try:
        if len(Chem.GetMolFrags(mol)) > 1:
            mol = _DESALT_LFC.choose(mol)
        return _DESALT_UN.uncharge(mol)
    except Exception:
        return mol


RDLogger.DisableLog("rdApp.*")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────────────
ECFP_RADIUS  = 2
ECFP_NBITS   = 2048
API_DELAY_S  = 0.08
CHEMBL_BASE  = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
RCSB_SEARCH  = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GQL     = "https://data.rcsb.org/graphql"

# ── EC-number prefix → (l1, l2, l3) — longest prefix wins ────────────────────
_EC_MAP = [
    ("2.7.10.1",  ("Enzyme", "Kinase", "Receptor Tyrosine Kinase")),
    ("2.7.10.2",  ("Enzyme", "Kinase", "Non-receptor Tyrosine Kinase")),
    ("2.7.11.1",  ("Enzyme", "Kinase", "Serine/Threonine Kinase")),
    ("2.7.11.2",  ("Enzyme", "Kinase", "Serine/Threonine Kinase")),
    ("2.7.12",    ("Enzyme", "Kinase", "Dual-specificity Kinase")),
    ("2.7.13",    ("Enzyme", "Kinase", "Histidine Kinase")),
    ("2.7.4",     ("Enzyme", "Kinase", "Nucleoside Phosphate Kinase")),
    ("2.7.1",     ("Enzyme", "Kinase", "Kinase")),
    ("2.7",       ("Enzyme", "Kinase", "Kinase")),
    ("3.4.21",    ("Enzyme", "Protease", "Serine Protease")),
    ("3.4.22",    ("Enzyme", "Protease", "Cysteine Protease")),
    ("3.4.23",    ("Enzyme", "Protease", "Aspartic Protease")),
    ("3.4.24",    ("Enzyme", "Protease", "Metalloprotease")),
    ("3.4.25",    ("Enzyme", "Protease", "Threonine Protease")),
    ("3.4",       ("Enzyme", "Protease", "Protease")),
    ("3.1.4.17",  ("Enzyme", "Phosphodiesterase", "Cyclic Nucleotide PDE")),
    ("3.1.4",     ("Enzyme", "Phosphodiesterase", "Phosphodiesterase")),
    ("3.1.3.48",  ("Enzyme", "Phosphatase", "Tyrosine Phosphatase")),
    ("3.1.3.16",  ("Enzyme", "Phosphatase", "Ser/Thr Phosphatase")),
    ("3.1.3",     ("Enzyme", "Phosphatase", "Phosphatase")),
    ("3.5.1.98",  ("Enzyme", "Deacetylase", "HDAC")),
    ("3.5.1",     ("Enzyme", "Hydrolase", "Amide Hydrolase")),
    ("2.1.1.356", ("Enzyme", "Methyltransferase", "Histone Methyltransferase")),
    ("2.1.1",     ("Enzyme", "Methyltransferase", "Methyltransferase")),
    ("2.3.1",     ("Enzyme", "Acetyltransferase", "Acetyltransferase")),
    ("3.6.5",     ("Enzyme", "GTPase", "GTPase")),
    ("3.6.1",     ("Enzyme", "ATPase", "ATPase")),
    ("4.6.1.1",   ("Enzyme", "Cyclase", "Adenylyl Cyclase")),
    ("4.6.1.2",   ("Enzyme", "Cyclase", "Guanylyl Cyclase")),
    ("4.6.1",     ("Enzyme", "Cyclase", "Nucleotidyl Cyclase")),
    ("1.14.14",   ("Enzyme", "Oxidoreductase", "Cytochrome P450")),
    ("1.14.13",   ("Enzyme", "Oxidoreductase", "Cytochrome P450")),
    ("1.14",      ("Enzyme", "Oxidoreductase", "Oxygenase")),
    ("1.1.1",     ("Enzyme", "Oxidoreductase", "Dehydrogenase")),
    ("1",         ("Enzyme", "Oxidoreductase", "Oxidoreductase")),
    ("2.4",       ("Enzyme", "Glycosyltransferase", "Glycosyltransferase")),
    ("5.2.1.8",   ("Enzyme", "Isomerase", "Peptidyl-Prolyl Isomerase (FKBP)")),
    ("2",         ("Enzyme", "Transferase", "Transferase")),
    ("3.2",       ("Enzyme", "Glycosylase", "Glycosylase")),
    ("3",         ("Enzyme", "Hydrolase", "Hydrolase")),
    ("4",         ("Enzyme", "Lyase", "Lyase")),
    ("5",         ("Enzyme", "Isomerase", "Isomerase")),
    ("6",         ("Enzyme", "Ligase", "Ligase")),
]

_PFAM_MAP = {
    "PF00001": ("Membrane receptor", "GPCR", "Family A GPCR"),
    "PF00002": ("Membrane receptor", "GPCR", "Family B GPCR"),
    "PF00003": ("Membrane receptor", "GPCR", "Family C GPCR"),
    "PF10292": ("Membrane receptor", "GPCR", "Family A GPCR"),
    "PF00104": ("Nuclear receptor",  "Nuclear Hormone Receptor", "NHR ligand-binding domain"),
    "PF00105": ("Nuclear receptor",  "Nuclear Hormone Receptor", "NHR DNA-binding domain"),
    "PF00520": ("Ion channel",       "Ion Channel", "Voltage-gated ion channel"),
    "PF07714": ("Enzyme",            "Kinase", "Protein Tyrosine Kinase"),
    "PF00069": ("Enzyme",            "Kinase", "Protein Kinase"),
    "PF13955": ("Enzyme",            "Kinase", "Protein Kinase"),
    "PF00089": ("Enzyme",            "Protease", "Serine Protease (Trypsin)"),
    "PF00082": ("Enzyme",            "Protease", "Serine Protease (Subtilisin)"),
    "PF00026": ("Enzyme",            "Protease", "Aspartic Protease"),
    "PF00112": ("Enzyme",            "Protease", "Cysteine Protease"),
    "PF01471": ("Enzyme",            "Protease", "Metalloprotease"),
    "PF00233": ("Enzyme",            "Phosphodiesterase", "Cyclic Nucleotide PDE"),
    "PF01663": ("Enzyme",            "Phosphodiesterase", "Phosphodiesterase"),
    "PF00102": ("Enzyme",            "Phosphatase", "Tyrosine Phosphatase"),
    "PF00782": ("Enzyme",            "Phosphatase", "Dual-specificity Phosphatase"),
    "PF02132": ("Enzyme",            "Deacetylase", "HDAC"),
    "PF00856": ("Enzyme",            "Methyltransferase", "SET domain Methyltransferase"),
    "PF00439": ("Enzyme",            "Bromodomain", "Bromodomain"),
    "PF00628": ("Enzyme",            "Bromodomain", "Bromodomain"),
    "PF00664": ("Transporter",       "ABC Transporter", "ABC transmembrane domain"),
    "PF00005": ("Transporter",       "ABC Transporter", "ABC ATPase domain"),
    "PF00179": ("Enzyme",            "Ubiquitin ligase", "E2 ubiquitin-conjugating enzyme"),
    "PF00167": ("Enzyme",            "Ubiquitin ligase", "RING domain E3 ligase"),
}


def _ec_to_class(ec: str):
    """Longest-prefix match of an EC number in _EC_MAP."""
    ec = ec.strip()
    for prefix, cls in _EC_MAP:
        if ec.startswith(prefix):
            return cls
    return None


# ── SIFTS classification file loaders ────────────────────────────────────────

def load_uniprot_pdb(path: str) -> dict:
    """
    Load uniprot_pdb.csv (compact UniProt->PDB summary).
    Format: SP_PRIMARY,PDB_IDS  where PDB_IDS is semicolon-separated.
    Returns {uniprot_accession -> [PDB_ID, ...]} (uppercase, sorted).
    """
    logger.info(f"  Loading UniProt->PDB compact map from {os.path.basename(path)}...")
    df = pd.read_csv(path, comment='#', low_memory=False)
    # Normalise column names (may be SP_PRIMARY,PDB_IDS or similar)
    df.columns = [c.strip().upper() for c in df.columns]
    uid_col = next((c for c in df.columns if 'PRIMARY' in c or 'UNIPROT' in c
                    or c in ('SP_PRIMARY', 'ACCESSION')), df.columns[0])
    pdb_col = next((c for c in df.columns if 'PDB' in c), df.columns[1])

    result: dict = {}
    for _, row in df.iterrows():
        uid  = str(row[uid_col]).strip()
        pdbs = str(row[pdb_col]).strip()
        if uid and uid != 'nan' and pdbs and pdbs != 'nan':
            result[uid] = sorted(p.upper() for p in pdbs.split(';') if p.strip())
    logger.info(f"    {len(result):,} UniProt accessions loaded")
    return result


def load_pdb_enzyme(path: str) -> dict:
    """
    Load pdb_chain_enzyme.csv  (EC numbers per PDB chain).
    Returns {PDB_ID_upper -> set(ec_numbers)}.
    """
    logger.info(f"  Loading PDB->EC from {os.path.basename(path)}...")
    df = pd.read_csv(path, comment='#', low_memory=False)
    df.columns = [c.strip().upper() for c in df.columns]
    pdb_col = next(c for c in df.columns if c == 'PDB')
    ec_col  = next(c for c in df.columns if 'ACCESSION' in c or 'EC' in c)
    result: dict = {}
    for _, row in df.iterrows():
        pdb = str(row[pdb_col]).strip().upper()
        ec  = str(row[ec_col]).strip()
        if pdb and ec and ec != 'nan':
            result.setdefault(pdb, set()).add(ec)
    logger.info(f"    {len(result):,} PDB entries with EC numbers")
    return result


def load_pdb_pfam(path: str) -> dict:
    """
    Load pdb_chain_pfam.csv  (Pfam domain IDs per PDB chain).
    Returns {PDB_ID_upper -> set(pfam_ids)}.
    """
    logger.info(f"  Loading PDB->Pfam from {os.path.basename(path)}...")
    df = pd.read_csv(path, comment='#', low_memory=False)
    df.columns = [c.strip().upper() for c in df.columns]
    pdb_col  = next(c for c in df.columns if c == 'PDB')
    pfam_col = next(c for c in df.columns if 'PFAM' in c or 'ACCESSION' in c)
    result: dict = {}
    for _, row in df.iterrows():
        pdb  = str(row[pdb_col]).strip().upper()
        pfam = str(row[pfam_col]).strip().upper()
        if pdb and pfam and pfam != 'NAN':
            result.setdefault(pdb, set()).add(pfam)
    logger.info(f"    {len(result):,} PDB entries with Pfam domains")
    return result


def classify_uniprot_from_sifts(
    uniprot_id: str,
    uniprot_pdb: dict,
    pdb_enzyme:  dict,
    pdb_pfam:    dict,
) -> dict:
    """
    Derive protein_class_l1/l2/l3 and protein_family for a UniProt accession
    using only local SIFTS files — no network call required.

    Priority:
      1. EC number (most specific for enzymes — kinase, protease, PDE, ...)
      2. Pfam domain (covers non-enzymes: GPCR, NHR, ion channel, ...)
      3. Returns all-N/A when neither is found (target has no PDB structure
         or its PDB structure lacks EC/Pfam annotation).

    For pipe-separated multi-accession cells (e.g. 'P00533|Q9UJ70'),
    all accessions are tried and the first hit is used.
    """
    _na = {"protein_class_l1": "N/A", "protein_class_l2": "N/A",
           "protein_class_l3": "N/A", "protein_family":   "N/A"}

    accessions = [a.strip() for a in str(uniprot_id).split("|")
                  if a.strip() and a.strip() != "N/A"]
    if not accessions:
        return _na

    # Collect all PDB IDs for this target
    pdb_ids: set = set()
    for acc in accessions:
        pdb_ids.update(uniprot_pdb.get(acc, []))

    if not pdb_ids:
        return _na

    # 1. Try EC numbers
    ec_hits: list = []
    for pdb in pdb_ids:
        for ec in pdb_enzyme.get(pdb, set()):
            cls = _ec_to_class(ec)
            if cls:
                ec_hits.append(cls)

    if ec_hits:
        # Most specific = longest l3 label (breaks ties toward sub-class)
        best = max(ec_hits, key=lambda c: len(c[2]))
        return {"protein_class_l1": best[0], "protein_class_l2": best[1],
                "protein_class_l3": best[2], "protein_family":   best[1]}

    # 2. Try Pfam domains
    pfam_hits: list = []
    for pdb in pdb_ids:
        for pfam in pdb_pfam.get(pdb, set()):
            cls = _PFAM_MAP.get(pfam)
            if cls:
                pfam_hits.append(cls)

    if pfam_hits:
        best = max(pfam_hits, key=lambda c: len(c[2]))
        return {"protein_class_l1": best[0], "protein_class_l2": best[1],
                "protein_class_l3": best[2], "protein_family":   best[1]}

    return _na


def annotate_protein_class_sifts(
    df: pd.DataFrame,
    uniprot_pdb: dict,
    pdb_enzyme:  dict,
    pdb_pfam:    dict,
) -> pd.DataFrame:
    """Apply SIFTS-based protein class annotation to every row."""
    logger.info("  Classifying targets via SIFTS EC + Pfam (fully offline)...")
    unique_uids = df["uniprot_id"].dropna().unique().tolist()
    uid_cache: dict = {}
    for uid in tqdm(unique_uids, desc="SIFTS classify", unit="target"):
        uid_cache[uid] = classify_uniprot_from_sifts(
            uid, uniprot_pdb, pdb_enzyme, pdb_pfam
        )
    for field in ("protein_class_l1", "protein_class_l2",
                  "protein_class_l3", "protein_family"):
        df[field] = df["uniprot_id"].map(
            lambda u, f=field: uid_cache.get(u, {}).get(f, "N/A")
        )
    n_l2 = (df["protein_class_l2"] != "N/A").sum()
    logger.info(f"  protein_class_l2 coverage: {n_l2:,} / {len(df):,} rows "
                f"({100*n_l2/len(df):.1f}%)")
    return df

_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=ECFP_RADIUS, fpSize=ECFP_NBITS)
_SESSION = requests.Session()
_SESSION.headers.update({"Accept": "application/json",
                          "User-Agent": "cliff-enrich/1.0"})


# ═════════════════════════════════════════════════════════════════════════════
# (a) TANIMOTO ECFP4
# ═════════════════════════════════════════════════════════════════════════════

def compute_tanimoto_column(df: pd.DataFrame) -> pd.Series:
    """
    Compute ECFP4 Tanimoto similarity for every row.

    Uses a fingerprint cache keyed by canonical SMILES so that the same
    compound (appearing in many pairs) is only parsed once.
    For 268 K rows with ~20 K unique active SMILES this saves ~90 % of
    RDKit calls.
    """
    logger.info("Computing ECFP4 Tanimoto similarities...")

    # Collect all unique SMILES needing fingerprints
    all_smiles = set(df["active_smiles"].dropna()) | set(df["inactive_smiles"].dropna())
    logger.info(f"  {len(all_smiles):,} unique SMILES -> fingerprint cache")

    fp_cache: dict = {}
    for smi in tqdm(all_smiles, desc="Fingerprints", unit="mol"):
        mol = desalt_mol(Chem.MolFromSmiles(str(smi).strip())) if isinstance(smi, str) else None
        fp_cache[smi] = _GEN.GetFingerprint(mol) if mol else None

    def pair_tanimoto(row) -> float:
        fa = fp_cache.get(row["active_smiles"])
        fb = fp_cache.get(row["inactive_smiles"])
        if fa is None or fb is None:
            return float("nan")
        return round(DataStructs.TanimotoSimilarity(fa, fb), 4)

    logger.info("  Computing pairwise Tanimoto...")
    tans = df.apply(pair_tanimoto, axis=1)
    logger.info(f"  Done. mean={tans.mean():.4f}  min={tans.min():.4f}  max={tans.max():.4f}")
    return tans


# ═════════════════════════════════════════════════════════════════════════════
# SQLite HELPERS
# ═════════════════════════════════════════════════════════════════════════════

def find_sqlite_db(path: str) -> str:
    path = os.path.normpath(path)
    if os.path.isfile(path) and path.endswith(".db"):
        return path
    if os.path.isdir(path):
        candidates = (glob.glob(os.path.join(path, "*.db")) +
                      glob.glob(os.path.join(path, "**", "*.db"), recursive=False))
        candidates = [c for c in candidates if os.path.isfile(c)]
        if len(candidates) == 1:
            logger.info(f"  Found SQLite: {candidates[0]}")
            return candidates[0]
        if len(candidates) > 1:
            raise FileNotFoundError(
                f"Multiple .db files in {path}. Specify the exact file.")
    raise FileNotFoundError(f"No SQLite .db found at: {path}")


def _sqlite_conn(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA cache_size   = -65536")   # 64 MB
    conn.execute("PRAGMA temp_store   = MEMORY")
    conn.execute("PRAGMA mmap_size    = 1073741824")
    conn.execute("PRAGMA synchronous  = OFF")
    return conn


# ═════════════════════════════════════════════════════════════════════════════
# (c) PROTEIN CLASS
# ═════════════════════════════════════════════════════════════════════════════

_SQL_PROTEIN_CLASS = """
SELECT
    td.chembl_id                   AS target_chembl_id,
    pc.l1                          AS protein_class_l1,
    pc.l2                          AS protein_class_l2,
    pc.l3                          AS protein_class_l3,
    pc.pref_name                   AS protein_family
FROM target_dictionary td
JOIN target_components tc   ON td.tid          = tc.tid
JOIN component_class   cc   ON tc.component_id = cc.component_id
JOIN protein_class     pc   ON cc.protein_class_id = pc.protein_class_id
WHERE td.chembl_id IN ({ph})
GROUP BY td.chembl_id
"""

_SQL_TARGET_TYPE_FALLBACK = """
SELECT chembl_id AS target_chembl_id, target_type
FROM target_dictionary
WHERE chembl_id IN ({ph})
"""


def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    cur = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?", (name,)
    )
    return cur.fetchone() is not None


def _list_tables(conn: sqlite3.Connection) -> list:
    return [r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
    ).fetchall()]


def fetch_protein_class_sqlite(db_path: str, target_ids: list) -> dict:
    """
    Returns {target_chembl_id: {protein_class_l1, l2, l3, protein_family}}.

    Tries three approaches in order:
      1. protein_class table (full hierarchy — present in most ChEMBL 37 builds)
      2. target_class table  (denormalized table added in some builds)
      3. target_type column  (always present in target_dictionary — broad fallback)
    """
    if not target_ids:
        return {}
    db_path = find_sqlite_db(db_path)
    conn = _sqlite_conn(db_path)

    try:
        available = _list_tables(conn)
        ph = ", ".join(["?"] * len(target_ids))

        # ── Strategy 1: protein_class table ───────────────────────────────
        if "protein_class" in available and "component_class" in available:
            logger.info("  Using protein_class table (full hierarchy)")
            sql = _SQL_PROTEIN_CLASS.format(ph=ph)
            try:
                rows = conn.execute(sql, target_ids).fetchall()
                if rows:
                    return {
                        row["target_chembl_id"]: {
                            "protein_class_l1": row["protein_class_l1"] or "N/A",
                            "protein_class_l2": row["protein_class_l2"] or "N/A",
                            "protein_class_l3": row["protein_class_l3"] or "N/A",
                            "protein_family":   row["protein_family"]   or "N/A",
                        }
                        for row in rows
                    }
            except sqlite3.OperationalError as exc:
                logger.warning(f"  protein_class query failed: {exc}")

        # ── Strategy 2: target_class table (denormalized) ─────────────────
        if "target_class" in available:
            logger.info("  Using target_class table (denormalized)")
            sql_tc = """
                SELECT td.chembl_id AS target_chembl_id,
                       tc.l1, tc.l2, tc.l3, tc.l4
                FROM target_dictionary td
                JOIN target_class tc ON td.tid = tc.tid
                WHERE td.chembl_id IN ({ph})
                GROUP BY td.chembl_id
            """.format(ph=ph)
            try:
                rows = conn.execute(sql_tc, target_ids).fetchall()
                if rows:
                    return {
                        row["target_chembl_id"]: {
                            "protein_class_l1": dict(row).get("l1") or "N/A",
                            "protein_class_l2": dict(row).get("l2") or "N/A",
                            "protein_class_l3": dict(row).get("l3") or "N/A",
                            "protein_family":   dict(row).get("l4") or "N/A",
                        }
                        for row in rows
                    }
            except sqlite3.OperationalError as exc:
                logger.warning(f"  target_class query failed: {exc}")

        # ── Strategy 3: target_type from target_dictionary (always works) ──
        logger.info(
            "  protein_class / target_class tables not found.\n"
            f"  Available tables: {', '.join(available)}\n"
            "  Falling back to target_dictionary.target_type (broad classification only).\n"
            "  Detailed protein family will be filled from UniProt REST API."
        )
        sql_tt = _SQL_TARGET_TYPE_FALLBACK.format(ph=ph)
        rows = conn.execute(sql_tt, target_ids).fetchall()
        return {
            row["target_chembl_id"]: {
                "protein_class_l1": row["target_type"] or "N/A",
                "protein_class_l2": "N/A",
                "protein_class_l3": "N/A",
                "protein_family":   "N/A",   # filled later by UniProt
            }
            for row in rows
        }

    finally:
        conn.close()


def fetch_protein_family_uniprot(uniprot_ids: list) -> dict:
    """
    Query UniProt REST API for protein family information.
    Returns {uniprot_id: family_string}.
    Extracts the "SIMILARITY" comment which contains text like
    "Belongs to the protein kinase superfamily. Ser/Thr protein kinase family."
    """
    if not uniprot_ids:
        return {}
    result: dict = {}
    unique = [u for u in set(uniprot_ids) if u and u != "N/A"]
    logger.info(f"  Fetching UniProt family for {len(unique):,} unique accessions...")
    for uid in tqdm(unique, desc="UniProt", unit="acc"):
        # Handle pipe-separated multi-accession entries
        accession = str(uid).split("|")[0].strip()
        if not accession or accession == "N/A":
            result[uid] = "N/A"
            continue
        try:
            r = _SESSION.get(
                f"{UNIPROT_BASE}/{accession}?format=json",
                timeout=15,
            )
            if r.status_code != 200:
                result[uid] = "N/A"
                time.sleep(API_DELAY_S)
                continue
            data = r.json()
            family = "N/A"
            for comment in data.get("comments", []):
                if comment.get("commentType") == "SIMILARITY":
                    texts = comment.get("texts", [])
                    if texts:
                        family = texts[0].get("value", "N/A")
                    break
            result[uid] = family
            time.sleep(API_DELAY_S)
        except Exception as exc:
            logger.debug(f"    UniProt {uid}: {exc}")
            result[uid] = "N/A"
    return result


def _fetch_protein_class_by_id(pc_id: int, cache: dict) -> dict:
    """
    Fetch one protein_class record by ID from ChEMBL REST.
    Results are cached in `cache` to avoid duplicate calls.
    Returns {protein_class_l1, l2, l3, protein_family} or {}.
    """
    if pc_id in cache:
        return cache[pc_id]
    try:
        r = _SESSION.get(f"{CHEMBL_BASE}/protein_class/{pc_id}.json", timeout=15)
        if r.status_code == 200:
            d = r.json()
            result = {
                "protein_class_l1": d.get("l1") or "N/A",
                "protein_class_l2": d.get("l2") or "N/A",
                "protein_class_l3": d.get("l3") or "N/A",
                "protein_family":   d.get("pref_name") or d.get("short_name") or "N/A",
            }
        else:
            result = {}
        time.sleep(API_DELAY_S)
    except Exception:
        result = {}
    cache[pc_id] = result
    return result


def fetch_chembl_target_class_api(target_ids: list) -> dict:
    """
    Two-step ChEMBL REST protein classification.

    Step 1 — /target/{id}.json:
        Reads protein_classifications[].protein_classification_id.
        NOTE: l1/l2/l3 are NOT in the target response — only id + desc.

    Step 2 — /protein_class/{id}.json (one call per UNIQUE class ID):
        Fetches l1, l2, l3 for each unique protein_classification_id.
        Results are cached so 783 targets sharing ~100 class IDs
        only trigger ~100 extra calls, not 783.

    L2 gives the actionable label: Kinase, Protease, Phosphodiesterase,
    Family A G protein-coupled receptor, Nuclear receptor, Ion channel...
    """
    if not target_ids:
        return {}
    unique_targets = list(set(target_ids))
    logger.info(f"  Fetching ChEMBL target classification for {len(unique_targets):,} targets...")

    pc_cache: dict = {}          # {protein_class_id: {l1, l2, l3, protein_family}}
    result:   dict = {}

    for tid in tqdm(unique_targets, desc="ChEMBL target class", unit="tgt"):
        try:
            r = _SESSION.get(f"{CHEMBL_BASE}/target/{tid}.json", timeout=15)
            if r.status_code != 200:
                result[tid] = {}
                time.sleep(API_DELAY_S)
                continue
            data = r.json()

            cls_result: dict = {}
            for comp in data.get("target_components", []):
                for cls in comp.get("protein_classifications", []):
                    pc_id = cls.get("protein_classification_id")
                    desc  = cls.get("protein_classification_desc", "")

                    if pc_id is not None:
                        # Fetch and cache the full hierarchy for this class ID
                        pc_data = _fetch_protein_class_by_id(pc_id, pc_cache)
                        if pc_data:
                            cls_result = pc_data
                            break

                    # Fallback: use desc only when ID lookup fails
                    if desc and not cls_result:
                        cls_result = {
                            "protein_class_l1": "N/A",
                            "protein_class_l2": "N/A",
                            "protein_class_l3": "N/A",
                            "protein_family":   desc,
                        }
                if cls_result:
                    break   # first component with a valid classification
            result[tid] = cls_result
            time.sleep(API_DELAY_S)

        except Exception as exc:
            logger.debug(f"    {tid}: {exc}")
            result[tid] = {}

    n_resolved = sum(1 for v in result.values()
                     if v.get("protein_class_l2", "N/A") != "N/A")
    logger.info(f"  Resolved L2 class for {n_resolved:,} / {len(unique_targets):,} targets")
    logger.info(f"  Unique protein class IDs fetched: {len(pc_cache):,}")
    return result


# ═════════════════════════════════════════════════════════════════════════════
# (b) UNIPROT ID FILL  (gaps only)
# ═════════════════════════════════════════════════════════════════════════════

_SQL_UNIPROT_FILL = """
SELECT
    td.chembl_id                    AS target_chembl_id,
    GROUP_CONCAT(DISTINCT cs.accession) AS uniprot_ids
FROM target_dictionary td
JOIN target_components  tc  ON td.tid          = tc.tid
JOIN component_sequences cs ON tc.component_id = cs.component_id
WHERE td.chembl_id IN ({ph})
GROUP BY td.chembl_id
"""


def fill_uniprot_sqlite(db_path: str, target_ids: list) -> dict:
    if not target_ids:
        return {}
    db_path = find_sqlite_db(db_path)
    ph  = ", ".join(["?"] * len(target_ids))
    sql = _SQL_UNIPROT_FILL.format(ph=ph)
    conn = _sqlite_conn(db_path)
    try:
        rows = conn.execute(sql, target_ids).fetchall()
    finally:
        conn.close()
    return {row["target_chembl_id"]: row["uniprot_ids"] or "N/A" for row in rows}


# ═════════════════════════════════════════════════════════════════════════════
# (d) INACTIVE ASSAY ID  (pct_fallback pairs)
# ═════════════════════════════════════════════════════════════════════════════

_SQL_INACTIVE_ASSAY = """
SELECT DISTINCT
    md.chembl_id    AS molecule_chembl_id,
    td.chembl_id    AS target_chembl_id,
    ass.chembl_id   AS assay_chembl_id
FROM activities       act
JOIN assays           ass ON act.assay_id   = ass.assay_id
JOIN target_dictionary td  ON ass.tid        = td.tid
JOIN molecule_dictionary md ON act.molregno  = md.molregno
WHERE md.chembl_id  IN ({ph_mol})
  AND td.chembl_id  IN ({ph_tgt})
  AND act.standard_type  IN ('Inhibition', 'Activity')
  AND act.standard_units = '%'
  AND ass.confidence_score >= 8
"""


def fetch_inactive_assay_sqlite(db_path: str,
                                inactive_ids: list,
                                target_ids: list) -> dict:
    """
    Returns {(inactive_chembl_id, target_chembl_id): assay_chembl_id}.
    """
    if not inactive_ids or not target_ids:
        return {}
    db_path = find_sqlite_db(db_path)
    ph_m = ", ".join(["?"] * len(inactive_ids))
    ph_t = ", ".join(["?"] * len(target_ids))
    sql  = _SQL_INACTIVE_ASSAY.format(ph_mol=ph_m, ph_tgt=ph_t)
    params = inactive_ids + target_ids
    conn = _sqlite_conn(db_path)
    try:
        rows = conn.execute(sql, params).fetchall()
    finally:
        conn.close()
    # Keep first assay found per (molecule, target) pair
    result: dict = {}
    for row in rows:
        key = (row["molecule_chembl_id"], row["target_chembl_id"])
        if key not in result:
            result[key] = row["assay_chembl_id"]
    return result


# ═════════════════════════════════════════════════════════════════════════════
# MAIN ENRICHMENT
# ═════════════════════════════════════════════════════════════════════════════


# ═════════════════════════════════════════════════════════════════════════════
# SIFTS-BASED PDB ANNOTATION
# ═════════════════════════════════════════════════════════════════════════════

def load_sifts(path: str) -> dict:
    """
    Load the SIFTS pdb_chain_uniprot.csv file and return a dict:
        {uniprot_accession -> [PDB_ID, ...]}  (sorted, uppercase, deduplicated)

    The file has a comment header line starting with '#' followed by a
    CSV with columns: PDB, CHAIN, SP_PRIMARY, RES_BEG, RES_END, ...
    SP_PRIMARY is the UniProt/SwissProt accession.

    Multi-chain entries (same PDB, multiple chains / segments) are collapsed
    to a single PDB ID per UniProt.
    """
    logger.info(f"  Loading SIFTS mapping from {path} ...")
    df = pd.read_csv(path, comment='#', low_memory=False)
    df = df.dropna(subset=['SP_PRIMARY', 'PDB'])
    df['PDB'] = df['PDB'].str.upper().str.strip()

    sifts = (
        df.groupby('SP_PRIMARY')['PDB']
          .apply(lambda x: sorted(set(x)))
          .to_dict()
    )
    n_pdb = sum(len(v) for v in sifts.values())
    logger.info(f"  SIFTS: {len(sifts):,} UniProt accessions, {n_pdb:,} PDB entries total")
    return sifts


def _resolve_sifts_for_row(uniprot_cell: str, sifts: dict) -> list:
    """
    Handle pipe-separated multi-accession uniprot_id cells (e.g. 'P00533|Q15303')
    and return the merged, sorted PDB list for all accessions.
    """
    if not uniprot_cell or uniprot_cell == "N/A":
        return []
    accessions = [a.strip() for a in str(uniprot_cell).split("|") if a.strip()]
    pdbs = set()
    for acc in accessions:
        pdbs.update(sifts.get(acc, []))
    return sorted(pdbs)


def find_compound_in_target_pdb(
    inchikey: str,
    uniprot_id: str,
    target_pdb_ids: list,
) -> dict:
    """
    Search RCSB for PDB structures where:
      (a) the ligand matches `inchikey` exactly (graph-exact chemical search), AND
      (b) the protein chain is annotated to `uniprot_id`.

    Using both filters prevents false positives (same ligand in a different
    protein target) and dramatically reduces the search space compared to
    a global InChIKey search.

    Returns {"pdb_id": str, "resolution_A": float|None, "method": str}
    or {} if not found / network unavailable.
    """
    if not inchikey or not uniprot_id or uniprot_id == "N/A":
        return {}

    # Restrict UniProt filter to accessions actually in SIFTS for this target
    accessions = [a.strip() for a in str(uniprot_id).split("|") if a.strip()]

    body = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "chemical",
                    "parameters": {
                        "value": inchikey,
                        "type": "descriptor",
                        "descriptor_type": "InChI",
                        "match_type": "graph-exact",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": (
                            "rcsb_polymer_entity_container_identifiers"
                            ".reference_sequence_identifiers.database_accession"
                        ),
                        "operator": "in",
                        "value": accessions,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 5}},
    }

    try:
        r = requests.post(RCSB_SEARCH, json=body, timeout=15)
        if r.status_code != 200:
            return {}
        hits = r.json().get("result_set", [])
        if not hits:
            return {}
        pdb_id = hits[0].get("identifier", "")
        if not pdb_id:
            return {}

        # Fetch resolution via GraphQL
        gql = requests.post(
            RCSB_GQL,
            json={
                "query": """
                    query($id: String!) {
                        entry(entry_id: $id) {
                            refine { ls_d_res_high }
                            exptl  { method }
                        }
                    }
                """,
                "variables": {"id": pdb_id},
            },
            timeout=10,
        )
        if gql.status_code != 200:
            return {"pdb_id": pdb_id}
        entry  = gql.json().get("data", {}).get("entry", {})
        refine = (entry.get("refine") or [{}])[0]
        exptl  = (entry.get("exptl")  or [{}])[0]
        return {
            "pdb_id":          pdb_id.upper(),
            "pdb_resolution_A": refine.get("ls_d_res_high"),
            "pdb_method":      exptl.get("method", "N/A"),
        }
    except Exception as exc:
        logger.debug(f"    RCSB search failed ({inchikey}): {exc}")
        return {}


def annotate_pdb_from_sifts(
    df: pd.DataFrame,
    sifts: dict,
    lookup_active_compound: bool = True,
) -> pd.DataFrame:
    """
    Add three PDB annotation columns using the local SIFTS mapping:

      n_target_pdb_structures  (int)
          Number of distinct PDB structures for the target protein's
          UniProt accession.  Reflects how well-structurally characterised
          the target is.  Derived entirely from the local SIFTS file —
          no network call.

      target_pdb_sample  (str)
          Comma-separated list of up to 5 representative PDB IDs for the
          target protein.  Useful as quick links for structural inspection.

      active_compound_pdb  (str)
          PDB ID(s) where the exact active compound (by InChIKey) is
          co-crystallized WITH the same target protein.  Uses RCSB
          combined chemical + UniProt entity search so there are no
          cross-target false positives.  Empty string if not found.
          Requires network access; skipped if `lookup_active_compound=False`.
      active_pdb_resolution_A  (float)
          Resolution of the best active-compound structure.
      active_pdb_method  (str)
          Experimental method (X-RAY DIFFRACTION, ELECTRON MICROSCOPY, ...).
    """
    df = df.copy()

    # ── Columns 1+2: target PDB count + sample from SIFTS (offline) ──────
    logger.info("  Annotating target PDB structures from SIFTS (offline)...")
    pdb_lists = df["uniprot_id"].apply(lambda u: _resolve_sifts_for_row(u, sifts))
    df["n_target_pdb_structures"] = pdb_lists.apply(len)
    df["target_pdb_sample"]       = pdb_lists.apply(
        lambda ids: ",".join(ids[:5]) if ids else "N/A"
    )
    n_with_pdb = (df["n_target_pdb_structures"] > 0).sum()
    logger.info(f"  Targets with >= 1 PDB structure: {n_with_pdb:,} / {len(df):,} pairs")

    # ── Column 3: active compound co-crystallised with target (RCSB) ──────
    df["active_compound_pdb"]   = ""
    df["active_pdb_resolution_A"] = float("nan")
    df["active_pdb_method"]     = "N/A"

    if not lookup_active_compound:
        logger.info("  Active-compound PDB lookup skipped (--no-compound-pdb).")
        return df

    # Only attempt for pairs where the target has at least one PDB structure
    has_pdb_mask = df["n_target_pdb_structures"] > 0
    if not has_pdb_mask.any():
        logger.info("  No target has a PDB structure — skipping compound PDB lookup.")
        return df

    # Deduplicate by (active_smiles, uniprot_id) to minimise RCSB calls
    from rdkit import Chem, RDLogger as _RL
    from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
    _RL.DisableLog("rdApp.*")

    unique_pairs = (
        df[has_pdb_mask][["active_smiles", "uniprot_id"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    logger.info(
        f"  Searching RCSB for active compound co-crystallised with target "
        f"({len(unique_pairs):,} unique active+target pairs)..."
    )

    # Build InChIKey cache
    ikey_cache: dict = {}
    for smi in unique_pairs["active_smiles"].unique():
        mol = desalt_mol(Chem.MolFromSmiles(str(smi))) if isinstance(smi, str) else None
        if mol:
            try:
                ikey_cache[smi] = InchiToInchiKey(MolToInchi(mol))
            except Exception:
                ikey_cache[smi] = None
        else:
            ikey_cache[smi] = None

    # Per-(smiles, uniprot) RCSB combined search
    compound_pdb_cache: dict = {}  # {(smi, uid): result_dict}

    for _, row in tqdm(unique_pairs.iterrows(), total=len(unique_pairs),
                       desc="RCSB compound-in-target", unit="pair"):
        smi = row["active_smiles"]
        uid = row["uniprot_id"]
        key = (smi, uid)
        ikey = ikey_cache.get(smi)
        if not ikey:
            compound_pdb_cache[key] = {}
            continue
        target_pdbs = _resolve_sifts_for_row(uid, sifts)
        result = find_compound_in_target_pdb(ikey, uid, target_pdbs)
        compound_pdb_cache[key] = result
        time.sleep(0.05)

    def _apply_compound_pdb(row):
        key = (row["active_smiles"], row["uniprot_id"])
        return compound_pdb_cache.get(key, {})

    results = df[has_pdb_mask].apply(_apply_compound_pdb, axis=1)
    df.loc[has_pdb_mask, "active_compound_pdb"]    = results.apply(
        lambda r: r.get("pdb_id", ""))
    df.loc[has_pdb_mask, "active_pdb_resolution_A"] = results.apply(
        lambda r: r.get("pdb_resolution_A", float("nan")))
    df.loc[has_pdb_mask, "active_pdb_method"]       = results.apply(
        lambda r: r.get("pdb_method", "N/A"))

    n_found = (df["active_compound_pdb"] != "").sum()
    logger.info(f"  Active compound found in PDB: {n_found:,} pairs")
    return df


def enrich(df: pd.DataFrame,
           sqlite_path: Optional[str],
           sifts_path: Optional[str] = None,
           maps_dir: Optional[str] = None,
           skip_protein_class: bool = False,
           use_uniprot_family: bool = True,
           lookup_active_compound_pdb: bool = True) -> pd.DataFrame:

    df = df.copy()
    n = len(df)
    logger.info(f"Enriching {n:,} cliff pairs...")

    # ── (a) Tanimoto ─────────────────────────────────────────────────────────
    logger.info("\n[Step 1] Tanimoto ECFP4 (RDKit, vectorised)")
    df["tanimoto_ecfp4"] = compute_tanimoto_column(df)

    # ── (d) Assay IDs ─────────────────────────────────────────────────────────
    logger.info("\n[Step 2] Assay ID columns")
    # Rename existing column → active_assay
    df = df.rename(columns={"assay_chembl_id": "active_assay_chembl_id"})

    # inactive_assay: same as active for XC50-XC50; SQLite lookup for pct_fallback
    df["inactive_assay_chembl_id"] = df["active_assay_chembl_id"]  # default

    fallback_mask = df.get("cliff_type", pd.Series(dtype=str)) == "pct_fallback"
    n_fallback = fallback_mask.sum()
    logger.info(f"  XC50-XC50 pairs : {(~fallback_mask).sum():,}  (inactive assay = active assay)")
    logger.info(f"  pct_fallback     : {n_fallback:,}  (inactive assay needs SQLite lookup)")

    if n_fallback > 0:
        if sqlite_path:
            fb_df = df[fallback_mask]
            inactive_ids = fb_df["inactive_chembl_id"].dropna().unique().tolist()
            target_ids_fb = fb_df["target_chembl_id"].dropna().unique().tolist()
            logger.info(f"  Querying inactive assay IDs for {len(inactive_ids):,} molecules...")
            ina_assay_map = fetch_inactive_assay_sqlite(sqlite_path, inactive_ids, target_ids_fb)
            logger.info(f"  Found {len(ina_assay_map):,} inactive assay mappings")

            def resolve_inactive_assay(row):
                key = (row["inactive_chembl_id"], row["target_chembl_id"])
                return ina_assay_map.get(key, "N/A")

            df.loc[fallback_mask, "inactive_assay_chembl_id"] = \
                df[fallback_mask].apply(resolve_inactive_assay, axis=1)
        else:
            logger.warning(
                "  pct_fallback pairs need --sqlite to resolve inactive assay IDs.\n"
                "  Set to 'N/A (needs SQLite)' for now."
            )
            df.loc[fallback_mask, "inactive_assay_chembl_id"] = "N/A (needs SQLite)"

    # ── (b) UniProt ID fill ───────────────────────────────────────────────────
    logger.info("\n[Step 3] UniProt ID verification / fill")
    if "uniprot_id" not in df.columns:
        df["uniprot_id"] = "N/A"

    missing_uid = df["uniprot_id"].isna() | (df["uniprot_id"] == "N/A")
    n_missing   = missing_uid.sum()
    logger.info(f"  UniProt ID present   : {(~missing_uid).sum():,}")
    logger.info(f"  UniProt ID missing   : {n_missing:,}")

    if n_missing > 0:
        gap_targets = df.loc[missing_uid, "target_chembl_id"].dropna().unique().tolist()
        if sqlite_path and gap_targets:
            uid_map = fill_uniprot_sqlite(sqlite_path, gap_targets)
            df.loc[missing_uid, "uniprot_id"] = \
                df.loc[missing_uid, "target_chembl_id"].map(uid_map).fillna("N/A")
            still_missing = (df["uniprot_id"] == "N/A").sum()
            logger.info(f"  After SQLite fill     : {still_missing:,} still missing")

    # ── (c) Protein class ─────────────────────────────────────────────────────
    logger.info("\n[Step 4] Protein class")
    df["protein_class_l1"] = "N/A"
    df["protein_class_l2"] = "N/A"
    df["protein_class_l3"] = "N/A"
    df["protein_family"]   = "N/A"

    if skip_protein_class:
        logger.info("  Skipped (--skip-protein-class)")
    else:
        # ── SIFTS offline (EC number + Pfam) — best source ───────────────
        sifts_classified = False
        if maps_dir and os.path.isdir(maps_dir):
            # Auto-discover files — accept both compressed and plain
            def _find(names):
                for n in names:
                    for suffix in ("", ".gz"):
                        p = os.path.join(maps_dir, n + suffix)
                        if os.path.exists(p):
                            return p
                return None

            f_updb   = _find(["uniprot_pdb.csv", "pdb_chain_uniprot.csv"])
            f_enzyme = _find(["pdb_chain_enzyme.csv"])
            f_pfam   = _find(["pdb_chain_pfam.csv"])

            if f_updb and f_enzyme and f_pfam:
                logger.info(f"  SIFTS offline classification from {maps_dir}")
                updb_map   = load_uniprot_pdb(f_updb)
                enzyme_map = load_pdb_enzyme(f_enzyme)
                pfam_map   = load_pdb_pfam(f_pfam)
                df = annotate_protein_class_sifts(df, updb_map, enzyme_map, pfam_map)
                sifts_classified = True
            else:
                missing = [n for n, f in [
                    ("uniprot_pdb.csv / pdb_chain_uniprot.csv", f_updb),
                    ("pdb_chain_enzyme.csv", f_enzyme),
                    ("pdb_chain_pfam.csv",   f_pfam),
                ] if not f]
                logger.warning(
                    f"  --maps-dir set but missing files: {missing}\n"
                    "  Falling back to ChEMBL REST."
                )

        unique_targets = df["target_chembl_id"].dropna().unique().tolist()
        pc_map: dict   = {}

        if not sifts_classified:
            # ── SQLite protein_class table ────────────────────────────────
            if sqlite_path:
                logger.info(f"  Querying protein_class from SQLite for {len(unique_targets):,} targets...")
                pc_map = fetch_protein_class_sqlite(sqlite_path, unique_targets)
                _useless = {"N/A", "SINGLE PROTEIN", "PROTEIN COMPLEX",
                            "PROTEIN FAMILY", "SELECTIVITY GROUP", "CHIMERIC PROTEIN", ""}
                n_from_sqlite = sum(
                    1 for v in pc_map.values()
                    if v.get("protein_class_l2", "N/A") not in _useless
                )
                logger.info(f"  Resolved {n_from_sqlite:,} targets with meaningful L2 from SQLite")

            # ── ChEMBL REST for remaining gaps ────────────────────────────
            _useless = {"N/A", "SINGLE PROTEIN", "PROTEIN COMPLEX",
                        "PROTEIN FAMILY", "SELECTIVITY GROUP", "CHIMERIC PROTEIN", ""}
            missing_l2 = [t for t in unique_targets
                          if pc_map.get(t, {}).get("protein_class_l2", "N/A") in _useless]
            if missing_l2:
                logger.info(f"  {len(missing_l2):,} targets need ChEMBL REST classification...")
                try:
                    rest_map = fetch_chembl_target_class_api(missing_l2)
                    pc_map.update({k: v for k, v in rest_map.items() if v})
                except Exception as exc:
                    logger.warning(f"  ChEMBL REST target class failed: {exc}")

            for field in ("protein_class_l1", "protein_class_l2",
                          "protein_class_l3", "protein_family"):
                df[field] = df["target_chembl_id"].map(
                    lambda t, f=field: pc_map.get(t, {}).get(f, "N/A")
                )

        n_l2 = (df["protein_class_l2"] != "N/A").sum()
        logger.info(f"  protein_class_l2 coverage: {n_l2:,} / {n:,} rows")

        # UniProt SIMILARITY text — only for targets where ChEMBL REST also returned nothing
        no_fam = df["protein_family"] == "N/A"
        if use_uniprot_family and no_fam.any():
            uids = df.loc[no_fam, "uniprot_id"].dropna().unique().tolist()
            uids = [u for u in uids if u and u != "N/A"]
            if uids:
                logger.info(f"  UniProt family fallback for {len(uids):,} accessions "
                            f"not resolved by ChEMBL REST...")
                uniprot_fam = fetch_protein_family_uniprot(uids)
                df.loc[no_fam, "protein_family"] = \
                    df.loc[no_fam, "uniprot_id"].map(uniprot_fam).fillna("N/A")
                n_filled = (df["protein_family"] != "N/A").sum()
                logger.info(f"  protein_family after UniProt fill: {n_filled:,} / {n:,}")

    # ── SIFTS: target PDB count + active compound co-crystal lookup ───────
    if sifts_path or maps_dir:
        logger.info("\n[SIFTS] PDB structure annotation...")
        # Resolve the chain-level or compact uniprot->pdb file
        _sifts_file = sifts_path
        if not _sifts_file and maps_dir:
            for name in ("uniprot_pdb.csv", "pdb_chain_uniprot.csv"):
                candidate = os.path.join(maps_dir, name)
                if os.path.exists(candidate):
                    _sifts_file = candidate
                    break

        if _sifts_file:
            # Determine if this is the compact uniprot_pdb or chain-level file
            with open(_sifts_file) as _f:
                for _l in _f:
                    if _l.startswith('#'):
                        continue
                    _hdr = _l.strip()
                    break
            _is_compact = 'PDB_IDS' in _hdr.upper() or _hdr.upper().startswith('SP_PRIMARY,PDB')
            if _is_compact:
                sifts = load_uniprot_pdb(_sifts_file)
            else:
                sifts = load_sifts(_sifts_file)

            df = annotate_pdb_from_sifts(
                df, sifts,
                lookup_active_compound=lookup_active_compound_pdb,
            )
        else:
            logger.warning("  No SIFTS UniProt->PDB file found — skipping PDB annotation.")
    else:
        logger.info("  SIFTS not provided — skipping PDB annotation "
                    "(use --sifts path/to/pdb_chain_uniprot.csv)")

    return df


# ═════════════════════════════════════════════════════════════════════════════
# OUTPUT COLUMN ORDER
# ═════════════════════════════════════════════════════════════════════════════

# All original columns are preserved; new columns are inserted at logical positions.
NEW_COLS_AFTER = {
    "inactive_std_type":          ["inactive_pct_activity"],
    "inactive_pct_activity":      ["tanimoto_ecfp4"],
    "delta_pXC50":                ["active_assay_chembl_id", "inactive_assay_chembl_id"],
    "assay_confidence_score":     ["protein_class_l1", "protein_class_l2",
                                   "protein_class_l3", "protein_family"],
    "protein_family":             ["n_target_pdb_structures", "target_pdb_sample",
                                   "active_compound_pdb", "active_pdb_resolution_A",
                                   "active_pdb_method"],
}


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Place new annotation columns immediately after their logical anchor columns.
    Works correctly even when new columns already exist in df (no duplicates).
    """
    # Build the desired order: start from original cols but move any
    # NEW_COLS_AFTER entries to their anchored positions.
    anchor_new: dict = {}
    for anchor, news in NEW_COLS_AFTER.items():
        anchor_new[anchor] = [c for c in news if c in df.columns]

    # Columns that will be repositioned (removed from natural position)
    repositioned = {c for news in anchor_new.values() for c in news}

    ordered = []
    seen = set()

    for c in df.columns:
        if c in repositioned:
            continue           # will be placed after its anchor
        if c in seen:
            continue           # deduplicate (shouldn't happen, but be safe)
        ordered.append(c)
        seen.add(c)
        # Insert repositioned columns right after their anchor
        for nc in anchor_new.get(c, []):
            if nc not in seen and nc in df.columns:
                ordered.append(nc)
                seen.add(nc)

    # Any remaining columns not yet placed (e.g. if anchor missing)
    for c in df.columns:
        if c not in seen:
            ordered.append(c)
            seen.add(c)

    # Drop accidental duplicate column names before selecting
    df = df.loc[:, ~df.columns.duplicated()]
    return df[[c for c in ordered if c in df.columns]]


# ═════════════════════════════════════════════════════════════════════════════
# CLI
# ═════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",  "-i", required=True,
                   help="Input all_cliffs.csv from activity_cliff_scanner.py")
    p.add_argument("--output", "-o", default=None,
                   help="Output CSV (default: <input>_enriched.csv)")
    p.add_argument("--sqlite", metavar="PATH", default=None,
                   help="ChEMBL SQLite file or directory.  Required for protein "
                        "class and inactive assay IDs.")
    p.add_argument("--skip-protein-class", action="store_true",
                   help="Skip protein class lookup (steps 4).  "
                        "Useful for a fast Tanimoto-only run.")
    p.add_argument("--skip-uniprot-family", action="store_true",
                   help="Skip UniProt SIMILARITY comment lookup "
                        "(used only as fallback when SQLite protein_class is missing).")
    p.add_argument("--chunksize", type=int, default=0,
                   help="Process input in chunks of N rows (0 = all at once).  "
                        "Use for very large files (>500 K rows) if RAM is tight.")
    p.add_argument("--maps-dir", metavar="DIR", default=None,
                   help=(
                       "Directory containing SIFTS flat files downloaded from "
                       "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/ .\n"
                       "Expected files (plain or .gz):\n"
                       "  uniprot_pdb.csv          (or pdb_chain_uniprot.csv)\n"
                       "  pdb_chain_enzyme.csv      -> protein_class via EC number\n"
                       "  pdb_chain_pfam.csv         -> protein_class via Pfam domain\n"
                       "When present, SIFTS-based classification replaces the ChEMBL "
                       "REST approach entirely — fully offline, >90%% coverage for "
                       "drug targets."
                   ))
    p.add_argument("--sifts", metavar="FILE", default=None,
                   help=(
                       "Path to pdb_chain_uniprot.csv (or uniprot_pdb.csv) for PDB "
                       "structure annotation only.  For full protein classification "
                       "use --maps-dir instead."
                   ))
    p.add_argument("--no-compound-pdb", action="store_true",
                   help=(
                       "With --sifts: annotate target PDB counts only; skip the "
                       "per-compound RCSB search for active_compound_pdb. "
                       "Much faster (no network calls)."
                   ))
    return p.parse_args()


def main():
    args = parse_args()
    output = args.output or args.input.replace(".csv", "_enriched.csv")

    logger.info("=" * 60)
    logger.info("  Activity Cliff Enricher")
    logger.info("=" * 60)
    logger.info(f"  Input    : {args.input}")
    logger.info(f"  Output   : {output}")
    logger.info(f"  SQLite   : {args.sqlite or 'not provided -- REST API fallback'}")

    # ── Load ─────────────────────────────────────────────────────────────────
    logger.info(f"\nLoading {args.input}...")
    df = pd.read_csv(args.input, low_memory=False)
    logger.info(f"  {len(df):,} rows, {len(df.columns)} columns")

    # ── Enrich ───────────────────────────────────────────────────────────────
    df_enriched = enrich(
        df,
        sqlite_path=args.sqlite,
        sifts_path=args.sifts,
        maps_dir=args.maps_dir,
        skip_protein_class=args.skip_protein_class,
        use_uniprot_family=not args.skip_uniprot_family,
        lookup_active_compound_pdb=(not args.no_compound_pdb),
    )

    # ── Reorder columns ───────────────────────────────────────────────────────
    df_enriched = reorder_columns(df_enriched)

    # ── Export ────────────────────────────────────────────────────────────────
    logger.info(f"\nWriting {output}...")
    df_enriched.to_csv(output, index=False, float_format="%.4f")

    # ── Summary ───────────────────────────────────────────────────────────────
    logger.info("\n" + "=" * 60)
    logger.info(f"  Rows written         : {len(df_enriched):,}")
    logger.info(f"  Columns              : {len(df_enriched.columns)}")
    logger.info(f"  New columns added    : tanimoto_ecfp4, active_assay_chembl_id,")
    logger.info(f"                         inactive_assay_chembl_id, protein_class_l1,")
    logger.info(f"                         protein_class_l2, protein_class_l3, protein_family")
    if "tanimoto_ecfp4" in df_enriched.columns:
        t = df_enriched["tanimoto_ecfp4"].dropna().to_numpy()
        if len(t):
            import numpy as _np
            logger.info(f"  Tanimoto mean/min    : {_np.mean(t):.4f} / {_np.min(t):.4f}")
    if "protein_class_l1" in df_enriched.columns:
        col = df_enriched["protein_class_l1"].astype(str).fillna("N/A")
        top = col.value_counts().head(5)
        logger.info("  Top protein classes  :")
        for cls, cnt in top.items():
            logger.info(f"    {str(cls):<30} {int(cnt):>8,}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
