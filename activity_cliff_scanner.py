#!/usr/bin/env python3
"""
=============================================================================
Activity Cliff Scanner for ChEMBL  (v3 ? SQLite + REST API)
=============================================================================
Scans ChEMBL assays (confidence_score >= 8) for activity cliffs defined as:

  [Criterion 1]  Tanimoto(ECFP4, r=2, 2048 bits) >= TANIMOTO_CUTOFF (0.95)
  [Criterion 2]  |pXC50_active - pXC50_inactive|  >= DELTA_PXIC50   (2.0)
  [Criterion 3]  No structural AND activity intermediate exists in the assay
                 (disable with --no-intermediate-check for speed)

Fallback (no XC50 for one compound):
  Active   : pXC50 >= 5.0  (IC50 <= 10 ?M)
  Inactive : % inhibition < 50 % at ~10 ?M single-point screen

DATA SOURCES (in order of preference)
--------------------------------------
  1. --sqlite  PATH   Local ChEMBL SQLite file or directory ? FASTEST,
                      no network, full dataset.  Recommended.
  2. --input-csv CSV  Pre-exported CSV (e.g. from sqlite3 CLI or pandas).
  3. ChEMBL REST API  Automatic, paginated.  Needs www.ebi.ac.uk reachable.
  4. --demo           Built-in toy data.  Validates pipeline, no network.

USAGE EXAMPLES
--------------
  # Recommended: full ChEMBL 37, all confident assays, EGFR only
  python activity_cliff_scanner.py \
      --sqlite "C:/path/to/chembl_37_sqlite" \
      --target CHEMBL279 \
      --output EGFR_cliffs.csv

  # All targets in the SQLite (large ? use --max-records to cap)
  python activity_cliff_scanner.py \
      --sqlite "C:/path/to/chembl_37_sqlite" \
      --max-records 500000 \
      --output all_cliffs.csv

  # From a pre-exported CSV
  python activity_cliff_scanner.py \
      --input-csv my_activities.csv \
      --output cliffs.csv

  # REST API (single target; network required)
  python activity_cliff_scanner.py \
      --target CHEMBL279 \
      --output EGFR_cliffs.csv

  # Demo (no network, no SQLite)
  python activity_cliff_scanner.py --demo

REQUIRED PACKAGES
-----------------
  pip install rdkit pandas numpy requests tqdm
  (chembl-webresource-client is NO LONGER required)

OUTPUT CSV COLUMNS
------------------
  active_chembl_id, active_smiles, active_pXC50, active_value_nM,
  active_std_type,
  inactive_chembl_id, inactive_smiles, inactive_pXC50, inactive_value_nM,
  inactive_std_type, [inactive_pct_activity ? fallback only],
  tanimoto_ecfp4, delta_pXC50,
  assay_chembl_id, target_chembl_id, target_name, target_organism,
  uniprot_id, assay_confidence_score,
  active_n_tautomers, active_tautomer_important,
  inactive_n_tautomers, inactive_tautomer_important,
  active_charge_ph7, inactive_charge_ph7,
  cliff_type
=============================================================================
"""

import os
import sys
import glob
import time
import sqlite3
import logging
import argparse
import warnings
from itertools import combinations
from typing import Optional

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

# ?? RDKit ????????????????????????????????????????????????????????????????????
from rdkit import Chem, RDLogger
from rdkit.Chem import DataStructs, rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

RDLogger.DisableLog("rdApp.*")          # silence ALL RDKit warnings/info/errors
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ?????????????????????????????????????????????????????????????????????????????
# CONSTANTS
# ?????????????????????????????????????????????????????????????????????????????

CONFIDENCE_CUTOFF   = 8
TANIMOTO_CUTOFF     = 0.95
DELTA_PXIC50_CUTOFF = 2.0
ECFP_RADIUS         = 2
ECFP_NBITS          = 2048
PXIC50_ACTIVE_MIN   = 5.0
PCT_INACTIVE_THRESH = 50.0
PCT_CONC_UM_LO      = 3.0
PCT_CONC_UM_HI      = 30.0
INTER_SIM_THRESH    = 0.85
API_DELAY_S         = 0.06
MAX_TAUTOMERS       = 64

# "XC50" is used throughout the code as a conceptual umbrella for any
# dose-response or binding affinity measurement that yields a molar potency.
# None of these are a literal ChEMBL standard_type called "XC50".
# Units that cannot be converted to nM (e.g. MIC in ?g/mL) are silently
# dropped during preprocessing by the to_nM() NaN filter.
XC50_TYPES = (
    "IC50",     # inhibitory concentration 50 %          ? most common
    "EC50",     # effective concentration 50 %
    "DC50",     # degradation concentration 50 %         ? PROTACs / degraders
    "AC50",     # activity concentration 50 %            ? functional assays
    "CC50",     # cytotoxic concentration 50 %           ? cell viability
    "GI50",     # growth inhibition 50 %                 ? antiproliferative
    "Ki",       # equilibrium inhibition constant
    "Kd",       # equilibrium dissociation constant
    "Kb",       # binding constant (inverse convention of Kd, some SPR papers)
    "ED50",     # effective dose 50 %
    "Potency",  # generic label used in some HTS deposits
)

PCT_TYPES = ("Inhibition", "Activity")

CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"

# ?? ChEMBL SQLite activity quality filters ???????????????????????????????????
# These mirror what the ChEMBL team recommends for curated, publication-quality
# binding data.  Both conditions must hold:
#   potential_duplicate = 0  ? removes flagged duplicates
#   data_validity_comment IS NULL  ? removes manually flagged outliers
#   standard_relation = '='  ? excludes ">" / "<" threshold-only values
SQLITE_QUALITY_FILTER = """
    AND (act.potential_duplicate = 0 OR act.potential_duplicate IS NULL)
    AND (act.data_validity_comment IS NULL
         OR act.data_validity_comment = 'Manually validated')
"""


# ?????????????????????????????????????????????????????????????????????????????
# IONIZABLE GROUP SMARTS  (pH 7 charge annotation)
# ?????????????????????????????????????????????????????????????????????????????

IONIZABLE_GROUPS = {
    # name: (smarts, dominant_form_pH7, typical_pKa)
    "primary_amine":   ("[NX3H2;!$(NC=O);!$(Nc)]",                    "+",        "9?11"),
    "secondary_amine": ("[NX3H1;!$(NC=O);!$([NH1]c)]",               "+",        "9?11"),
    "tertiary_amine":  ("[NX3H0;$(N([CX4])[CX4][CX4]);!$(NC=O)]",    "+",        "8?10"),
    "guanidine":       ("NC(=N)N",                                     "+",        "12?13"),
    "amidine":         ("[CX3](=N)N",                                  "+",        "10?12"),
    "pyridine_N":      ("[n;H0;+0]",                                   "partial+", "~5"),
    "imidazole_N":     ("[nH]1ccnc1",                                  "partial+", "~6"),
    "piperazine_N":    ("N1CCNCC1",                                    "+",        "5?9"),
    "carboxylic_acid": ("[CX3](=O)[OX2H1]",                           "-",        "3?5"),
    "tetrazole":       ("[nH]1nnnc1",                                  "-",        "~5"),
    "sulfonamide_NH":  ("[NH1]S(=O)(=O)",                             "partial-", "~10"),
    "phosphate":       ("[PX4](=O)([OX2H])[OX2H]",                   "-/-",      "1/6"),
    "sulfonic_acid":   ("S(=O)(=O)[OH]",                              "-",        "<2"),
    "phenol":          ("[c][OX2H]",                                   "neutral",  "~10"),
    "thiol":           ("[SX2H]",                                      "neutral",  "~8"),
}


# ?????????????????????????????????????????????????????????????????????????????
# UNIT CONVERSION
# ?????????????????????????????????????????????????????????????????????????????

_UNIT_TO_NM = {
    "nm": 1.0, "pm": 1e-3,
    "um": 1e3, "?m": 1e3, "?m": 1e3,
    "mm": 1e6, "m":  1e9,
    "mol/l": 1e9, "nmol/l": 1.0, "umol/l": 1e3, "mmol/l": 1e6,
}

def to_nM(value: float, unit: str) -> float:
    return value * _UNIT_TO_NM.get(str(unit).strip().lower(), float("nan"))

def to_pXC50(value_nM: float) -> float:
    if not np.isfinite(value_nM) or value_nM <= 0:
        return float("nan")
    return -np.log10(value_nM * 1e-9)


# ?????????????????????????????????????????????????????????????????????????????
# MOLECULAR CHEMISTRY UTILITIES
# ?????????????????????????????????????????????????????????????????????????????

_MG  = rdFingerprintGenerator.GetMorganGenerator(radius=ECFP_RADIUS, fpSize=ECFP_NBITS)
_TE  = rdMolStandardize.TautomerEnumerator()
_TE.SetMaxTautomers(MAX_TAUTOMERS)
_LFC = rdMolStandardize.LargestFragmentChooser()
_UN  = rdMolStandardize.Uncharger()          # for salt-form normalisation


def smi_to_mol(smiles: str) -> Optional[Chem.Mol]:
    if not smiles or not isinstance(smiles, str):
        return None
    return Chem.MolFromSmiles(smiles.strip())

def mol_to_fp(mol: Optional[Chem.Mol]):
    return _MG.GetFingerprint(mol) if mol is not None else None

def tanimoto(fp1, fp2) -> float:
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def bulk_tanimoto(fp_query, fp_list) -> np.ndarray:
    if fp_query is None:
        return np.zeros(len(fp_list))
    valid = [(i, fp) for i, fp in enumerate(fp_list) if fp is not None]
    out = np.zeros(len(fp_list))
    if valid:
        idxs, fps = zip(*valid)
        sims = DataStructs.BulkTanimotoSimilarity(fp_query, list(fps))
        for i, s in zip(idxs, sims):
            out[i] = s
    return out

def count_tautomers(mol: Optional[Chem.Mol]) -> int:
    if mol is None:
        return 1
    try:
        return len(list(_TE.Enumerate(mol)))
    except Exception:
        return 1

def inchikey_from_mol(mol: Optional[Chem.Mol]) -> Optional[str]:
    if mol is None:
        return None
    try:
        inchi = MolToInchi(mol)
        return InchiToInchiKey(inchi) if inchi else None
    except Exception:
        return None

def parent_inchikey_from_mol(mol: Optional[Chem.Mol]) -> Optional[str]:
    """
    InChIKey of the *parent* structure after:
      1. Largest-fragment selection  ? strips counterions / co-crystals
      2. Charge neutralisation       ? collapses free-base / salt / zwitterion
                                       to the same neutral parent

    Used exclusively for identity filtering; the original registered SMILES
    is always preserved in the output.
    """
    if mol is None:
        return None
    try:
        parent = _LFC.choose(mol)   # e.g. compound?HCl ? compound
        parent = _UN.uncharge(parent)  # e.g. ?COO? ? ?COOH, ?NH?? ? ?NH?
        return inchikey_from_mol(parent)
    except Exception:
        return inchikey_from_mol(mol)   # fall back to raw InChIKey

def get_charge_ph7(mol: Optional[Chem.Mol]) -> str:
    if mol is None:
        return "N/A"
    hits = []
    for name, (smarts, form, pka) in IONIZABLE_GROUPS.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            continue
        matches = mol.GetSubstructMatches(pat)
        if matches:
            hits.append(f"{name}(n={len(matches)}, pKa?{pka}?{form})")
    return "; ".join(hits) if hits else "neutral"


# ?????????????????????????????????????????????????????????????????????????????
# SQLITE DATA FETCHING  ? recommended for ChEMBL 37 local database
# ?????????????????????????????????????????????????????????????????????????????

def find_sqlite_db(path: str) -> str:
    """
    Accept either:
      - Path to the .db file directly, or
      - Path to a directory containing a single .db file.
    Returns the resolved .db file path.
    """
    path = os.path.normpath(path)

    if os.path.isfile(path) and path.endswith(".db"):
        return path

    if os.path.isdir(path):
        # Look for .db files inside (also one level deep for extracted archives)
        candidates = (
            glob.glob(os.path.join(path, "*.db")) +
            glob.glob(os.path.join(path, "**", "*.db"), recursive=False)
        )
        candidates = [c for c in candidates if os.path.isfile(c)]
        if len(candidates) == 1:
            logger.info(f"  Found SQLite database: {candidates[0]}")
            return candidates[0]
        if len(candidates) > 1:
            raise FileNotFoundError(
                f"Multiple .db files found in {path}:\n"
                + "\n".join(f"  {c}" for c in candidates)
                + "\nPlease specify the exact file with --sqlite path/to/chembl_37.db"
            )
        raise FileNotFoundError(
            f"No .db file found in directory: {path}\n"
            "Expected something like chembl_37.db"
        )

    raise FileNotFoundError(f"Path does not exist or is not a .db file: {path}")


def _sqlite_conn(db_path: str) -> sqlite3.Connection:
    """Open SQLite connection with performance pragmas."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    # WAL mode + large cache for read-heavy scans
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA cache_size   = -131072")   # 128 MB cache
    conn.execute("PRAGMA temp_store   = MEMORY")
    conn.execute("PRAGMA mmap_size    = 2147483648") # 2 GB mmap
    conn.execute("PRAGMA synchronous  = OFF")        # read-only, safe
    return conn


def _target_clause(target_id: Optional[str]) -> tuple:
    """Return (SQL snippet, params) for optional target filter."""
    if target_id:
        return "AND td.chembl_id = ?", (target_id,)
    return "", ()


# ?? XC50 activities from SQLite ???????????????????????????????????????????????

_SQL_XC50 = """
SELECT
    md.chembl_id                        AS molecule_chembl_id,
    cs.canonical_smiles,
    CAST(act.standard_value  AS REAL)   AS standard_value,
    act.standard_units,
    act.standard_type,
    CAST(act.pchembl_value   AS REAL)   AS pchembl_value,
    ass.chembl_id                       AS assay_chembl_id,
    ass.assay_type,
    td.chembl_id                        AS target_chembl_id,
    td.pref_name                        AS target_name,
    td.organism                         AS target_organism,
    CAST(ass.confidence_score AS INTEGER) AS assay_confidence_score
FROM       activities          act
JOIN       assays              ass  ON act.assay_id  = ass.assay_id
JOIN       target_dictionary   td   ON ass.tid       = td.tid
JOIN       molecule_dictionary md   ON act.molregno  = md.molregno
JOIN       compound_structures cs   ON act.molregno  = cs.molregno
WHERE  ass.confidence_score >= {conf}
  AND  act.standard_relation = '='
  AND  act.standard_type     IN ({placeholders})
  AND  act.standard_value    IS NOT NULL
  AND  cs.canonical_smiles   IS NOT NULL
  AND  act.standard_units    IN ('nM','uM','?M','?M','pM','nM','nm','um',
                                  'nmol/l','umol/l','mmol/l','mol/l','mM','M')
  {quality}
  {target_clause}
"""


def fetch_xc50_from_sqlite(
    db_path:   str,
    target_id: Optional[str] = None,
    max_records: int = 0,          # 0 = no limit
) -> pd.DataFrame:
    """
    Pull all XC50-type binding measurements (IC50, EC50, Ki, Kd, AC50?)
    from a local ChEMBL SQLite database.

    Parameters
    ----------
    db_path     : path to chembl_XX.db
    target_id   : ChEMBL target ID to restrict (e.g. 'CHEMBL279'), or None
    max_records : cap on total rows returned (0 = all)
    """
    db_path = find_sqlite_db(db_path)
    logger.info(f"Querying XC50 activities from: {db_path}")

    ph = ", ".join(["?"] * len(XC50_TYPES))
    target_sql, target_params = _target_clause(target_id)
    limit_clause = f"LIMIT {max_records}" if max_records > 0 else ""

    sql = (
        _SQL_XC50.format(
            conf=CONFIDENCE_CUTOFF,
            placeholders=ph,
            quality=SQLITE_QUALITY_FILTER,
            target_clause=target_sql,
        )
        + limit_clause
    )
    params = list(XC50_TYPES) + list(target_params)

    conn = _sqlite_conn(db_path)
    try:
        logger.info("  Running XC50 SQL query (may take 15?60 s on first run)?")
        t0 = time.time()
        df = pd.read_sql_query(sql, conn, params=params)
        logger.info(f"  Retrieved {len(df):,} XC50 rows in {time.time()-t0:.1f} s")
    finally:
        conn.close()

    return df


# ?? % activity / inhibition from SQLite ??????????????????????????????????????

_SQL_PCT = """
SELECT
    md.chembl_id                        AS molecule_chembl_id,
    cs.canonical_smiles,
    CAST(act.standard_value  AS REAL)   AS standard_value,
    act.standard_units,
    act.standard_type,
    ass.chembl_id                       AS assay_chembl_id,
    ass.assay_type,
    td.chembl_id                        AS target_chembl_id,
    CAST(ass.confidence_score AS INTEGER) AS assay_confidence_score
FROM       activities          act
JOIN       assays              ass  ON act.assay_id  = ass.assay_id
JOIN       target_dictionary   td   ON ass.tid       = td.tid
JOIN       molecule_dictionary md   ON act.molregno  = md.molregno
JOIN       compound_structures cs   ON act.molregno  = cs.molregno
WHERE  ass.confidence_score >= {conf}
  AND  act.standard_relation  = '='
  AND  act.standard_type      IN ({placeholders})
  AND  act.standard_units     = '%'
  AND  act.standard_value     IS NOT NULL
  AND  cs.canonical_smiles    IS NOT NULL
  {quality}
  {target_clause}
"""


def fetch_pct_from_sqlite(
    db_path:   str,
    target_id: Optional[str] = None,
    max_records: int = 0,
) -> pd.DataFrame:
    db_path = find_sqlite_db(db_path)
    ph = ", ".join(["?"] * len(PCT_TYPES))
    target_sql, target_params = _target_clause(target_id)
    limit_clause = f"LIMIT {max_records}" if max_records > 0 else ""

    sql = (
        _SQL_PCT.format(
            conf=CONFIDENCE_CUTOFF,
            placeholders=ph,
            quality=SQLITE_QUALITY_FILTER,
            target_clause=target_sql,
        )
        + limit_clause
    )
    params = list(PCT_TYPES) + list(target_params)

    conn = _sqlite_conn(db_path)
    try:
        logger.info("  Running % activity SQL query?")
        df = pd.read_sql_query(sql, conn, params=params)
        logger.info(f"  Retrieved {len(df):,} % activity rows")
    finally:
        conn.close()

    return df


# ?? Target metadata from SQLite ???????????????????????????????????????????????

_SQL_TARGET_META = """
SELECT
    td.chembl_id    AS target_chembl_id,
    td.pref_name    AS target_name,
    td.organism,
    td.target_type,
    GROUP_CONCAT(DISTINCT cseq.accession) AS uniprot_ids
FROM target_dictionary  td
LEFT JOIN target_components  tc   ON td.tid          = tc.tid
LEFT JOIN component_sequences cseq ON tc.component_id = cseq.component_id
WHERE td.chembl_id IN ({placeholders})
GROUP BY td.chembl_id, td.pref_name, td.organism, td.target_type
"""

def fetch_target_meta_from_sqlite(db_path: str, target_ids: list) -> dict:
    if not target_ids:
        return {}
    db_path = find_sqlite_db(db_path)
    ph  = ", ".join(["?"] * len(target_ids))
    sql = _SQL_TARGET_META.format(placeholders=ph)

    conn = _sqlite_conn(db_path)
    try:
        rows = conn.execute(sql, target_ids).fetchall()
    finally:
        conn.close()

    meta = {}
    for row in rows:
        tid = row["target_chembl_id"]
        meta[tid] = {
            "uniprot_id":  row["uniprot_ids"] or "N/A",
            "target_name": row["target_name"] or "N/A",
            "organism":    row["organism"]    or "N/A",
            "target_type": row["target_type"] or "N/A",
        }
    return meta


# ?????????????????????????????????????????????????????????????????????????????
# CHEMBL REST API (fallback when no SQLite / no CSV)
# ?????????????????????????????????????????????????????????????????????????????

_SESSION = requests.Session()
_SESSION.headers.update({
    "Accept":     "application/json",
    "User-Agent": "activity-cliff-scanner/3.0",
})

_ACTIVITY_FIELDS = ",".join([
    "molecule_chembl_id", "canonical_smiles",
    "standard_value", "standard_units", "standard_type",
    "assay_chembl_id", "assay_type",
    "target_chembl_id", "target_organism",
    "assay_confidence_score", "pchembl_value",
])

def _chembl_paginate(endpoint, params, record_key, max_records=10_000, page_size=1000):
    params = {**params, "limit": page_size, "offset": 0, "format": "json"}
    records = []
    while True:
        try:
            r = _SESSION.get(f"{CHEMBL_API_BASE}/{endpoint}.json", params=params, timeout=30)
            r.raise_for_status()
            data = r.json()
        except Exception as exc:
            logger.warning(f"  REST error at offset {params['offset']}: {exc}")
            break
        batch = data.get(record_key, [])
        if not batch:
            break
        records.extend(batch)
        if len(records) >= max_records or not data.get("page_meta", {}).get("next"):
            break
        params["offset"] += page_size
        time.sleep(API_DELAY_S)
    return records[:max_records]

def fetch_xc50_activities(target_id=None, max_per_type=10_000):
    base = {
        "assay_confidence_score__gte": CONFIDENCE_CUTOFF,
        "standard_relation": "=",
        "standard_value__isnull": "false",
        "fields": _ACTIVITY_FIELDS,
    }
    if target_id:
        base["target_chembl_id"] = target_id
    records = []
    for stype in XC50_TYPES:
        logger.info(f"  REST: fetching {stype}?")
        batch = _chembl_paginate("activity", {**base, "standard_type": stype},
                                 "activities", max_per_type)
        logger.info(f"    {len(batch):>6,} {stype} records")
        records.extend(batch)
    return pd.DataFrame(records)

def fetch_pct_activities(target_id=None, max_per_type=5_000):
    base = {
        "assay_confidence_score__gte": CONFIDENCE_CUTOFF,
        "standard_relation": "=",
        "standard_value__isnull": "false",
        "standard_units": "%",
        "fields": _ACTIVITY_FIELDS,
    }
    if target_id:
        base["target_chembl_id"] = target_id
    records = []
    for stype in PCT_TYPES:
        logger.info(f"  REST: fetching {stype}%?")
        batch = _chembl_paginate("activity", {**base, "standard_type": stype},
                                 "activities", max_per_type)
        logger.info(f"    {len(batch):>6,} records")
        records.extend(batch)
    return pd.DataFrame(records)

def fetch_target_metadata(target_ids: list) -> dict:
    meta = {}
    for tid in tqdm(target_ids, desc="Target metadata (REST)", unit="tgt"):
        try:
            r = _SESSION.get(f"{CHEMBL_API_BASE}/target/{tid}.json", timeout=15)
            r.raise_for_status()
            data = r.json()
            uids = [
                xref.get("xref_id", "")
                for comp in data.get("target_components", [])
                for xref in comp.get("target_component_xrefs", [])
                if xref.get("xref_src_db") == "UniProt"
            ]
            meta[tid] = {
                "uniprot_id":  "|".join(uids) or "N/A",
                "target_name": data.get("pref_name",   "N/A"),
                "organism":    data.get("organism",     "N/A"),
            }
            time.sleep(API_DELAY_S)
        except Exception as exc:
            logger.debug(f"Target {tid}: {exc}")
            meta[tid] = {}
    return meta


# ?????????????????????????????????????????????????????????????????????????????



# ?????????????????????????????????????????????????????????????????????????????
# DATA PREPROCESSING
# ?????????????????????????????????????????????????????????????????????????????


# ═════════════════════════════════════════════════════════════════════════════
# EXTERNAL CSV LOADER  (non-ChEMBL data with SMILES + pChEMBL value)
# ═════════════════════════════════════════════════════════════════════════════

# Column name candidates tried in priority order for auto-detection
_SMILES_CANDIDATES   = ["smiles", "canonical_smiles", "structure", "mol",
                         "molecule", "smi", "Smiles", "SMILES"]
_ACTIVITY_CANDIDATES = ["pchembl_value", "pxc50", "pic50", "pki", "pkd",
                         "pec50", "pac50", "pactivity", "p_activity",
                         "activity", "potency", "value", "p_value"]
_ID_CANDIDATES       = ["id", "compound_id", "molecule_id", "chembl_id",
                         "name", "compound_name", "mol_id", "cmpd_id",
                         "molregno", "reg_id"]
_ASSAY_CANDIDATES    = ["assay_id", "assay", "assay_chembl_id", "group",
                         "series", "project", "batch", "dataset"]
_TARGET_CANDIDATES   = ["target_id", "target_chembl_id", "target", "protein",
                         "gene", "receptor"]


def _detect_col(df: pd.DataFrame, candidates: list, explicit: Optional[str]) -> Optional[str]:
    """Return `explicit` if given; otherwise the first candidate found in df.columns."""
    if explicit:
        if explicit in df.columns:
            return explicit
        raise ValueError(
            f"Column '{explicit}' not found. Available: {list(df.columns)}"
        )
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def load_external_csv(
    path: str,
    smiles_col:   Optional[str] = None,
    activity_col: Optional[str] = None,
    id_col:       Optional[str] = None,
    assay_col:    Optional[str] = None,
    target_col:   Optional[str] = None,
    default_assay:  str = "EXTERNAL_ASSAY",
    default_target: str = "EXTERNAL_TARGET",
    default_type:   str = "IC50",
) -> pd.DataFrame:
    """
    Load an external CSV containing SMILES and pChEMBL values and translate it
    into the internal dataframe format expected by preprocess_xc50().

    Column detection
    ----------------
    Column names are auto-detected from common synonyms.  Use the explicit
    override flags (--smiles-col, --activity-col, etc.) when auto-detection
    picks the wrong column.

    Required
    --------
    smiles_col    SMILES string for each compound.
    activity_col  pChEMBL value = -log10(IC50/EC50/Ki/Kd in M).
                  Values must be in pChEMBL units (typically 5-12).
                  NOT nanomolar concentrations — those go through --input-csv.

    Optional
    --------
    id_col        Compound identifier.  Auto-generated (EXT_00001 ...) if absent.
    assay_col     Grouping column; compounds in different groups are not paired.
                  If absent, all compounds are treated as one virtual assay.
    target_col    Target identifier (propagated to output for filtering).

    Internal mapping
    ----------------
    The function back-converts pChEMBL to nM for standard_value so the existing
    preprocess_xc50() pipeline works unchanged:
        standard_value (nM) = 10^(9 - pchembl_value)
    """
    logger.info(f"Loading external CSV: {path}")
    df = pd.read_csv(path, low_memory=False)
    logger.info(f"  {len(df):,} rows, {len(df.columns)} columns: {list(df.columns)}")

    # ── Detect columns ────────────────────────────────────────────────────
    smi_c  = _detect_col(df, _SMILES_CANDIDATES,   smiles_col)
    act_c  = _detect_col(df, _ACTIVITY_CANDIDATES, activity_col)
    id_c   = _detect_col(df, _ID_CANDIDATES,       id_col)
    assay_c= _detect_col(df, _ASSAY_CANDIDATES,    assay_col)
    tgt_c  = _detect_col(df, _TARGET_CANDIDATES,   target_col)

    if smi_c is None:
        raise ValueError(
            f"Cannot detect a SMILES column in {list(df.columns)}.\n"
            f"Use --smiles-col to specify it explicitly."
        )
    if act_c is None:
        raise ValueError(
            f"Cannot detect a pChEMBL activity column in {list(df.columns)}.\n"
            f"Use --activity-col to specify it.\n"
            f"The column must contain pChEMBL values (-log10 molar), "
            f"NOT raw IC50 concentrations."
        )

    _assay_label  = assay_c  or f"(single virtual assay: {default_assay})"
    _target_label = tgt_c    or f"(default: {default_target})"
    logger.info(f"  SMILES      : '{smi_c}'")
    logger.info(f"  pChEMBL     : '{act_c}'")
    logger.info(f"  Compound ID : '{id_c or '(auto-generated)'}'")
    logger.info(f"  Assay group : '{_assay_label}'")
    logger.info(f"  Target      : '{_target_label}'")
    # ── Build internal dataframe ──────────────────────────────────────────
    out = pd.DataFrame()
    out["canonical_smiles"] = df[smi_c].astype(str).str.strip()

    # pChEMBL value
    out["pchembl_value"] = pd.to_numeric(df[act_c], errors="coerce")
    invalid = out["pchembl_value"].isna().sum()
    if invalid:
        logger.warning(f"  {invalid:,} rows with non-numeric activity dropped")

    # Back-convert pChEMBL → nM so preprocess_xc50() standard_value path works
    # pChEMBL = -log10(M)  →  M = 10^(-pChEMBL)  →  nM = 10^(9-pChEMBL)
    out["standard_value"] = 10.0 ** (9.0 - out["pchembl_value"])
    out["standard_units"] = "nM"
    out["standard_type"]  = default_type

    # Compound ID
    if id_c:
        out["molecule_chembl_id"] = df[id_c].astype(str).str.strip()
    else:
        width = len(str(len(df)))
        out["molecule_chembl_id"] = [f"EXT_{i:0{width}d}" for i in range(len(df))]

    # ── Assay group ───────────────────────────────────────────────────────
    # Priority:
    #   1. Explicit --assay-col           → group by that column
    #   2. No assay col, but target col   → group by target (one "assay" per target)
    #   3. Neither                        → single virtual assay (all vs all)
    if assay_c:
        out["assay_chembl_id"] = df[assay_c].astype(str).str.strip()
        _grp_src = f"assay column '{assay_c}'"
    elif tgt_c:
        out["assay_chembl_id"] = df[tgt_c].astype(str).str.strip()
        _grp_src = (f"target column '{tgt_c}' (no --assay-col given; "
                    f"compounds grouped per target — cross-target pairs excluded)")
    else:
        out["assay_chembl_id"] = default_assay
        _grp_src = f"single virtual assay '{default_assay}' (no grouping column)"

    logger.info(f"  Grouping by  : {_grp_src}")

    # Target
    out["target_chembl_id"] = (
        df[tgt_c].astype(str).str.strip() if tgt_c else default_target
    )
    out["target_name"]     = out["target_chembl_id"]
    out["target_organism"] = "N/A"

    # Confidence score — external data treated as fully trusted
    out["assay_confidence_score"] = 9

    # Preserve any extra columns the user may have (for reference in output)
    extra = [c for c in df.columns
             if c not in (smi_c, act_c, id_c, assay_c, tgt_c)]
    for c in extra:
        if c not in out.columns:
            out[f"ext_{c}"] = df[c].values

    logger.info(
        f"  External CSV loaded: {len(out):,} rows, "
        f"{out['assay_chembl_id'].nunique():,} assay group(s), "
        f"{out['target_chembl_id'].nunique():,} target(s)"
    )
    return out


def preprocess_xc50(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value", "canonical_smiles"])
    df = df[df["standard_value"] > 0]

    df["value_nM"] = df.apply(
        lambda r: to_nM(r["standard_value"], r.get("standard_units", "nM")), axis=1
    )
    df = df[df["value_nM"] > 0].dropna(subset=["value_nM"])

    # Prefer ChEMBL's pre-computed pchembl_value
    if "pchembl_value" in df.columns:
        df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
        df["pXC50"] = df["pchembl_value"].where(
            df["pchembl_value"].notna(), df["value_nM"].map(to_pXC50)
        )
    else:
        df["pXC50"] = df["value_nM"].map(to_pXC50)
    df = df.dropna(subset=["pXC50"])

    # Chemistry objects
    df["mol"]           = df["canonical_smiles"].map(smi_to_mol)
    df = df[df["mol"].notna()]
    df["fp"]            = df["mol"].map(mol_to_fp)
    df["inchikey"]      = df["mol"].map(inchikey_from_mol)
    df["parent_inchikey"] = df["mol"].map(parent_inchikey_from_mol)
    df = df[df["fp"].notna()]

    # Keep most potent per molecule?assay pair.
    # Three dedup passes in order of specificity:
    #   1. Same ChEMBL ID + assay          (exact duplicate rows)
    #   2. Same raw InChIKey + assay        (different ID, same structure)
    #   3. Same parent InChIKey + assay     (salt / free-base / zwitterion forms)
    df = df.sort_values("pXC50", ascending=False)
    df = df.drop_duplicates(subset=["molecule_chembl_id", "assay_chembl_id"])
    df = df.drop_duplicates(subset=["inchikey",           "assay_chembl_id"])
    df = df.drop_duplicates(subset=["parent_inchikey",    "assay_chembl_id"])
    logger.info(f"  Preprocessed XC50 records: {len(df):,}")
    return df


def preprocess_pct(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value", "canonical_smiles"])

    if "assay_concentration_value" in df.columns:
        df["conc_val"] = pd.to_numeric(df["assay_concentration_value"], errors="coerce")
        df = df[df["conc_val"].between(PCT_CONC_UM_LO, PCT_CONC_UM_HI) | df["conc_val"].isna()]

    df["mol"] = df["canonical_smiles"].map(smi_to_mol)
    df = df[df["mol"].notna()]
    df["fp"]  = df["mol"].map(mol_to_fp)
    df = df[df["fp"].notna()]
    df = (
        df.sort_values("standard_value")
          .drop_duplicates(subset=["molecule_chembl_id", "assay_chembl_id"])
          .reset_index(drop=True)
    )
    logger.info(f"  Preprocessed pct records:  {len(df):,}")
    return df


# ?????????????????????????????????????????????????????????????????????????????
# INTERMEDIATE CHECK
# ?????????????????????????????????????????????????????????????????????????????

def has_no_intermediate(fp_a, fp_b, p_a, p_b, all_fps, all_p, sim_thresh=INTER_SIM_THRESH):
    p_lo, p_hi = sorted([p_a, p_b])
    mask = (all_p > p_lo) & (all_p < p_hi)
    if not mask.any():
        return True
    sims_a = bulk_tanimoto(fp_a, all_fps)
    sims_b = bulk_tanimoto(fp_b, all_fps)
    return not (mask & (sims_a >= sim_thresh) & (sims_b >= sim_thresh)).any()


# ?????????????????????????????????????????????????????????????????????????????
# CLIFF DETECTION
# ?????????????????????????????????????????????????????????????????????????????

def find_cliffs_xc50(
    df: pd.DataFrame,
    tanimoto_thresh: float = TANIMOTO_CUTOFF,
    delta_p: float = DELTA_PXIC50_CUTOFF,
    check_intermediate: bool = True,
    output_path: Optional[str] = None,
    checkpoint_path: Optional[str] = None,
    resume: bool = False,
) -> int:
    """
    Find activity cliff pairs within each assay.

    Performance
    -----------
    Compounds in each assay group are sorted by pXC50 descending.
    np.searchsorted locates the first index j where
    pXC50[i] - pXC50[j] >= delta_p in O(log n) per compound.
    BulkTanimotoSimilarity is then called ONLY on those candidates,
    not on the full O(n2) pair set.  For assays where the activity
    range is < delta_p (the majority), zero Tanimoto calls are made.

    Resilience
    ----------
    Results are appended to output_path after every assay.
    Completed assay IDs are written to checkpoint_path.
    Pass resume=True to skip already-finished assays on restart.

    Returns
    -------
    Total number of cliff pairs found.
    """
    # ?? Checkpoint: load done assays ????????????????????????????????????
    done_assays: set = set()
    if resume and checkpoint_path and os.path.exists(checkpoint_path):
        with open(checkpoint_path) as fh:
            done_assays = {ln.strip() for ln in fh if ln.strip()}
        logger.info(f"  Resuming: {len(done_assays):,} assays already completed")

    # ?? Output file header ???????????????????????????????????????????????
    # Write CSV header once; subsequent batches append without header.
    write_header = True
    if resume and output_path and os.path.exists(output_path):
        write_header = False   # file already has a header from previous run

    # ?? Assay scan ???????????????????????????????????????????????????????
    by_assay    = df.groupby("assay_chembl_id")
    total_assays = len(by_assay)
    n_pairs     = 0
    FLUSH_EVERY = 200   # write to disk every N assays

    logger.info(
        f"Scanning {total_assays:,} assays for activity cliffs "
        f"({len(done_assays):,} already done)?"
    )

    rows_buffer: list = []
    ckpt_fh = open(checkpoint_path, "a") if checkpoint_path else None

    try:
        for assay_id, grp in tqdm(by_assay, desc="Assay scan", unit="assay",
                                   total=total_assays):
            if assay_id in done_assays:
                continue

            grp = grp.reset_index(drop=True)
            n   = len(grp)
            if n < 2:
                # Mark done even if trivially skipped
                if ckpt_fh:
                    ckpt_fh.write(assay_id + "\n")
                continue

            # ?? Sort group by pXC50 descending for searchsorted trick ????
            order   = np.argsort(grp["pXC50"].values)[::-1]
            pvals   = grp["pXC50"].values[order].astype(float)
            fps     = [grp["fp"].iloc[k]            for k in order]
            mids    = [grp["molecule_chembl_id"].iloc[k] for k in order]
            ikeys   = [grp["inchikey"].iloc[k]      for k in order]
            pikeys  = [grp["parent_inchikey"].iloc[k] for k in order]
            neg_p   = -pvals   # ascending array for searchsorted

            for i in range(n):
                # Binary search: first j where pvals[i] - pvals[j] >= delta_p
                # i.e., -pvals[j] >= -pvals[i] + delta_p
                lo = int(np.searchsorted(neg_p, neg_p[i] + delta_p, side="left"))
                lo = max(lo, i + 1)   # avoid self-pair and double-counting
                if lo >= n:
                    continue           # no j satisfies delta_pXC50 criterion

                # ?? Compute Tanimoto ONLY for delta_pXC50-qualifying js ??
                cand_range = range(lo, n)
                cand_fps   = [fps[j] for j in cand_range]
                sims       = bulk_tanimoto(fps[i], cand_fps)

                for k, j in enumerate(cand_range):
                    sim_ij = sims[k]

                    # Identity guards
                    if mids[i] == mids[j]:
                        continue
                    if sim_ij >= 1.0 - 1e-9:
                        continue
                    if (ikeys[i]  is not None and ikeys[j]  is not None
                            and ikeys[i] == ikeys[j]):
                        continue
                    if (pikeys[i] is not None and pikeys[j] is not None
                            and pikeys[i] == pikeys[j]):
                        continue

                    if sim_ij < tanimoto_thresh:
                        continue

                    # Intermediate check (optional)
                    dp = pvals[i] - pvals[j]   # always >= delta_p here
                    if check_intermediate:
                        other_fps = fps[:i] + fps[i+1:j] + fps[j+1:]
                        other_p   = np.concatenate([pvals[:i],
                                                    pvals[i+1:j],
                                                    pvals[j+1:]])
                        if not has_no_intermediate(fps[i], fps[j],
                                                   pvals[i], pvals[j],
                                                   other_fps, other_p):
                            continue

                    ra = grp.iloc[order[i]]
                    ri = grp.iloc[order[j]]
                    rows_buffer.append({
                        "active_chembl_id":   ra["molecule_chembl_id"],
                        "active_smiles":      ra["canonical_smiles"],
                        "active_pXC50":       round(float(pvals[i]), 4),
                        "active_value_nM":    round(float(ra.get("value_nM", np.nan)), 3),
                        "active_std_type":    ra.get("standard_type", ""),
                        "inactive_chembl_id": ri["molecule_chembl_id"],
                        "inactive_smiles":    ri["canonical_smiles"],
                        "inactive_pXC50":     round(float(pvals[j]), 4),
                        "inactive_value_nM":  round(float(ri.get("value_nM", np.nan)), 3),
                        "inactive_std_type":  ri.get("standard_type", ""),
                        "delta_pXC50":        round(float(dp), 4),
                        "assay_chembl_id":    assay_id,
                        "target_chembl_id":   ra.get("target_chembl_id", "N/A"),
                        "target_name":        ra.get("target_name", "N/A"),
                        "target_organism":    ra.get("target_organism", "N/A"),
                        "assay_confidence_score": ra.get("assay_confidence_score", np.nan),
                        "_mol_active":        grp.iloc[order[i]]["mol"],
                        "_mol_inactive":      grp.iloc[order[j]]["mol"],
                        "cliff_type":         "XC50",
                    })

            # ?? Mark assay done ??????????????????????????????????????????
            if ckpt_fh:
                ckpt_fh.write(assay_id + "\n")

            # ?? Flush buffer to disk every FLUSH_EVERY assays ????????????
            if output_path and len(rows_buffer) >= FLUSH_EVERY:
                _flush_rows(rows_buffer, output_path, write_header)
                n_pairs     += len(rows_buffer)
                write_header = False
                rows_buffer  = []

        # ?? Final flush ??????????????????????????????????????????????????
        if output_path and rows_buffer:
            _flush_rows(rows_buffer, output_path, write_header)
            n_pairs += len(rows_buffer)
        elif not output_path:
            n_pairs = len(rows_buffer)

    finally:
        if ckpt_fh:
            ckpt_fh.flush()
            ckpt_fh.close()

    logger.info(f"Found {n_pairs:,} XC50 cliff pairs")
    return n_pairs


def _flush_rows(rows: list, path: str, write_header: bool) -> None:
    """Append a batch of cliff-pair dicts to the output CSV."""
    if not rows:
        return
    pd.DataFrame(rows).to_csv(
        path, mode="a", header=write_header, index=False, float_format="%.4f"
    )


def find_cliffs_pct_fallback(df_xc50, df_pct,
                              tanimoto_thresh=TANIMOTO_CUTOFF,
                              active_min=PXIC50_ACTIVE_MIN,
                              same_assay_type=False):
    """
    Pair active XC50 compounds with structurally similar % inhibition inactives.

    same_assay_type=True (--same-assay-type flag):
      Only pair compounds whose assays share the same assay_type code.
      Avoids cross-biology comparisons:
        B (binding: Ki/Kd/SPR)  paired only with B % inhibition screens
        F (functional: IC50/EC50 cell-based) paired only with F % inhibition
      This is a scientifically sound middle ground between the strict
      same-assay-id requirement and an unrestricted cross-assay comparison.
    """
    if df_xc50.empty or df_pct.empty:
        return pd.DataFrame()
    actives   = df_xc50[df_xc50["pXC50"] >= active_min].copy()
    inactives = df_pct[df_pct["standard_value"].astype(float) < PCT_INACTIVE_THRESH].copy()
    rows = []
    n_skipped_type = 0
    for target_id, ag in actives.groupby("target_chembl_id"):
        ig = inactives[inactives.get("target_chembl_id",
                        pd.Series(dtype=str)) == target_id]
        if ig.empty:
            continue
        ig_fps = ig["fp"].tolist()
        for _, ra in ag.iterrows():
            active_atype = str(ra.get("assay_type", "")).strip().upper()
            sims = bulk_tanimoto(ra["fp"], ig_fps)
            for k, (sim, (_, ri)) in enumerate(zip(sims, ig.iterrows())):
                if sim < tanimoto_thresh:
                    continue
                # Assay-type filter: skip cross-biology pairs when requested
                if same_assay_type:
                    inactive_atype = str(ri.get("assay_type", "")).strip().upper()
                    if (active_atype and inactive_atype
                            and active_atype != inactive_atype):
                        n_skipped_type += 1
                        continue
                rows.append({
                    "active_chembl_id":    ra["molecule_chembl_id"],
                    "active_smiles":       ra["canonical_smiles"],
                    "active_pXC50":        round(float(ra["pXC50"]), 4),
                    "active_value_nM":     round(float(ra.get("value_nM", np.nan)), 3),
                    "active_std_type":     ra.get("standard_type", ""),
                    "active_assay_type":   active_atype,
                    "inactive_chembl_id":  ri["molecule_chembl_id"],
                    "inactive_smiles":     ri["canonical_smiles"],
                    "inactive_pXC50":      float("nan"),
                    "inactive_value_nM":   float("nan"),
                    "inactive_std_type":   ri.get("standard_type", ""),
                    "inactive_pct_activity": float(ri["standard_value"]),
                    "inactive_assay_type": str(ri.get("assay_type", "")).strip().upper(),
                    "delta_pXC50":         float("nan"),
                    "assay_chembl_id":     ra.get("assay_chembl_id", "N/A"),
                    "target_chembl_id":    target_id,
                    "target_name":         ra.get("target_name", "N/A"),
                    "target_organism":     ra.get("target_organism", "N/A"),
                    "assay_confidence_score": ra.get("assay_confidence_score", np.nan),
                    "_mol_active":         ra["mol"],
                    "_mol_inactive":       ri["mol"],
                    "cliff_type":          "pct_fallback",
                })
    if same_assay_type:
        logger.info(f"  Assay-type filter removed {n_skipped_type:,} cross-type pairs")
    logger.info(f"Found {len(rows):,} fallback cliff pairs")
    return pd.DataFrame(rows)


# ?????????????????????????????????????????????????????????????????????????????
# ENRICHMENT  (tautomers, pH=7 charges)
# ?????????????????????????????????????????????????????????????????????????????

def enrich(df: pd.DataFrame, target_meta: dict) -> pd.DataFrame:
    """
    Attach tautomer flags and pH=7 charge summaries to cliff pairs.
    All mol-level ops are computed once per unique SMILES and merged back.
    """
    if df.empty:
        return df
    df = df.copy()

    # ?? Target metadata ???????????????????????????????????????????????????
    for field, key in [("uniprot_id",      "uniprot_id"),
                       ("target_name",     "target_name"),
                       ("target_organism", "organism")]:
        if field not in df.columns:
            df[field] = df["target_chembl_id"].map(
                lambda t, k=key: target_meta.get(t, {}).get(k, "N/A")
            )
        else:
            mask = df[field].isna() | (df[field] == "N/A")
            df.loc[mask, field] = df.loc[mask, "target_chembl_id"].map(
                lambda t, k=key: target_meta.get(t, {}).get(k, "N/A")
            )

    # ?? Per-unique-SMILES mol cache ???????????????????????????????????????
    all_active_smiles   = df["active_smiles"].dropna().unique().tolist()
    all_inactive_smiles = df["inactive_smiles"].dropna().unique().tolist()
    unique_smiles       = list(set(all_active_smiles + all_inactive_smiles))
    n_unique_active     = len(all_active_smiles)

    logger.info(f"  Mol-level enrichment: {n_unique_active} unique actives, "
                f"{len(unique_smiles)} total unique structures "
                f"({len(df):,} rows)?")

    mol_cache:    dict = {}
    taut_cache:   dict = {}
    charge_cache: dict = {}

    logger.info("  Tautomer enumeration + pH=7 charge states?")
    for smi in tqdm(unique_smiles, desc="Per-structure", unit="mol"):
        mol = smi_to_mol(smi)
        mol_cache[smi]    = mol
        taut_cache[smi]   = count_tautomers(mol)
        charge_cache[smi] = get_charge_ph7(mol)

    df["active_n_tautomers"]          = df["active_smiles"].map(taut_cache)
    df["inactive_n_tautomers"]        = df["inactive_smiles"].map(taut_cache)
    df["active_tautomer_important"]   = df["active_n_tautomers"]   >= 2
    df["inactive_tautomer_important"] = df["inactive_n_tautomers"] >= 2
    df["active_charge_ph7"]           = df["active_smiles"].map(charge_cache)
    df["inactive_charge_ph7"]         = df["inactive_smiles"].map(charge_cache)

    df.drop(columns=["_mol_active", "_mol_inactive"], inplace=True, errors="ignore")
    return df


# ?????????????????????????????????????????????????????????????????????????????
# OUTPUT
# ?????????????????????????????????????????????????????????????????????????????

FINAL_COLS = [
    "active_chembl_id", "active_smiles", "active_pXC50",
    "active_value_nM", "active_std_type",
    "inactive_chembl_id", "inactive_smiles", "inactive_pXC50",
    "inactive_value_nM", "inactive_std_type", "inactive_pct_activity",
    "delta_pXC50",
    "assay_chembl_id", "target_chembl_id", "target_name",
    "target_organism", "uniprot_id", "assay_confidence_score",
    "active_n_tautomers", "active_tautomer_important",
    "inactive_n_tautomers", "inactive_tautomer_important",
    "active_charge_ph7", "inactive_charge_ph7",
    "cliff_type",
]

def format_output(df: pd.DataFrame) -> pd.DataFrame:
    # Final safety net: drop any surviving same-molecule pair
    if "active_chembl_id" in df.columns and "inactive_chembl_id" in df.columns:
        same = df["active_chembl_id"] == df["inactive_chembl_id"]
        if same.any():
            logger.warning(f"  Safety filter removed {same.sum()} same-ChEMBL-ID pair(s)")
            df = df[~same]
    cols = [c for c in FINAL_COLS if c in df.columns]
    out = df[cols].copy().sort_values("delta_pXC50", ascending=False, na_position="last")
    return out.reset_index(drop=True)


# ?????????????????????????????????????????????????????????????????????????????
# DEMO DATA
# ?????????????????????????????????????????????????????????????????????????????

DEMO_DATA = [
    {"molecule_chembl_id": "CHEMBL_DEMO_A1",
     "canonical_smiles": "Cc1ccc(Nc2ncc3cc(-c4ccc(CC(=O)N5CCCC5)cc4)ccc3n2)cc1",
     "standard_value": 0.8, "standard_units": "nM", "standard_type": "IC50",
     "pchembl_value": 9.10, "assay_chembl_id": "CHEMBL_ASSAY_D1",
     "target_chembl_id": "CHEMBL279", "target_name": "EGFR",
     "target_organism": "Homo sapiens", "assay_confidence_score": 9},
    {"molecule_chembl_id": "CHEMBL_DEMO_A2",
     "canonical_smiles": "Cc1ccc(Nc2ncc3cc(-c4ccc(CC(=O)N5CCCC5)cc4)cnc3n2)cc1",
     "standard_value": 5100, "standard_units": "nM", "standard_type": "IC50",
     "pchembl_value": 5.29, "assay_chembl_id": "CHEMBL_ASSAY_D1",
     "target_chembl_id": "CHEMBL279", "target_name": "EGFR",
     "target_organism": "Homo sapiens", "assay_confidence_score": 9},
    {"molecule_chembl_id": "CHEMBL_DEMO_B1",
     "canonical_smiles": "O=C(Nc1ccc(-c2ccccc2)cc1)c1ccc(NC(=O)c2cccc(Cl)c2)cc1",
     "standard_value": 12, "standard_units": "nM", "standard_type": "IC50",
     "pchembl_value": 7.92, "assay_chembl_id": "CHEMBL_ASSAY_D2",
     "target_chembl_id": "CHEMBL2842", "target_name": "JAK1",
     "target_organism": "Homo sapiens", "assay_confidence_score": 9},
    {"molecule_chembl_id": "CHEMBL_DEMO_B2",
     "canonical_smiles": "O=C(Nc1ccc(-c2ccccc2)cc1)c1ccc(NC(=O)c2cccc(F)c2)cc1",
     "standard_value": 18500, "standard_units": "nM", "standard_type": "IC50",
     "pchembl_value": 4.73, "assay_chembl_id": "CHEMBL_ASSAY_D2",
     "target_chembl_id": "CHEMBL2842", "target_name": "JAK1",
     "target_organism": "Homo sapiens", "assay_confidence_score": 9},
    {"molecule_chembl_id": "CHEMBL_DEMO_C1",
     "canonical_smiles": "Cc1ccc(Nc2ncc3cc(-c4ccccc4)ccc3n2)cc1",
     "standard_value": 500, "standard_units": "nM", "standard_type": "IC50",
     "pchembl_value": 6.30, "assay_chembl_id": "CHEMBL_ASSAY_D1",
     "target_chembl_id": "CHEMBL279", "target_name": "EGFR",
     "target_organism": "Homo sapiens", "assay_confidence_score": 9},
]


def run_demo(args) -> pd.DataFrame:
    logger.info("=== DEMO MODE ? built-in example data ===")
    df_proc = preprocess_xc50(pd.DataFrame(DEMO_DATA))
    # Demo writes to a temp file; no checkpoint needed
    tmp_raw = "_demo_raw.csv"
    n = find_cliffs_xc50(
        df_proc,
        tanimoto_thresh=args.tanimoto,
        delta_p=args.delta_pxic50,
        check_intermediate=(not args.no_intermediate_check),
        output_path=tmp_raw,
        checkpoint_path=None,
        resume=False,
    )
    if n == 0 or not os.path.exists(tmp_raw):
        logger.warning("Demo: no cliffs found ? try --tanimoto 0.40 --delta-pxic50 1.5")
        return pd.DataFrame()
    df_cliffs = pd.read_csv(tmp_raw)
    df_cliffs["_mol_active"]   = df_cliffs["active_smiles"].map(smi_to_mol)
    df_cliffs["_mol_inactive"] = df_cliffs["inactive_smiles"].map(smi_to_mol)
    os.remove(tmp_raw)
    df_enriched = enrich(df_cliffs, {})
    return format_output(df_enriched)


# ?????????????????????????????????????????????????????????????????????????????
# CLI
# ?????????????????????????????????????????????????????????????????????????????

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = p.add_mutually_exclusive_group()
    src.add_argument(
        "--sqlite", metavar="PATH",
        help=(
            "Path to ChEMBL SQLite file (.db) or the directory containing it.\n"
            "Example: --sqlite C:\\path\\to\\chembl_37_sqlite\n"
            "         --sqlite C:\\path\\to\\chembl_37_sqlite\\chembl_37.db"
        ),
    )
    src.add_argument(
        "--input-csv", metavar="CSV",
        help="Pre-exported ChEMBL activity CSV (skips API/SQLite).",
    )
    src.add_argument(
        "--external-csv", metavar="CSV",
        help=(
            "External CSV with SMILES and pChEMBL values (non-ChEMBL data).\n"
            "Columns are auto-detected; use --smiles-col / --activity-col to override.\n"
            "The activity column must contain pChEMBL values (-log10 molar, e.g. 8.3),\n"
            "NOT raw IC50 concentrations in nM."
        ),
    )
    src.add_argument(
        "--demo", action="store_true",
        help="Run on built-in example data (no network, no SQLite required).",
    )
    # External CSV column mapping (only relevant when --external-csv is used)
    p.add_argument("--smiles-col",   metavar="COL", default=None,
                   help="SMILES column name in --external-csv (auto-detected if omitted).")
    p.add_argument("--activity-col", metavar="COL", default=None,
                   help="pChEMBL activity column in --external-csv (auto-detected if omitted).")
    p.add_argument("--id-col",       metavar="COL", default=None,
                   help="Compound ID column in --external-csv (auto-generated if omitted).")
    p.add_argument("--assay-col",    metavar="COL", default=None,
                   help=(
                       "Grouping column in --external-csv.  Compounds in different "
                       "groups are never paired.  If omitted, all compounds are treated "
                       "as one virtual assay."
                   ))
    p.add_argument("--target-col",   metavar="COL", default=None,
                   help="Target column in --external-csv (propagated to output).")
    p.add_argument("--target", "-t", metavar="CHEMBL_ID",
                   help="Restrict scan to a single ChEMBL target (e.g. CHEMBL279).")
    p.add_argument("--output", "-o", default="activity_cliffs.csv",
                   help="Output CSV path (default: activity_cliffs.csv).")
    p.add_argument("--tanimoto", type=float, default=TANIMOTO_CUTOFF,
                   help=f"ECFP4 Tanimoto cutoff (default {TANIMOTO_CUTOFF}).")
    p.add_argument("--delta-pxic50", type=float, default=DELTA_PXIC50_CUTOFF,
                   help=f"|?pXC50| cutoff (default {DELTA_PXIC50_CUTOFF}).")
    p.add_argument("--no-intermediate-check", action="store_true",
                   help="Skip Criterion 3 (faster, returns more pairs).")
    p.add_argument("--max-records", type=int, default=0,
                   help="Cap total rows fetched (0 = unlimited). Useful for first tests.")
    p.add_argument("--resume", action="store_true",
                   help=(
                       "Resume an interrupted run. Reads <output>.ckpt to skip "
                       "already-processed assays and appends to the existing CSV."
                   ))
    p.add_argument("--xc50-only", action="store_true",
                   help=(
                       "Skip pct_fallback pairs entirely. Only report pairs where "
                       "BOTH compounds have a measured XC50 in the SAME assay. "
                       "Strictest filter. Supersedes --same-assay-type."
                   ))
    p.add_argument("--same-assay-type", action="store_true",
                   help=(
                       "For pct_fallback pairs, require that the active XC50 assay "
                       "and the inactive %% inhibition assay share the same assay_type "
                       "(B=binding, F=functional, etc.). "
                       "Prevents cross-biology comparisons (e.g. Ki from SPR paired "
                       "with %% inhibition from a cell viability screen). "
                       "Recommended middle ground between --xc50-only and unrestricted."
                   ))
    return p.parse_args()


def print_summary(df: pd.DataFrame, path: str):
    logger.info("\n" + "?" * 60)
    logger.info(f"  Cliff pairs written          : {len(df):,}")
    logger.info(f"  Output                       : {path}")
    if not df.empty:
        d = df["delta_pXC50"].dropna()
        if len(d):
            logger.info(f"  ?pXC50  mean / max           : {d.mean():.2f} / {d.max():.2f}")
        logger.info(f"  Unique targets               : {df['target_chembl_id'].nunique()}")
        logger.info(f"  Unique assays                : {df['assay_chembl_id'].nunique()}")
        n_t = df.get("active_tautomer_important", pd.Series(False)).sum()
        logger.info(f"  Actives with >1 tautomer     : {n_t}")
    logger.info("?" * 60)


def main():
    args = parse_args()

    logger.info("?" * 60)
    logger.info("  Activity Cliff Scanner  v3  (ChEMBL/ECFP4)")
    logger.info("?" * 60)
    logger.info(f"  Confidence cutoff  : >= {CONFIDENCE_CUTOFF}")
    logger.info(f"  Tanimoto           : >= {args.tanimoto}")
    logger.info(f"  |?pXC50|           : >= {args.delta_pxic50}")
    logger.info(f"  Intermediate check : {not args.no_intermediate_check}")
    logger.info(f"  XC50-only mode     : {args.xc50_only}"
                + (" (pct_fallback pairs skipped)" if args.xc50_only else
                   " (pct_fallback pairs included)"))
    if not args.xc50_only:
        logger.info(f"  Same assay-type    : {args.same_assay_type}"
                    + (" (B-B / F-F only)" if args.same_assay_type else
                       " (cross-type allowed; use --same-assay-type to restrict)"))
    if args.target:
        logger.info(f"  Target filter      : {args.target}")

    # ?? Checkpoint paths ??????????????????????????????????????????????????
    checkpoint_path = args.output.replace(".csv", ".ckpt")
    if args.resume:
        logger.info(f"  Resume mode ON ? checkpoint: {checkpoint_path}")

    # ?? Demo ?????????????????????????????????????????????????????????????
    if args.demo:
        df_final = run_demo(args)
        df_final.to_csv(args.output, index=False, float_format="%.4f")
        print_summary(df_final, args.output)
        return

    # ?? Data source ???????????????????????????????????????????????????????
    target_meta: dict = {}

    if args.sqlite:
        # ?? SQLite path (recommended) ?????????????????????????????????????
        logger.info(f"\n[Step 1] Loading from SQLite: {args.sqlite}")
        db_path = find_sqlite_db(args.sqlite)

        df_raw = fetch_xc50_from_sqlite(
            db_path,
            target_id=args.target,
            max_records=args.max_records,
        )
        df_pct_raw = fetch_pct_from_sqlite(
            db_path,
            target_id=args.target,
            max_records=args.max_records // 4 if args.max_records > 0 else 0,
        )
        # Fetch target metadata from SQLite (no REST needed)
        logger.info("\n[Step 1b] Fetching target metadata from SQLite?")
        unique_targets = list(
            set(df_raw["target_chembl_id"].dropna().unique().tolist() +
                df_pct_raw["target_chembl_id"].dropna().unique().tolist()
                if not df_pct_raw.empty else [])
        )
        target_meta = fetch_target_meta_from_sqlite(db_path, unique_targets)

    elif args.input_csv:
        # Pre-exported ChEMBL CSV
        logger.info(f"\n[Step 1] Loading from CSV: {args.input_csv}")
        df_raw     = pd.read_csv(args.input_csv, low_memory=False)
        df_pct_raw = pd.DataFrame()

    elif args.external_csv:
        # External non-ChEMBL data (SMILES + pChEMBL)
        logger.info(f"\n[Step 1] Loading external CSV: {args.external_csv}")
        df_raw = load_external_csv(
            args.external_csv,
            smiles_col   = args.smiles_col,
            activity_col = args.activity_col,
            id_col       = args.id_col,
            assay_col    = args.assay_col,
            target_col   = args.target_col,
        )
        df_pct_raw = pd.DataFrame()   # no % inhibition fallback for external data

    else:
        # ?? ChEMBL REST API ???????????????????????????????????????????????
        logger.info("\n[Step 1] Querying ChEMBL REST API?")
        try:
            ping = _SESSION.get(f"{CHEMBL_API_BASE}/status.json", timeout=10)
            if ping.status_code != 200:
                raise ConnectionError(f"ChEMBL API HTTP {ping.status_code}")
            logger.info(f"  API status: {ping.json()}")
        except Exception as exc:
            logger.error(
                f"ChEMBL REST API unreachable: {exc}\n"
                "  Solutions:\n"
                "    1. --sqlite  C:\\path\\to\\chembl_37_sqlite\n"
                "    2. --input-csv  my_exported_activities.csv\n"
                "    3. --demo  (validate pipeline without data)\n"
                "  Try again later if EBI is temporarily down."
            )
            sys.exit(1)
        df_raw     = fetch_xc50_activities(args.target, args.max_records or 10_000)
        df_pct_raw = fetch_pct_activities(args.target, 5_000)

    if df_raw.empty:
        logger.error("No XC50 data loaded. Check your filters or data source.")
        sys.exit(1)

    # ?? Preprocess ????????????????????????????????????????????????????????
    logger.info("\n[Step 2] Preprocessing?")
    df_proc     = preprocess_xc50(df_raw)
    df_pct_proc = preprocess_pct(df_pct_raw) if not df_pct_raw.empty else pd.DataFrame()

    # ?? Cliff detection ???????????????????????????????????????????????????
    logger.info("\n[Step 3] Detecting activity cliffs?")
    # XC50 cliffs are written incrementally to disk; returns total pair count.
    n_xc50 = find_cliffs_xc50(
        df_proc,
        tanimoto_thresh=args.tanimoto,
        delta_p=args.delta_pxic50,
        check_intermediate=(not args.no_intermediate_check),
        output_path=checkpoint_path.replace(".ckpt", "_raw.csv"),
        checkpoint_path=checkpoint_path,
        resume=args.resume,
    )
    df_fallback = find_cliffs_pct_fallback(df_proc, df_pct_proc,
                                            tanimoto_thresh=args.tanimoto,
                                            same_assay_type=args.same_assay_type) \
                  if not args.xc50_only else pd.DataFrame()

    # Re-load raw XC50 results for enrichment.
    # IMPORTANT: load the file whenever it exists, not only when n_xc50 > 0.
    # On a --resume run that finishes the last batch, n_xc50 returns 0 (no NEW
    # pairs this session) but all previous pairs are already in _raw.csv.
    raw_csv = checkpoint_path.replace(".ckpt", "_raw.csv")
    if os.path.exists(raw_csv):
        df_xc50_cliffs = pd.read_csv(raw_csv, low_memory=False)
        logger.info(f"  Loaded {len(df_xc50_cliffs):,} XC50 cliff pairs from {raw_csv}")
        df_xc50_cliffs["_mol_active"]   = df_xc50_cliffs["active_smiles"].map(smi_to_mol)
        df_xc50_cliffs["_mol_inactive"] = df_xc50_cliffs["inactive_smiles"].map(smi_to_mol)
    else:
        df_xc50_cliffs = pd.DataFrame()
        logger.info("  No XC50 cliff pairs found (raw CSV absent).")

    # Summary before enrichment ? shows exactly how many of each type
    n_xc50_pairs     = len(df_xc50_cliffs)
    n_fallback_pairs = len(df_fallback)
    logger.info(
        f"\n  Pair type breakdown:\n"
        f"    XC50-XC50 pairs (both pXC50, same assay) : {n_xc50_pairs:>8,}\n"
        f"    Fallback pairs  (active pXC50 + %inh)    : {n_fallback_pairs:>8,}\n"
        f"    Total                                     : {n_xc50_pairs + n_fallback_pairs:>8,}\n"
        + (
            "\n  --xc50-only: fallback pairs skipped. All pairs are same-assay XC50-XC50.\n"
            if args.xc50_only else
            f"\n  Tip: re-run with --xc50-only to restrict to same-assay XC50-XC50 pairs only\n"
            f"  (strictest filter: both measurements under identical conditions).\n"
            f"  Current Tanimoto: {args.tanimoto} -- try 0.85 if XC50-XC50 = 0.\n"
        )
    )

    df_all = pd.concat([df_xc50_cliffs, df_fallback], ignore_index=True)

    if df_all.empty:
        logger.warning("No cliffs found. Consider --tanimoto 0.85 --delta-pxic50 1.5")
        return

    # ?? Fetch REST target metadata only if SQLite didn't supply it ????????
    if not args.sqlite and not args.input_csv and not args.external_csv:
        logger.info("\n[Step 4] Fetching target metadata (REST)?")
        unique_targets = df_all["target_chembl_id"].dropna().unique().tolist()
        target_meta = fetch_target_metadata(unique_targets)

    # ?? Enrich ????????????????????????????????????????????????????????????
    logger.info("\n[Step 5] Enriching (tautomers / pH=7 charges)?")
    df_enriched = enrich(df_all, target_meta)

    # ?? Export ????????????????????????????????????????????????????????????
    logger.info("\n[Step 6] Writing output?")
    df_final = format_output(df_enriched)
    df_final.to_csv(args.output, index=False, float_format="%.4f")
    print_summary(df_final, args.output)


if __name__ == "__main__":
    main()
