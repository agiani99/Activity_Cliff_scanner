#!/usr/bin/env python3
"""
cliff_dashboard.py  --  Activity Cliff Visual Dashboard
========================================================
Reads an activity_cliffs CSV (output of activity_cliff_scanner.py) and
generates a self-contained HTML dashboard showing:

  * Side-by-side 2D depictions with MCS scaffold (light blue) and
    structurally differing atoms (green = active, red = inactive) highlighted
  * Lowest-energy MMFF conformer alignment on the MCS scaffold, RMSD reported
  * Diff-atom description: element symbols and fragment SMILES
  * Sortable / filterable DataTable with per-row detail modal

Usage
-----
  python cliff_dashboard.py --input all_cliffs.csv --output dashboard.html

  # Show only top 500 pairs by delta_pXC50
  python cliff_dashboard.py --input all_cliffs.csv --top-n 500

  # Filter to one target and skip 3D (faster)
  python cliff_dashboard.py --input all_cliffs.csv \
      --target CHEMBL279 --skip-3d --output EGFR_dashboard.html

  # Relax minimum delta_pXC50 to include fallback pairs
  python cliff_dashboard.py --input all_cliffs.csv --min-delta 0

Requirements
------------
  pip install rdkit pandas numpy tqdm
"""

import argparse
import html as html_mod
import os
import sys
import warnings
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
from tqdm import tqdm

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------

SVG_W, SVG_H         = 340, 230      # thumbnail dimensions (px)
SVG_W_LG, SVG_H_LG  = 560, 380      # modal (large) dimensions

# RGB tuples (0-1)
COLOR_ACTIVE_DIFF   = (0.10, 0.62, 0.10)   # green  -- active-unique atoms
COLOR_INACTIVE_DIFF = (0.82, 0.14, 0.14)   # red    -- inactive-unique atoms
COLOR_SCAFFOLD      = (0.72, 0.82, 0.95)   # pale blue -- MCS scaffold

N_CONFS  = 1        # conformers to generate per molecule (1 is fast and sufficient)
MCS_TIMEOUT = 30    # seconds for rdFMCS (large peptides / macrocycles need > 8 s)


# ---------------------------------------------------------------------------
# CHEMISTRY UTILITIES
# ---------------------------------------------------------------------------

def smi_to_mol(smi: str) -> Optional[Chem.Mol]:
    if not smi or not isinstance(smi, str):
        return None
    return Chem.MolFromSmiles(smi.strip())


def find_mcs(mol_a: Chem.Mol, mol_b: Chem.Mol):
    """
    Return rdFMCS result or None on failure / tiny MCS.

    Strategy (three passes in order of specificity):
      1. Strict  : completeRingsOnly + ringMatchesRingOnly — ideal for same-ring scaffolds
      2. Relaxed : no ring constraints — catches ring-size changes (cyclohexyl -> cyclopentyl)
      3. Both passes accept timed-out (canceled) results when numAtoms >= 3;
         rdFMCS returns the best partial MCS found so far on timeout, which is
         usually very close to the true MCS for drug-like molecules.
    """
    common_kwargs = dict(
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        timeout=MCS_TIMEOUT,
    )
    try:
        # Pass 1 — strict ring matching
        res = rdFMCS.FindMCS([mol_a, mol_b],
                             completeRingsOnly=True,
                             ringMatchesRingOnly=True,
                             **common_kwargs)
        if res.numAtoms >= 3:
            return res          # good result, canceled or not

        # Pass 2 — relaxed: handles ring-size changes (e.g. 6- vs 5-membered rings)
        res = rdFMCS.FindMCS([mol_a, mol_b],
                             completeRingsOnly=False,
                             ringMatchesRingOnly=False,
                             **common_kwargs)
        if res.numAtoms >= 3:
            return res

        return None
    except Exception:
        return None


def align_2d_coords(mol_a: Chem.Mol, mol_b: Chem.Mol, mcs_pat) -> bool:
    """
    Assign 2D coordinates to mol_a (standard layout) and mol_b aligned
    to mol_a via the MCS scaffold.  Modifies both molecules in-place.

    Uses the refPatt overload of GenerateDepictionMatching2DStructure
    (signature: mol, ref, confId, refPatt, acceptFailure) which lets
    RDKit handle the substructure matching internally.  This is safer
    than the atomMap overload, which triggers a C++ assertion error when:
      - the MCS search was canceled (partial result, wrong index ordering)
      - the relaxed MCS (completeRingsOnly=False) produces a pattern that
        matches a different atom subset than expected
      - mol_a and mol_b have different atom counts (very common)

    Returns True if scaffold alignment succeeded, False if independent
    2D coordinates were generated instead.
    """
    rdDepictor.Compute2DCoords(mol_a)

    if mcs_pat is None:
        rdDepictor.Compute2DCoords(mol_b)
        return False

    try:
        # refPatt overload: (mol, reference, confId=-1, refPatt=None, acceptFailure=False)
        rdDepictor.GenerateDepictionMatching2DStructure(
            mol_b, mol_a, -1, mcs_pat, True   # acceptFailure=True — never crash
        )
        return True
    except Exception:
        rdDepictor.Compute2DCoords(mol_b)
        return False
    """
    Returns (scaffold_atoms, diff_atoms) given a molecule and an MCS pattern.
    scaffold_atoms: indices in mol matching the MCS
    diff_atoms:     indices NOT in the MCS match
    """
    if mcs_pattern is None:
        return [], list(range(mol.GetNumAtoms()))
    match = mol.GetSubstructMatch(mcs_pattern)
    if not match:
        return [], list(range(mol.GetNumAtoms()))
    scaffold = list(match)
    diff     = [i for i in range(mol.GetNumAtoms()) if i not in match]
    return scaffold, diff


def get_atom_sets(mol: Chem.Mol, mcs_pattern) -> Tuple[List[int], List[int]]:
    """
    Returns (scaffold_atoms, diff_atoms) given a molecule and an MCS pattern.
    scaffold_atoms: indices in mol matching the MCS
    diff_atoms:     indices NOT in the MCS match
    """
    if mcs_pattern is None:
        return [], list(range(mol.GetNumAtoms()))
    match = mol.GetSubstructMatch(mcs_pattern)
    if not match:
        return [], list(range(mol.GetNumAtoms()))
    scaffold = list(match)
    diff     = [i for i in range(mol.GetNumAtoms()) if i not in match]
    return scaffold, diff


def diff_description(mol: Chem.Mol, diff_atoms: List[int]) -> str:
    """
    Short human-readable description of the atoms that differ.
    For 1-3 atoms: element symbols (e.g. 'N', 'Cl').
    For larger fragments: fragment SMILES.
    """
    if not diff_atoms:
        return "identical scaffold"
    if len(diff_atoms) <= 3:
        syms = sorted({mol.GetAtomWithIdx(i).GetSymbol() for i in diff_atoms})
        return ", ".join(syms)
    try:
        frag = Chem.MolFragmentToSmiles(mol, atomsToUse=diff_atoms)
        return frag if frag else f"{len(diff_atoms)} atoms"
    except Exception:
        return f"{len(diff_atoms)} atoms"


def mol_to_svg(mol: Chem.Mol,
               scaffold_atoms: List[int],
               diff_atoms: List[int],
               diff_color: tuple,
               w: int = SVG_W,
               h: int = SVG_H) -> str:
    """
    Generate an SVG of mol with:
      scaffold highlighted in pale blue, diff atoms in diff_color.
    Returns raw SVG string (no XML declaration).
    """
    if mol is None:
        return (f'<svg width="{w}" height="{h}" xmlns="http://www.w3.org/2000/svg">'
                f'<text x="10" y="20" fill="#999">N/A</text></svg>')
    try:
        drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
        opts   = drawer.drawOptions()
        opts.addStereoAnnotation = True
        opts.padding = 0.12

        atom_cols = {}
        atom_rads = {}
        bond_cols = {}
        hi_bonds  = []

        for idx in scaffold_atoms:
            atom_cols[idx] = COLOR_SCAFFOLD
            atom_rads[idx] = 0.25

        for idx in diff_atoms:
            atom_cols[idx] = diff_color
            atom_rads[idx] = 0.45

        for bond in mol.GetBonds():
            bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bi in diff_atoms and ei in diff_atoms:
                hi_bonds.append(bond.GetIdx())
                bond_cols[bond.GetIdx()] = diff_color
            elif bi in scaffold_atoms and ei in scaffold_atoms:
                hi_bonds.append(bond.GetIdx())
                bond_cols[bond.GetIdx()] = COLOR_SCAFFOLD

        rdMolDraw2D.PrepareMolForDrawing(mol, forceCoords=False)
        drawer.DrawMolecule(
            mol,
            highlightAtoms=list(atom_cols.keys()),
            highlightAtomColors=atom_cols,
            highlightBonds=hi_bonds,
            highlightBondColors=bond_cols,
            highlightAtomRadii=atom_rads,
        )
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        # Strip XML declaration for inline embedding
        if svg.startswith("<?xml"):
            svg = svg[svg.index("<svg"):]
        return svg
    except Exception as exc:
        return (f'<svg width="{w}" height="{h}" xmlns="http://www.w3.org/2000/svg">'
                f'<text x="10" y="20" fill="red">Draw error: {html_mod.escape(str(exc))}</text></svg>')


def compute_rmsd(mol_a: Chem.Mol,
                 mol_b: Chem.Mol,
                 mcs_pattern) -> Optional[float]:
    """
    Generate one ETKDGv3 + MMFF conformer per molecule and return the
    RMSD after alignment on MCS scaffold heavy atoms.
    Returns None if conformer generation or alignment fails.
    """
    try:
        ma_h = Chem.AddHs(Chem.RWMol(mol_a))
        mb_h = Chem.AddHs(Chem.RWMol(mol_b))

        ps = AllChem.ETKDGv3()
        ps.randomSeed  = 42
        ps.numThreads  = 0

        if AllChem.EmbedMolecule(ma_h, ps) == -1:
            return None
        if AllChem.EmbedMolecule(mb_h, ps) == -1:
            return None

        AllChem.MMFFOptimizeMolecule(ma_h, maxIters=500)
        AllChem.MMFFOptimizeMolecule(mb_h, maxIters=500)

        ma_3d = Chem.RemoveHs(ma_h)
        mb_3d = Chem.RemoveHs(mb_h)

        if mcs_pattern is not None:
            match_a = mol_a.GetSubstructMatch(mcs_pattern)
            match_b = mol_b.GetSubstructMatch(mcs_pattern)
            if match_a and match_b and len(match_a) >= 3:
                atom_map = list(zip(match_b, match_a))
                rmsd = rdMolAlign.AlignMol(mb_3d, ma_3d, atomMap=atom_map)
                return round(float(rmsd), 3)

        # Fallback: best-RMS without constraint
        rmsd = rdMolAlign.GetBestRMS(mb_3d, ma_3d)
        return round(float(rmsd), 3)

    except Exception:
        return None


# ---------------------------------------------------------------------------
# HTML GENERATION
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Activity Cliff Dashboard</title>
<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css" rel="stylesheet">
<link href="https://cdn.datatables.net/1.13.7/css/dataTables.bootstrap5.min.css" rel="stylesheet">
<style>
  body {{ font-size: .85rem; background: #f5f6fa; }}
  h1   {{ font-size: 1.4rem; font-weight: 700; }}
  .mol-thumb svg  {{ width: 220px; height: 148px; }}
  .mol-large svg  {{ width: 520px; height: 355px; }}
  td, th {{ vertical-align: middle !important; }}
  .badge-type {{ font-size: .68rem; }}
  .rmsd-pill  {{ font-size: .78rem; font-family: monospace; }}
  .diff-label {{ font-size: .72rem; color: #555; }}
  .stat-card  {{ border-left: 4px solid; padding: .7rem 1rem; }}
  .stat-val   {{ font-size: 1.5rem; font-weight: 700; line-height: 1; }}
  .stat-lbl   {{ font-size: .72rem; color: #666; text-transform: uppercase; }}
  .c-blue  {{ border-color: #0d6efd; }}
  .c-green {{ border-color: #198754; }}
  .c-amber {{ border-color: #ffc107; }}
  .c-red   {{ border-color: #dc3545; }}
  .legend-dot {{ display: inline-block; width: 12px; height: 12px;
                 border-radius: 50%; margin-right: 4px; }}
  .ld-scaffold {{ background: #b8d0f0; }}
  .ld-active   {{ background: #1a9e1a; }}
  .ld-inactive {{ background: #d12424; }}
  .table-hover tbody tr:hover {{ cursor: pointer; background: #eef3ff; }}
  #clifftable_filter input {{ border-radius: 20px; padding: .25rem .75rem; }}
</style>
</head>
<body>
<div class="container-fluid py-3">

  <div class="d-flex align-items-center mb-2 gap-3">
    <div>
      <h1>Activity Cliff Dashboard</h1>
      <p class="text-muted mb-0" style="font-size:.8rem">
        Source: <code>{input_file}</code> &nbsp;|&nbsp;
        Showing <strong>{n_shown}</strong> of <strong>{n_total:,}</strong> pairs
        &nbsp;|&nbsp; Tanimoto &ge; {tanimoto_min:.2f} &nbsp;|&nbsp; |&Delta;pXC50| &ge; {delta_min:.1f}
      </p>
    </div>
  </div>

  <!-- Summary stats -->
  <div class="row g-2 mb-3">
    <div class="col-6 col-md-3">
      <div class="card stat-card c-blue">
        <div class="stat-val">{n_shown}</div>
        <div class="stat-lbl">Pairs shown</div>
      </div>
    </div>
    <div class="col-6 col-md-3">
      <div class="card stat-card c-green">
        <div class="stat-val">{n_targets}</div>
        <div class="stat-lbl">Unique targets</div>
      </div>
    </div>
    <div class="col-6 col-md-3">
      <div class="card stat-card c-amber">
        <div class="stat-val">{mean_delta:.2f}</div>
        <div class="stat-lbl">Mean |&Delta;pXC50|</div>
      </div>
    </div>
    <div class="col-6 col-md-3">
      <div class="card stat-card c-red">
        <div class="stat-val">{pct_3d:.0f}%</div>
        <div class="stat-lbl">Pairs with RMSD</div>
      </div>
    </div>
  </div>

  <!-- Legend -->
  <div class="mb-2 text-muted" style="font-size:.78rem">
    <span class="legend-dot ld-scaffold"></span>MCS scaffold (common)&ensp;
    <span class="legend-dot ld-active"></span>Active-unique atoms (green)&ensp;
    <span class="legend-dot ld-inactive"></span>Inactive-unique atoms (red)
  </div>

  <!-- Table -->
  <div class="card shadow-sm">
    <div class="card-body p-0">
      <table id="clifftable" class="table table-hover table-bordered mb-0"
             style="width:100%">
        <thead class="table-dark">
          <tr>
            <th>#</th>
            <th>Active</th>
            <th>Inactive</th>
            <th data-bs-toggle="tooltip" title="log units difference in potency">|&Delta;pXC50|</th>
            <th data-bs-toggle="tooltip" title="RMSD (Angstrom) after scaffold alignment of MMFF conformers">RMSD (&Aring;)</th>
            <th>Active diff</th>
            <th>Inactive diff</th>
            <th>Target</th>
            <th>Conf.</th>
            <th>Type</th>
          </tr>
        </thead>
        <tbody>
{table_rows}
        </tbody>
      </table>
    </div>
  </div>
</div>

<!-- Detail modal -->
<div class="modal fade" id="detailModal" tabindex="-1">
  <div class="modal-dialog modal-xl">
    <div class="modal-content">
      <div class="modal-header py-2">
        <h5 class="modal-title" id="modalTitle">Cliff Pair Detail</h5>
        <button type="button" class="btn-close" data-bs-dismiss="modal"></button>
      </div>
      <div class="modal-body">
        <div class="row">
          <div class="col-6 text-center">
            <div class="fw-bold text-success mb-1" id="modalActiveLabel"></div>
            <div id="modalActiveSVG" class="mol-large"></div>
            <div class="mt-1 text-muted" id="modalActiveDiff"></div>
          </div>
          <div class="col-6 text-center">
            <div class="fw-bold text-danger mb-1" id="modalInactiveLabel"></div>
            <div id="modalInactiveSVG" class="mol-large"></div>
            <div class="mt-1 text-muted" id="modalInactiveDiff"></div>
          </div>
        </div>
        <hr>
        <div id="modalMeta" class="row g-2 text-center"></div>
      </div>
    </div>
  </div>
</div>

<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js"></script>
<script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
<script src="https://cdn.datatables.net/1.13.7/js/dataTables.bootstrap5.min.js"></script>
<script>
var PAIR_DATA = {pair_data_json};

$(document).ready(function() {{
  // Enable tooltips
  $('[data-bs-toggle="tooltip"]').each(function() {{
    new bootstrap.Tooltip(this);
  }});

  var dt = $('#clifftable').DataTable({{
    pageLength: 25,
    lengthMenu: [10, 25, 50, 100],
    order: [[3, 'desc']],
    columnDefs: [
      {{ targets: [1, 2], orderable: false }},
    ],
    language: {{ search: 'Filter:', zeroRecords: 'No matching cliff pairs found.' }},
  }});

  // Row click opens detail modal
  $('#clifftable tbody').on('click', 'tr', function() {{
    var rowIdx = dt.row(this).index();
    var d = PAIR_DATA[rowIdx];
    if (!d) return;

    $('#modalTitle').text('Pair #' + (rowIdx + 1) + '  \u2014  ' + d.target);
    $('#modalActiveLabel').html(
      d.active_id + '  pXC50 = ' + d.active_pxc50
    );
    $('#modalInactiveLabel').html(
      d.inactive_id + '  ' +
      (d.inactive_pxc50 !== null
        ? 'pXC50 = ' + d.inactive_pxc50
        : '%inh = ' + d.inactive_pct)
    );
    $('#modalActiveSVG').html(d.svg_active_lg);
    $('#modalInactiveSVG').html(d.svg_inactive_lg);
    $('#modalActiveDiff').html('<span class="text-success fw-bold">Active unique: </span>' + d.diff_active);
    $('#modalInactiveDiff').html('<span class="text-danger fw-bold">Inactive unique: </span>' + d.diff_inactive);

    var meta = [
      ['|&Delta;pXC50|', d.delta],
      ['RMSD (&Aring;)', d.rmsd !== null ? d.rmsd : 'n/a'],
      ['MCS atoms', d.mcs_atoms],
      ['Assay', d.assay],
      ['Conf. score', d.conf],
      ['Type', d.cliff_type],
    ];
    var html = '';
    meta.forEach(function(m) {{
      html += '<div class="col-6 col-md-2"><div class="border rounded p-1">' +
              '<div class="stat-lbl">' + m[0] + '</div>' +
              '<div class="fw-bold">' + m[1] + '</div></div></div>';
    }});
    $('#modalMeta').html(html);

    new bootstrap.Modal(document.getElementById('detailModal')).show();
  }});
}});
</script>
</body>
</html>
"""

_ROW_TEMPLATE = """\
          <tr>
            <td data-order="{rank}">{rank}</td>
            <td class="mol-thumb p-1">{svg_active}</td>
            <td class="mol-thumb p-1">{svg_inactive}</td>
            <td data-order="{delta_order}"><span class="badge {delta_badge_cls} fs-6">{delta_str}</span></td>
            <td>{rmsd_html}</td>
            <td class="diff-label text-success">{diff_active}</td>
            <td class="diff-label text-danger">{diff_inactive}</td>
            <td style="max-width:140px;word-wrap:break-word;">{target}</td>
            <td>{conf}</td>
            <td><span class="badge badge-type {type_cls}">{cliff_type}</span></td>
          </tr>"""


# ---------------------------------------------------------------------------
# PROCESSING
# ---------------------------------------------------------------------------

def process_pair(row: pd.Series, skip_3d: bool) -> Optional[dict]:
    """
    Process one cliff pair row.
    Returns a dict with SVGs, RMSD, diffs, or None on complete failure.
    """
    smi_a = row.get("active_smiles",   "")
    smi_b = row.get("inactive_smiles", "")

    mol_a = smi_to_mol(smi_a)
    mol_b = smi_to_mol(smi_b)
    if mol_a is None or mol_b is None:
        return None

    # MCS
    mcs     = find_mcs(mol_a, mol_b)
    mcs_pat = Chem.MolFromSmarts(mcs.smartsString) if mcs else None
    n_mcs   = mcs.numAtoms if mcs else 0

    scaffold_a, diff_a = get_atom_sets(mol_a, mcs_pat)
    scaffold_b, diff_b = get_atom_sets(mol_b, mcs_pat)

    diff_desc_a = diff_description(mol_a, diff_a)
    diff_desc_b = diff_description(mol_b, diff_b)

    # Align 2D layouts so the scaffold sits at the same position in both images
    align_2d_coords(mol_a, mol_b, mcs_pat)

    # Thumbnail SVGs
    svg_a_sm = mol_to_svg(mol_a, scaffold_a, diff_a, COLOR_ACTIVE_DIFF,   SVG_W,    SVG_H)
    svg_b_sm = mol_to_svg(mol_b, scaffold_b, diff_b, COLOR_INACTIVE_DIFF, SVG_W,    SVG_H)

    # Large SVGs for modal
    svg_a_lg = mol_to_svg(mol_a, scaffold_a, diff_a, COLOR_ACTIVE_DIFF,   SVG_W_LG, SVG_H_LG)
    svg_b_lg = mol_to_svg(mol_b, scaffold_b, diff_b, COLOR_INACTIVE_DIFF, SVG_W_LG, SVG_H_LG)

    # 3D RMSD
    rmsd = None
    if not skip_3d:
        rmsd = compute_rmsd(mol_a, mol_b, mcs_pat)

    return {
        "svg_active_sm":   svg_a_sm,
        "svg_inactive_sm": svg_b_sm,
        "svg_active_lg":   svg_a_lg,
        "svg_inactive_lg": svg_b_lg,
        "diff_active":     diff_desc_a,
        "diff_inactive":   diff_desc_b,
        "n_diff_active":   len(diff_a),
        "n_diff_inactive": len(diff_b),
        "n_mcs":           n_mcs,
        "rmsd":            rmsd,
    }


def build_dashboard(df: pd.DataFrame,
                    input_file: str,
                    output_file: str,
                    skip_3d: bool = False) -> None:

    n_total = len(df)
    import json

    table_rows  = []
    pair_data   = []   # JSON payload for modal

    for rank, (_, row) in enumerate(
        tqdm(df.iterrows(), total=len(df), desc="Processing pairs", unit="pair"),
        start=1
    ):
        result = process_pair(row, skip_3d)
        if result is None:
            continue

        delta      = row.get("delta_pXC50", np.nan)
        delta_str  = f"{delta:.2f}" if pd.notna(delta) else "n/a"
        delta_ord  = float(delta) if pd.notna(delta) else 0.0
        delta_badge_cls = "bg-primary" if pd.notna(delta) else "bg-secondary"

        rmsd       = result["rmsd"]
        rmsd_html  = (
            f'<span class="rmsd-pill badge bg-secondary">{rmsd:.3f}</span>'
            if rmsd is not None else
            '<span class="text-muted">—</span>'
        )

        target = str(row.get("target_name", row.get("target_chembl_id", "?")))
        if len(target) > 28:
            target = target[:26] + "…"

        cliff_type = str(row.get("cliff_type", "XC50"))
        type_cls   = "bg-info text-dark" if "pct" in cliff_type else "bg-success"

        conf  = str(row.get("assay_confidence_score", ""))
        assay = str(row.get("assay_chembl_id", ""))

        act_id  = str(row.get("active_chembl_id",   ""))
        ina_id  = str(row.get("inactive_chembl_id", ""))

        act_pxc = row.get("active_pXC50", np.nan)
        ina_pxc = row.get("inactive_pXC50", np.nan)
        ina_pct = row.get("inactive_pct_activity", np.nan)

        act_pxc_str = f"{act_pxc:.2f}" if pd.notna(act_pxc) else "n/a"
        ina_pxc_str = f"{ina_pxc:.2f}" if pd.notna(ina_pxc) else None
        ina_pct_str = f"{ina_pct:.1f}" if pd.notna(ina_pct) else "n/a"

        # Active cell label
        act_label = (f'<div style="font-size:.7rem">'
                     f'<b>{html_mod.escape(act_id)}</b><br>'
                     f'pXC50 = <b>{act_pxc_str}</b></div>')
        # Inactive cell label
        if ina_pxc_str:
            ina_potency = f'pXC50 = <b>{ina_pxc_str}</b>'
        else:
            ina_potency = f'%inh = <b>{ina_pct_str}</b>'
        ina_label = (f'<div style="font-size:.7rem">'
                     f'<b>{html_mod.escape(ina_id)}</b><br>'
                     f'{ina_potency}</div>')

        svg_active_cell   = result["svg_active_sm"]   + act_label
        svg_inactive_cell = result["svg_inactive_sm"] + ina_label

        table_rows.append(_ROW_TEMPLATE.format(
            rank            = rank,
            svg_active      = svg_active_cell,
            svg_inactive    = svg_inactive_cell,
            delta_str       = delta_str,
            delta_order     = delta_ord,
            delta_badge_cls = delta_badge_cls,
            rmsd_html       = rmsd_html,
            diff_active = html_mod.escape(result["diff_active"]),
            diff_inactive=html_mod.escape(result["diff_inactive"]),
            target      = html_mod.escape(target),
            conf        = conf,
            cliff_type  = cliff_type,
            type_cls    = type_cls,
        ))

        # JSON payload for modal (no SVG in JS — re-inject from hidden data attrs below)
        pair_data.append({
            "active_id":    act_id,
            "inactive_id":  ina_id,
            "active_pxc50": act_pxc_str,
            "inactive_pxc50": ina_pxc_str,
            "inactive_pct": ina_pct_str,
            "delta":        delta_str,
            "rmsd":         rmsd,
            "mcs_atoms":    result["n_mcs"],
            "diff_active":  result["diff_active"],
            "diff_inactive":result["diff_inactive"],
            "svg_active_lg":  result["svg_active_lg"],
            "svg_inactive_lg":result["svg_inactive_lg"],
            "target":  row.get("target_name", row.get("target_chembl_id", "")),
            "assay":   assay,
            "conf":    conf,
            "cliff_type": cliff_type,
        })

    # Summary stats
    delta_vals = df["delta_pXC50"].dropna()
    n_with_rmsd = sum(1 for p in pair_data if p["rmsd"] is not None)

    html_out = _HTML_TEMPLATE.format(
        input_file      = html_mod.escape(os.path.basename(input_file)),
        n_shown         = len(pair_data),
        n_total         = n_total,
        n_targets       = df["target_chembl_id"].nunique() if "target_chembl_id" in df.columns else "?",
        mean_delta      = float(delta_vals.mean()) if len(delta_vals) else 0.0,
        pct_3d          = 100 * n_with_rmsd / max(len(pair_data), 1),
        tanimoto_min    = df.get("tanimoto_ecfp4", pd.Series([0.95])).min() if "tanimoto_ecfp4" in df.columns else 0.95,
        delta_min       = float(delta_vals.min()) if len(delta_vals) else 0.0,
        table_rows      = "\n".join(table_rows),
        pair_data_json  = json.dumps(pair_data),
    )

    with open(output_file, "w", encoding="utf-8") as fh:
        fh.write(html_out)

    print(f"\n  Dashboard written: {output_file}")
    print(f"  Pairs rendered  : {len(pair_data):,}")
    print(f"  With 3D RMSD    : {n_with_rmsd:,}")
    print(f"  Open in browser : file://{os.path.abspath(output_file)}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",  "-i", required=True,
                   help="Activity cliff CSV (from activity_cliff_scanner.py).")
    p.add_argument("--output", "-o", default=None,
                   help="Output HTML file. Default: <input>_dashboard.html")
    p.add_argument("--top-n",  type=int, default=200,
                   help="Show top N pairs by |delta_pXC50| (default 200). Use 0 for all.")
    p.add_argument("--min-delta", type=float, default=2.0,
                   help="Minimum |delta_pXC50| to include (default 2.0).")
    p.add_argument("--target", default=None,
                   help="Filter to a single ChEMBL target ID (e.g. CHEMBL279).")
    p.add_argument("--cliff-type", choices=["XC50", "pct_fallback", "all"],
                   default="all",
                   help="Which cliff types to include (default: all).")
    p.add_argument("--skip-3d", action="store_true",
                   help="Skip 3D conformer generation / RMSD (faster, no conformer data).")
    return p.parse_args()


def main():
    args = parse_args()

    print(f"Reading {args.input} …")
    df = pd.read_csv(args.input, low_memory=False)
    print(f"  Loaded {len(df):,} pairs")

    # Filters
    if args.target:
        df = df[df["target_chembl_id"] == args.target]
        print(f"  After target filter ({args.target}): {len(df):,}")

    if args.cliff_type != "all" and "cliff_type" in df.columns:
        df = df[df["cliff_type"] == args.cliff_type]
        print(f"  After cliff_type filter ({args.cliff_type}): {len(df):,}")

    if args.min_delta > 0 and "delta_pXC50" in df.columns:
        # Keep rows where delta meets the threshold OR where delta is NaN
        # (NaN = pct_fallback pairs that have no computable delta_pXC50)
        df = df[(df["delta_pXC50"] >= args.min_delta) | df["delta_pXC50"].isna()]
        print(f"  After min_delta filter ({args.min_delta}): {len(df):,}")

    # Sort by delta_pXC50 descending (NaN last)
    if "delta_pXC50" in df.columns:
        df = df.sort_values("delta_pXC50", ascending=False, na_position="last")

    # Cap
    n_total = len(df)
    if args.top_n and args.top_n > 0:
        df = df.head(args.top_n)
        if len(df) < n_total:
            print(f"  Capped to top {args.top_n} by |delta_pXC50|")

    if df.empty:
        print("No pairs remain after filtering. Adjust --min-delta or --target.")
        sys.exit(1)

    output = args.output or args.input.replace(".csv", "_dashboard.html")
    print(f"\nGenerating dashboard for {len(df):,} pairs -> {output}")
    print(f"  3D RMSD: {'disabled' if args.skip_3d else 'enabled (ETKDGv3 + MMFF)'}")

    build_dashboard(df, args.input, output, skip_3d=args.skip_3d)


if __name__ == "__main__":
    main()
