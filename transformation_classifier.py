#!/usr/bin/env python3
"""
transformation_classifier.py -- rule-based structural-transformation label for each cliff pair.

Used by cliff_enrich.py (--classify-transformations) and runnable stand-alone:

    python transformation_classifier.py pairs.csv out_prefix

Method (all labels from the ACTIVE compound's perspective; every pair gets exactly ONE label):
  1. MCS of the two molecules (cliff_dashboard.find_mcs).  "Diff atoms" = atoms outside the MCS.
     dn = (#diff atoms in active) - (#diff atoms in inactive).
  2. Rules are applied in the fixed priority order below; the first that matches wins.

  P1  Ring gain / ring loss          number of rings differs
  P2  Stereochemistry change         identical flat SMILES, different isomeric SMILES
  P3  dn == 0 (same size)
        Halogen exchange             diff atoms are halogens on both sides, different halogen
        Heteroatom swap              diff atoms differ in element (e.g. C<->N, O<->S, C<->O)
        Other same-size change       same elements outside the MCS (isomer / position / bond-order)
  P4  dn != 0, exactly one extra atom (|dn| == 1)
        Methyl addition / removal    extra atom is a terminal carbon and the other side has no diff atoms (H <-> CH3)
        Ring expansion / contraction extra atom is a carbon and ring atoms lie outside the MCS on both sides
        Chain extension / shortening other extra carbon (CH2 in a chain)
        Halogen gain / loss          extra atom is a halogen
        Heteroatom gain / loss       extra atom is N, O, S, P, ...
  P5  dn != 0, |dn| >= 2
        Halogen gain / loss          only halogens differ
        Larger / smaller substituent otherwise
  Unclassifiable (MCS failed)        "Unclassified (no MCS)"

The classes are descriptive and approximate (MCS-based, heavy atoms only, not a formal MMP fragmentation).
"""
import os
import re
import sys
import multiprocessing as mp
from collections import Counter

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

HALOGENS = {"F", "Cl", "Br", "I"}


def _features(args):
    idx, active_smiles, inactive_smiles = args
    import cliff_dashboard as cd
    ma, mi = cd.smi_to_mol(active_smiles), cd.smi_to_mol(inactive_smiles)
    if ma is None or mi is None:
        return idx, {"mcs_ok": False}
    mcs = cd.find_mcs(ma, mi)
    pat = Chem.MolFromSmarts(mcs.smartsString) if mcs else None
    if pat is None:
        return idx, {"mcs_ok": False}
    _, da = cd.get_atom_sets(ma, pat)
    _, di = cd.get_atom_sets(mi, pat)
    el = lambda m, d: "".join(sorted(m.GetAtomWithIdx(k).GetSymbol() for k in d))
    ring = lambda m, d: sum(m.GetAtomWithIdx(k).IsInRing() for k in d)
    term = lambda m, d: sum(1 for k in d if m.GetAtomWithIdx(k).GetSymbol() == "C"
                            and m.GetAtomWithIdx(k).GetDegree() == 1)
    flat = lambda m: Chem.MolToSmiles(m, isomericSmiles=False)
    return idx, dict(mcs_ok=True, na=len(da), ni=len(di), el_a=el(ma, da), el_i=el(mi, di),
                     ring_a_diff=ring(ma, da), ring_i_diff=ring(mi, di),
                     term_a=term(ma, da), term_i=term(mi, di),
                     ra=rdMolDescriptors.CalcNumRings(ma), ri=rdMolDescriptors.CalcNumRings(mi),
                     same_flat=flat(ma) == flat(mi), same_iso=Chem.MolToSmiles(ma) == Chem.MolToSmiles(mi))


def _sym(s):
    return Counter(re.findall("[A-Z][a-z]?", s if isinstance(s, str) else ""))


def classify_row(r):
    """r: mapping with the feature fields.  Returns a single transformation label."""
    if not r.get("mcs_ok", False):
        return "Unclassified (no MCS)"
    if r["ra"] > r["ri"]:
        return "Ring gain in active"
    if r["ra"] < r["ri"]:
        return "Ring loss in active"
    if r["same_flat"] and not r["same_iso"]:
        return "Stereochemistry change"
    dn = r["na"] - r["ni"]
    xa, xi = _sym(r["el_a"]) - _sym(r["el_i"]), _sym(r["el_i"]) - _sym(r["el_a"])
    if dn == 0:
        sa, si = set(xa), set(xi)
        if sa and si and sa <= HALOGENS and si <= HALOGENS:
            return "Halogen exchange"
        if sa or si:
            return "Heteroatom swap"
        return "Other same-size change"
    sign_add = dn > 0
    extra = xa if sign_add else xi
    big_diff, small_n = (r["na"], r["ni"]) if sign_add else (r["ni"], r["na"])
    term_big = r["term_a"] if sign_add else r["term_i"]
    side = "in active"
    if abs(dn) == 1:
        e = next(iter(extra), None)
        if e == "C":
            if small_n == 0 and term_big == 1:
                return "Methyl addition in active" if sign_add else "Methyl removal in active"
            if r["ring_a_diff"] > 0 and r["ring_i_diff"] > 0:
                return "Ring expansion in active" if sign_add else "Ring contraction in active"
            return "Chain extension in active" if sign_add else "Chain shortening in active"
        if e in HALOGENS:
            return "Halogen gain in active" if sign_add else "Halogen loss in active"
        return "Heteroatom gain in active" if sign_add else "Heteroatom loss in active"
    if extra and set(extra) <= HALOGENS:
        return "Halogen gain in active" if sign_add else "Halogen loss in active"
    return "Larger substituent in active" if sign_add else "Smaller substituent in active"


def classify_pairs(df, active_col="active_smiles", inactive_col="inactive_smiles", processes=None):
    """Return a Series of labels aligned to df.index (features computed in parallel)."""
    jobs = [(i, a, b) for i, a, b in zip(df.index, df[active_col], df[inactive_col])]
    n = processes or max(1, mp.cpu_count() - 1)
    if n > 1 and len(jobs) > 8:
        with mp.Pool(n) as p:
            res = dict(p.imap_unordered(_features, jobs, chunksize=4))
    else:
        res = dict(map(_features, jobs))
    feats = pd.DataFrame([res[i] for i in df.index], index=df.index)
    return feats.apply(classify_row, axis=1), feats


def summary_table(df, label_col, delta_col="delta_pXC50"):
    t = df.groupby(label_col).agg(n=(label_col, "size"), mean=(delta_col, "mean"), median=(delta_col, "median"))
    t["pct"] = 100 * t.n / len(df)
    return t.sort_values("n", ascending=False).round(3)


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    src, pre = sys.argv[1], sys.argv[2]
    df = pd.read_csv(src)
    labels, feats = classify_pairs(df)
    df["transformation_type"] = labels
    t = summary_table(df, "transformation_type")
    t.to_csv(pre + "_table.csv")
    pd.concat([df, feats], axis=1).to_csv(pre + "_pairs.csv", index=False)
    print(t)
