#!/usr/bin/env python3
"""
RECONSTRUCTION of the Fig. 5 transformation classifier (the original script was not in the project).
Validation on the old 229-pair set: ring-gain row reproduced exactly; other rows differ -> boundaries are approximate.

Per pair (MCS via cliff_dashboard.find_mcs; transformations from the ACTIVE compound's perspective):
  dn = (#atoms outside MCS in active) - (#atoms outside MCS in inactive)
  ring gain      : active has more rings
  stereo change  : same flat SMILES, different isomeric SMILES
  methyl +/-     : dn = +/-1 and the single extra atom is a carbon
  larger/smaller : any other dn > 0 / dn < 0 (ring loss kept apart)
  other          : dn = 0
Usage: python transformation_classifier_reconstruction.py pairs.csv out_prefix   (run from the act_cliff folder)
"""
import sys, re, os, multiprocessing as mp
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import cliff_dashboard as cd

def feat(a):
    i, r = a
    ma, mi = cd.smi_to_mol(r.active_smiles), cd.smi_to_mol(r.inactive_smiles)
    if ma is None or mi is None: return i, {}
    mcs = cd.find_mcs(ma, mi); pat = Chem.MolFromSmarts(mcs.smartsString) if mcs else None
    _, da = cd.get_atom_sets(ma, pat); _, di = cd.get_atom_sets(mi, pat)
    el = lambda m, d: "".join(sorted(m.GetAtomWithIdx(k).GetSymbol() for k in d))
    flat = lambda m: Chem.MolToSmiles(m, isomericSmiles=False)
    return i, dict(na=len(da), ni=len(di), el_a=el(ma, da), el_i=el(mi, di),
                   ra=rdMolDescriptors.CalcNumRings(ma), ri=rdMolDescriptors.CalcNumRings(mi),
                   same_flat=flat(ma) == flat(mi), same_iso=Chem.MolToSmiles(ma) == Chem.MolToSmiles(mi))

def classify(d):
    d = d.copy(); d["dn"] = d.na - d.ni
    sym = lambda s: Counter(re.findall("[A-Z][a-z]?", s if isinstance(s, str) else ""))
    def cls(r):
        if r.ra > r.ri: return "Ring gain in active"
        if r.same_flat and not r.same_iso: return "Stereochemistry change"
        xa, xi = sym(r.el_a) - sym(r.el_i), sym(r.el_i) - sym(r.el_a)
        if r.dn == 1:  return "Methyl addition in active" if dict(xa) == {"C": 1} else "Larger substituent in active"
        if r.dn == -1: return "Methyl removal in active" if dict(xi) == {"C": 1} else "Smaller substituent in active"
        if r.dn > 1:   return "Larger substituent in active"
        if r.dn < -1:  return "Smaller substituent in active" if r.ra >= r.ri else "Ring loss in active"
        return "Heteroatom / other change"
    d["tx"] = d.apply(cls, axis=1); return d

if __name__ == "__main__":
    src, pre = sys.argv[1], sys.argv[2]
    df = pd.read_csv(src)
    with mp.Pool(max(1, mp.cpu_count() - 1)) as p: res = dict(p.imap_unordered(feat, list(df.iterrows())))
    df = classify(pd.concat([df, pd.DataFrame([res[i] for i in df.index], index=df.index)], axis=1))
    tab = df.groupby("tx").agg(n=("tx", "size"), mean=("delta_pXC50", "mean"), median=("delta_pXC50", "median"),
                               pxc50_active=("active_pXC50", "mean"), pxc50_inactive=("inactive_pXC50", "mean"))
    tab["pct"] = 100 * tab.n / len(df)
    tab.sort_values("n", ascending=False).round(3).to_csv(pre + "_table.csv"); df.to_csv(pre + "_pairs.csv", index=False)
    print(tab.sort_values("n", ascending=False).round(3))
