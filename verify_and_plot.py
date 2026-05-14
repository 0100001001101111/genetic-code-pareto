#!/usr/bin/env python3
"""
Verification Script: Genetic Code Pareto Analysis
==================================================
Self-contained. All data defined inline. No external files loaded.
Reproduces every number in the paper from raw genetic code definitions.

Usage: python3 verify_and_plot.py
Output: Figure 1, Table 1, all claim-relevant numbers to stdout.

Brothers, B. (2026). Project Aletheia.
"""

import numpy as np
import os
import sys
import time
import hashlib
import platform
import subprocess
from pathlib import Path
from collections import defaultdict
from scipy import stats

START_TIME = time.time()
np.random.seed(42)
os.makedirs('output', exist_ok=True)

def file_sha256(path):
    path = Path(path)
    if not path.exists():
        return "MISSING"
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def git_commit():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"

def package_version(name):
    try:
        import importlib.metadata
        return importlib.metadata.version(name)
    except Exception:
        return "not-installed"

# ================================================================
# RAW DATA: Amino acid properties
# Sources: Woese 1966 (polar req), Kyte-Doolittle 1982 (hydropathy),
#          Grantham 1974 (volume), Akashi-Gojobori 2002 (ATP cost)
# Format: (polar_requirement, hydropathy, volume, cost_ATP)
# ================================================================
AA_PROPS = {
    'F': (5.0, 2.8, 132, 52), 'L': (4.9, 3.8, 111, 27), 'I': (4.9, 4.5, 111, 32),
    'M': (5.3, 1.9, 105, 34), 'V': (5.6, 4.2, 105, 23), 'S': (7.5, -0.8, 99, 11),
    'P': (6.6, -1.6, 112, 12.5), 'T': (6.6, -0.7, 119, 18.5), 'A': (7.0, 1.8, 31, 11.5),
    'Y': (5.4, -1.3, 141, 50), 'H': (8.4, -3.2, 96, 38.5), 'Q': (8.6, -3.5, 114, 26.5),
    'N': (10.0, -3.5, 56, 14.5), 'K': (10.1, -3.9, 119, 30.5), 'D': (13.0, -3.5, 54, 12.5),
    'E': (12.5, -3.5, 83, 15.5), 'C': (4.8, 2.5, 55, 24.5), 'W': (5.2, -0.9, 170, 74.5),
    'R': (9.1, -4.5, 148, 27.5), 'G': (7.9, -0.4, 3, 11.5), '*': (0, 0, 0, 0),
}

# ================================================================
# RAW DATA: Standard genetic code (NCBI translation table 1)
# ================================================================
STANDARD = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G','AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
}

# ================================================================
# RAW DATA: All 23 NCBI variant codes (differences from standard)
# Source: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
# ================================================================
VARIANTS = {
    2:  ('Vertebrate Mito', {'AGA':'*','AGG':'*','ATA':'M','TGA':'W'}, 'mito'),
    3:  ('Yeast Mito', {'ATA':'M','CTT':'T','CTC':'T','CTA':'T','CTG':'T','TGA':'W'}, 'mito'),
    4:  ('Mycoplasma', {'TGA':'W'}, 'mito'),
    5:  ('Invertebrate Mito', {'AGA':'S','AGG':'S','ATA':'M','TGA':'W'}, 'mito'),
    6:  ('Ciliate Nuclear', {'TAA':'Q','TAG':'Q'}, 'nuclear'),
    9:  ('Echinoderm Mito', {'AAA':'N','AGA':'S','AGG':'S','TGA':'W'}, 'mito'),
    10: ('Euplotid Nuclear', {'TGA':'C'}, 'nuclear'),
    12: ('Alt Yeast Nuclear', {'CTG':'S'}, 'nuclear'),
    13: ('Ascidian Mito', {'AGA':'G','AGG':'G','ATA':'M','TGA':'W'}, 'mito'),
    14: ('Alt Flatworm Mito', {'AAA':'N','AGA':'S','AGG':'S','TAA':'Y','TGA':'W'}, 'mito'),
    15: ('Blepharisma Nuclear', {'TAG':'Q'}, 'nuclear'),
    16: ('Chlorophycean Mito', {'TAG':'L'}, 'mito'),
    21: ('Trematode Mito', {'AGA':'S','AGG':'S','ATA':'M','TGA':'W','AAA':'N'}, 'mito'),
    22: ('Scenedesmus Mito', {'TCA':'*','TAG':'L'}, 'mito'),
    23: ('Thraustochytrium Mito', {'TTA':'*'}, 'mito'),
    24: ('Rhabdopleuridae Mito', {'AGA':'S','AGG':'K','TGA':'W'}, 'mito'),
    25: ('Cand. Div. SR1', {'TGA':'G'}, 'nuclear'),
    26: ('Pachysolen Nuclear', {'CTG':'A'}, 'nuclear'),
    27: ('Karyorelictid Nuclear', {'TAA':'Q','TAG':'Q','TGA':'W'}, 'nuclear'),
    29: ('Mesodinium Nuclear', {'TAA':'Y','TAG':'Y'}, 'nuclear'),
    30: ('Peritrich Nuclear', {'TAA':'E','TAG':'E'}, 'nuclear'),
    31: ('Blastocrithidia Nucl.', {'TGA':'W','TAG':'E','TAA':'E'}, 'nuclear'),
    33: ('Cephalodiscidae Mito', {'AGA':'S','AGG':'K','TGA':'W','TAA':'Y'}, 'mito'),
}

# ================================================================
# SCORING FUNCTIONS
# ================================================================
def make_code(changes):
    code = dict(STANDARD)
    for c, aa in changes.items():
        code[c] = aa
    return code

def error_tolerance(code, prop_idx=0):
    """Mean absolute property distance across all non-synonymous single-base mutations."""
    bases = 'TCAG'
    total = 0; n = 0
    for codon, aa in code.items():
        if aa == '*' or aa not in AA_PROPS: continue
        for pos in range(3):
            for b in bases:
                if b == codon[pos]: continue
                mut = codon[:pos] + b + codon[pos+1:]
                maa = code.get(mut, '*')
                if maa == '*' or maa not in AA_PROPS or aa == maa: continue
                total += abs(AA_PROPS[aa][prop_idx] - AA_PROPS[maa][prop_idx])
                n += 1
    return total / max(n, 1)

def cost_efficiency(code):
    """Spearman correlation: codon count per AA vs biosynthetic cost (ATP)."""
    ac = defaultdict(int)
    for c, aa in code.items():
        if aa != '*': ac[aa] += 1
    return stats.spearmanr([AA_PROPS[a][3] for a in ac], [ac[a] for a in ac])[0]

def wobble_synonymy(code):
    """Fraction of 3rd-position mutations that are synonymous."""
    bases = 'TCAG'
    syn = 0; total = 0
    for codon, aa in code.items():
        if aa == '*': continue
        for b in bases:
            if b == codon[2]: continue
            total += 1
            if code.get(codon[:2] + b, '*') == aa: syn += 1
    return syn / total * 100

def codon_wobble_synonymy(code, codon):
    """Third-position synonymy for a single sense codon."""
    bases = 'TCAG'
    aa = code[codon]
    if aa == '*':
        return None
    syn = 0; total = 0
    for b in bases:
        if b == codon[2]:
            continue
        total += 1
        if code.get(codon[:2] + b, '*') == aa:
            syn += 1
    return syn / total * 100

def codon_wobble_reassignment_correlation():
    """Across standard sense codons, compare wobble buffering with reassignment load."""
    variant_codes = [make_code(changes) for _, changes, _ in VARIANTS.values()]
    all_codes = [STANDARD] + variant_codes
    wobble_by_codon = []
    reassignments_by_codon = []
    for codon, aa in STANDARD.items():
        if aa == '*':
            continue
        values = [codon_wobble_synonymy(code, codon) for code in all_codes]
        values = [v for v in values if v is not None]
        wobble_by_codon.append(sum(values) / len(values))
        reassignments_by_codon.append(sum(1 for code in variant_codes if code[codon] != aa))
    return stats.spearmanr(wobble_by_codon, reassignments_by_codon)

def compatibility_burden(changes):
    """Summarize reassignment burden relative to the standard code."""
    stop_to_sense = sum(1 for codon, aa in changes.items() if STANDARD[codon] == '*' and aa != '*')
    sense_to_stop = sum(1 for codon, aa in changes.items() if STANDARD[codon] != '*' and aa == '*')
    sense_to_sense = sum(1 for codon, aa in changes.items() if STANDARD[codon] != '*' and aa != '*')
    return {
        "reassigned_codons": len(changes),
        "stop_to_sense": stop_to_sense,
        "sense_to_stop": sense_to_stop,
        "sense_to_sense": sense_to_sense,
        "stop_sense_boundary": stop_to_sense + sense_to_stop,
        "changes": "; ".join(f"{codon}:{STANDARD[codon]}>{aa}" for codon, aa in sorted(changes.items())),
    }

# ================================================================
# SCORE ALL 24 CODES
# ================================================================
print("=" * 60)
print("GENETIC CODE PARETO ANALYSIS — Step 21 Verification")
print("=" * 60)

results = []
std_et = error_tolerance(STANDARD, 0)
std_et_h = error_tolerance(STANDARD, 1)
std_et_v = error_tolerance(STANDARD, 2)
std_ce = cost_efficiency(STANDARD)
std_wb = wobble_synonymy(STANDARD)
results.append(('Standard', 1, 'standard', std_et, std_et_h, std_et_v, std_ce, std_wb, 0))

for vid, (name, changes, vtype) in VARIANTS.items():
    code = make_code(changes)
    results.append((name, vid, vtype, error_tolerance(code, 0),
                     error_tolerance(code, 1), error_tolerance(code, 2),
                     cost_efficiency(code), wobble_synonymy(code), len(changes)))

# ================================================================
# RANDOM BASELINE (1,000,000 codes, seed 42)
# ================================================================
print("\nGenerating 1,000,000 random codes (seed 42)...")
blocks = defaultdict(list)
for c, aa in STANDARD.items():
    if aa != '*': blocks[aa].append(c)
bl = list(blocks.values())
ba = list(blocks.keys())

N_RANDOM = 1_000_000
rand_et = []; rand_et_h = []; rand_et_v = []; rand_ce = []
for i in range(N_RANDOM):
    perm = np.random.permutation(ba)
    rc = dict(STANDARD)
    for j, block in enumerate(bl):
        for c in block: rc[c] = perm[j]
    rand_et.append(error_tolerance(rc, 0))
    rand_et_h.append(error_tolerance(rc, 1))
    rand_et_v.append(error_tolerance(rc, 2))
    rand_ce.append(cost_efficiency(rc))
    if (i+1) % 100000 == 0:
        print(f"  {i+1}/{N_RANDOM} done")

# ================================================================
# PRINT ALL CLAIM-RELEVANT NUMBERS
# ================================================================
print(f"\n{'='*60}")
print("ALL CLAIM-RELEVANT NUMBERS")
print(f"{'='*60}")

print(f"\nStandard code:")
print(f"  ET (polar req):  {std_et:.4f}")
print(f"  Cost (Spearman): {std_ce:.4f}")
print(f"  Wobble syn:      {std_wb:.1f}%")

et_pctile = sum(1 for x in rand_et if x > std_et) / N_RANDOM * 100
ce_pctile = sum(1 for x in rand_ce if x > std_ce) / N_RANDOM * 100
print(f"  ET percentile:   {et_pctile:.1f}%")
print(f"  Cost percentile: {ce_pctile:.1f}%")

print(f"  ET (hydropathy): {std_et_h:.4f}")
print(f"  ET (volume):     {std_et_v:.4f}")
std_et_h_pctile = sum(1 for x in rand_et_h if x > std_et_h) / N_RANDOM * 100
std_et_v_pctile = sum(1 for x in rand_et_v if x > std_et_v) / N_RANDOM * 100
print(f"  Hydropathy percentile: {std_et_h_pctile:.1f}%")
print(f"  Volume percentile:     {std_et_v_pctile:.1f}%")

# Dominators
doms = [(n, vid, et, ce) for n, vid, vt, et, eth, etv, ce, wb, nc in results
        if et < std_et and ce < std_ce and vid != 1]
print(f"\nDominators of standard: {len(doms)}")
for n, vid, et, ce in sorted(doms, key=lambda x: x[2]):
    print(f"  {n:25s} (table {vid:2d}): ET={et:.4f}, Cost={ce:.4f}")

# Pareto frontier
pareto = []
for r in results:
    if not any(r2[3] < r[3] and r2[6] < r[6] for r2 in results if r2[1] != r[1]):
        pareto.append(r)
print(f"\nPareto frontier ({len(pareto)} codes):")
for n, vid, vt, et, eth, etv, ce, wb, nc in sorted(pareto, key=lambda x: x[3]):
    print(f"  {n:25s}: ET={et:.4f}, Cost={ce:.4f}")

# Trade-off
et_ch = [et - std_et for _, vid, _, et, eth, etv, ce, wb, nc in results if vid != 1]
ce_ch = [ce - std_ce for _, vid, _, et, eth, etv, ce, wb, nc in results if vid != 1]
tau, p = stats.kendalltau(et_ch, ce_ch)
print(f"\nTrade-off: Kendall tau = {tau:.4f}, p = {p:.4f}")

# Wobble range
wbs = [wb for _, _, _, _, _, _, _, wb, _ in results]
print(f"\nWobble range: {min(wbs):.1f}% - {max(wbs):.1f}%")
wobble_codon_rho, wobble_codon_p = codon_wobble_reassignment_correlation()
print(f"Per-codon wobble vs reassignment count: Spearman rho = {wobble_codon_rho:.4f}, p = {wobble_codon_p:.4f}")
wobble_code_rho, wobble_code_p = stats.spearmanr(
    [wb for _, _, _, _, _, _, _, wb, _ in results],
    [nc for _, _, _, _, _, _, _, _, nc in results],
)
print(f"Code-level wobble vs changed-codon count: Spearman rho = {wobble_code_rho:.4f}, p = {wobble_code_p:.4f}")

# Mycoplasma
myco = [r for r in results if r[1] == 4][0]
print(f"\nMycoplasma: ET={myco[3]:.4f}, Cost={myco[6]:.4f}")

# Compatibility burden
compatibility = {
    1: {
        "name": "Standard",
        "type": "standard",
        "context": "nuclear/universal",
        "dominates": False,
        **compatibility_burden({}),
    }
}
for vid, (name, changes, vtype) in VARIANTS.items():
    compatibility[vid] = {
        "name": name,
        "type": vtype,
        "context": "organelle/reduced" if vtype == "mito" else "nuclear variant",
        "dominates": any(d[1] == vid for d in doms),
        **compatibility_burden(changes),
    }

dominator_burdens = [compatibility[vid] for _, vid, _, _ in doms]
frontier_ids = {vid for _, vid, _, _, _, _, _, _, _ in pareto}
frontier_burdens = [compatibility[vid] for vid in frontier_ids]
print("\nCompatibility burden:")
print(f"  Dominators in organelle/reduced context: {sum(1 for b in dominator_burdens if b['type'] == 'mito')}/{len(dominator_burdens)}")
print(f"  Pareto frontier in organelle/reduced context: {sum(1 for b in frontier_burdens if b['type'] == 'mito')}/{len(frontier_burdens)}")
print(f"  Dominators crossing stop/sense boundary: {sum(1 for b in dominator_burdens if b['stop_sense_boundary'] > 0)}/{len(dominator_burdens)}")
print(f"  Nuclear dominators with no stop/sense boundary change: {sum(1 for b in dominator_burdens if b['type'] == 'nuclear' and b['stop_sense_boundary'] == 0)}/{len(dominator_burdens)}")

print(f"\nInvocation: {' '.join(sys.argv)}")
print(f"Git commit: {git_commit()}")
print(f"Platform: {platform.platform()}")
print(f"Software: Python {sys.version.split()[0]}, scipy {package_version('scipy')}, numpy {np.__version__}, pandas {package_version('pandas')}, matplotlib {package_version('matplotlib')}")
print(f"Random codes: {N_RANDOM}, Seed: 42")

# ================================================================
# GENERATE FIGURE 1
# ================================================================
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Random cloud
    ax.scatter(rand_et[::10], rand_ce[::10], c='lightgray', s=2, alpha=0.3, zorder=1)

    # Biological codes
    colors = {'mito': '#2196F3', 'nuclear': '#FF9800', 'standard': '#E91E63'}
    for n, vid, vt, et, eth, etv, ce, wb, nc in results:
        c = colors[vt]
        on_f = any(r[1] == vid for r in pareto)
        ax.scatter(et, ce,
                   c=c,
                   marker='*' if vt == 'standard' else ('o' if vt == 'mito' else 's'),
                   s=200 if vt == 'standard' else 80,
                   zorder=10 if on_f else 5,
                   edgecolors='black' if on_f else 'none',
                   linewidths=2 if on_f else 0)

    # Pareto frontier line
    ps = sorted(pareto, key=lambda x: x[3])
    ax.plot([p[3] for p in ps], [p[6] for p in ps], 'k--', lw=1, alpha=0.5, zorder=4)

    # Labels
    for p in pareto:
        ax.annotate(p[0], (p[3], p[6]), xytext=(p[3]+0.015, p[6]+0.008),
                    fontsize=8, fontstyle='italic')
    ax.annotate('Standard', (std_et, std_ce), xytext=(std_et+0.015, std_ce-0.015),
                fontsize=9, fontweight='bold', color='#E91E63')

    ax.set_xlabel('Error tolerance (mean polar requirement distance)\n'
                  '\u2190 better                                                          worse \u2192',
                  fontsize=11)
    ax.set_ylabel('Biosynthetic cost efficiency (Spearman \u03c1)\n\u2190 more efficient',
                  fontsize=11)
    ax.set_title('Pareto Frontier of Biological Genetic Codes', fontsize=13, fontweight='bold')

    legend = [
        Line2D([0],[0], marker='*', color='w', markerfacecolor='#E91E63', markersize=14,
               label='Standard code'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#2196F3', markersize=10,
               label='Mitochondrial variants'),
        Line2D([0],[0], marker='s', color='w', markerfacecolor='#FF9800', markersize=10,
               label='Nuclear variants'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='lightgray', markersize=8,
               label=f'Random codes (n={N_RANDOM//1000}K)'),
        Line2D([0],[0], linestyle='--', color='black', alpha=0.5, label='Pareto frontier'),
    ]
    ax.legend(handles=legend, loc='upper right', fontsize=9)

    plt.tight_layout()
    fig.savefig('output/figure1.png', dpi=300, bbox_inches='tight')
    fig.savefig('output/figure1.pdf', bbox_inches='tight')
    fig.savefig('output/figure1.tiff', dpi=300, bbox_inches='tight')
    print("\nFigure saved: output/figure1.png, .pdf, .tiff (300 DPI)")
    plt.close()
except ImportError:
    print("\nmatplotlib not available; figure not generated")

# ================================================================
# GENERATE TABLE 1
# ================================================================
results_sorted = sorted(results, key=lambda x: x[3])
lines = [
    "# Table 1: Complete scoring of all 24 biological genetic codes\n",
    "| ID | Name | Type | #Chg | ET (PR) | Delta% | ET %ile | ET hydropathy %ile | ET volume %ile | Cost (rho) | Cost %ile | Wobble | Pareto |",
    "|---:|:-----|:-----|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|",
]
for n, vid, vt, et, eth, etv, ce, wb, nc in results_sorted:
    delta = (et - std_et) / std_et * 100
    etp = sum(1 for x in rand_et if x > et) / N_RANDOM * 100
    ethp = sum(1 for x in rand_et_h if x > eth) / N_RANDOM * 100
    etvp = sum(1 for x in rand_et_v if x > etv) / N_RANDOM * 100
    cep = sum(1 for x in rand_ce if x > ce) / N_RANDOM * 100
    pf = '**Yes**' if any(r[1] == vid for r in pareto) else ''
    lines.append(f"| {vid} | {n} | {vt} | {nc} | {et:.4f} | {delta:+.1f} | {etp:.1f} | {ethp:.1f} | {etvp:.1f} | {ce:.4f} | {cep:.1f} | {wb:.1f} | {pf} |")

with open('output/table1.md', 'w') as f:
    f.write('\n'.join(lines))
print("Table saved: output/table1.md")

# ================================================================
# GENERATE TABLE 2
# ================================================================
lines = [
    "# Table 2: Compatibility burden of biological genetic codes\n",
    "| ID | Name | Context | Dominates standard | Reassigned codons | Stop→sense | Sense→stop | Sense→sense | Stop/sense boundary | Changes |",
    "|---:|:-----|:--------|:---:|---:|---:|---:|---:|---:|:---|",
]
for n, vid, vt, et, eth, etv, ce, wb, nc in results_sorted:
    b = compatibility[vid]
    lines.append(
        f"| {vid} | {n} | {b['context']} | {'yes' if b['dominates'] else 'no'} | "
        f"{b['reassigned_codons']} | {b['stop_to_sense']} | {b['sense_to_stop']} | "
        f"{b['sense_to_sense']} | {b['stop_sense_boundary']} | {b['changes'] or 'none'} |"
    )

with open('output/table2_compatibility.md', 'w') as f:
    f.write('\n'.join(lines))
print("Table saved: output/table2_compatibility.md")

print("\nOutput checksums:")
for out_path in ("output/figure1.png", "output/figure1.pdf", "output/figure1.tiff", "output/table1.md", "output/table2_compatibility.md"):
    print(f"  {out_path}: {file_sha256(out_path)}")
print(f"Runtime seconds: {time.time() - START_TIME:.1f}")

print(f"\n{'='*60}")
print("VERIFICATION COMPLETE")
print(f"{'='*60}")
