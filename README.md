# Genetic Code Pareto Analysis

Companion repository for:

**"Mitochondrial Genetic Codes Dominate the Standard Code on Error Tolerance and Biosynthetic Cost: Evidence for a Hidden Optimality Axis"**

Bo Brothers, Project Aletheia (2026)

## Reproduce every number in the paper

```bash
pip install -r requirements.txt
python verify_and_plot.py
```

This single script:
- Defines all 24 NCBI genetic codes from raw codon tables (no external data files)
- Defines amino acid properties (polar requirement, hydropathy, volume, biosynthetic cost)
- Scores every code on error tolerance and cost efficiency
- Generates 1,000,000 random codes (seed 42) for percentile baselines
- Identifies Pareto-dominated codes
- Computes the trade-off correlation (Kendall tau)
- Computes wobble-position synonymy for every code and the per-codon wobble/reassignment correlation
- Computes compatibility burden: reassigned codons, stop-to-sense changes, sense-to-stop changes, and sense-to-sense changes
- Prints every claim-relevant number to stdout
- Generates Figure 1 (PNG + PDF + TIFF), Table 1 (Markdown), and Table 2 compatibility burden (Markdown)

No external data files. No pre-computed results. Raw definitions in, verified numbers + figure out.

## Output

- `output/figure1.png` — Pareto frontier of all 24 biological genetic codes (300 DPI)
- `output/figure1.pdf` — Vector format
- `output/figure1.tiff` — 300 DPI submission figure
- `output/table1.md` — Complete scoring table
- `output/table2_compatibility.md` — Reassignment and stop/sense compatibility burden table

## Key finding

Eight variant codes (six mitochondrial, two nuclear) dominate the standard genetic code on **both** error tolerance and biosynthetic cost efficiency simultaneously. The standard code sits inside the Pareto frontier, not on it. Improvements on both axes are positively correlated (Kendall tau = 0.51, p = 0.0004). The dominant alternatives are concentrated in mitochondrial/organelle contexts, and most cross the stop/sense boundary, supporting the interpretation that the standard code persists because of broader compatibility constraints.

## Software

Tested with Python 3.14.3, scipy 1.17.1, numpy 2.4.3, pandas 3.0.1, and matplotlib 3.10.8. Compatible with Python 3.10+.

## License

CC-BY 4.0

## Contact

Bo Brothers
Project Aletheia
bo@projectaletheia.org
https://project-aletheia.vercel.app/
