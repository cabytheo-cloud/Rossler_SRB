# Non-rigorous numerics (not part of the proof)

Double-precision exploration used to locate the tangencies and to draw Figures 1–2 of the paper.

- `num57.cpp` — scan and refinement of tangency parameters near c = 5.7 (critical point landing on p after n iterates).
- `pmap.py`, `fp.py`, `tang1.py` — Poincaré map, fixed point and tangency function in Python/SciPy.
- `figdata57.py`, `makefigs57.py` (and `figdata.py`, `makefigs.py` for the band-merging case) — figure data and plots.

These scripts were run from the original working tree; input/output paths (e.g. `near57/params_n16a.json`,
`paper/figs/`) may need to be adjusted (the parameter files are in `capd/`). Nothing in these scripts enters
the logic of the computer-assisted proofs.
