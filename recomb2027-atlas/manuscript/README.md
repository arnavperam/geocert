# Manuscript

`main.tex` builds the paper; `supplement.tex` builds the appendix. The paper has a separate title page, 10 body pages and two reference pages.

```sh
cd recomb2027-atlas/manuscript
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error supplement.tex
pdflatex -interaction=nonstopmode -halt-on-error supplement.tex
```

The geometry schematic is generated with `python scripts/geometry_pipeline.py` and uses illustrative coordinates. To regenerate the reach plot from the restored scientific results:

```sh
python scripts/reach_figures.py --data-root /path/to/extracted/recomb2027-atlas
```

These presentation scripts use NumPy and Matplotlib. The exact algorithms, inputs and recorded results remain in the pinned scientific source identified by `../release/PUBLICATION.json`. The abstract in `../ABSTRACT.txt` matches the current paper.
