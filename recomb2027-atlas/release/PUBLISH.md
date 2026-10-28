# Published files

The paper uses the scientific source, prepared inputs and checked results at [commit 6581f628](https://github.com/arnavperam/geocert/tree/6581f628a5ac8fb83629f9e0dee0ffa4d685375d/recomb2027-atlas). This immutable revision is the reproducibility reference.

The [main paper](../manuscript/main.pdf), [supplement](../manuscript/supplement.pdf) and [abstract](../ABSTRACT.txt) are stored beside their LaTeX sources in the current checkout. Paper and documentation revisions are separate from the frozen scientific archive. [PUBLICATION.json](PUBLICATION.json) records the source identity and document hashes.

Use the [archive instructions](../ARCHIVE.md) to reconstruct the scientific package and the [README](../README.md) to run it. The archive includes its own file manifest; the checked parts reconstruct the same archive without a network request. Prepared inputs and compressed result records are restored within the extracted package.

[Validation status](../docs/VALIDATION_STATUS.md) links the completed Python 3.12/3.13 CI runs and states what was recomputed. The main and supplemental PDFs can be rebuilt by running `pdflatex main.tex` and `pdflatex supplement.tex` twice each from `manuscript/`.
