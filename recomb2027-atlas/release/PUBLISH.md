# Public evaluation candidate

The canonical candidate is the `submission/atlas-1.1.0` branch of `arnavperam/geocert`, in `recomb2027-atlas/`. Its scientific source, prepared inputs and checked results are published together. The manuscript cites the immutable scientific source commit; a later document commit contains the rebuilt PDFs. No stable tag or conference submission is implied by the version label.

The browsable code matches the checked scientific archive. `ARCHIVE.md` in the public source folder explains how to reconstruct that archive and verify each part. Restore large records with `python release/restore_results.py`, then run the commands in the README. The read-only CI checks the exact candidate under Python 3.12 and 3.13. Consult the actual workflow result before claiming public CI success.

The final evaluation package contains the same scientific code and result records plus the corrected manuscript, supplement, aligned abstract and all LaTeX inputs. Build both PDFs with two pdflatex passes from `manuscript/`; the main body remains ten pages. `docs/FINAL_EVALUATION.txt` and the final publication record describe changes from Agent 1's original archive. The original archive identity is retained for provenance.

A specialist cold read and the author's understanding of the proofs remain part of evaluation. This work does not fill EasyChair, submit a paper or certify anatomical accuracy. The author-reviewed stable tag and conference submission are separate finalization steps.
