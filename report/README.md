# report

Paper draft, imported from the Overleaf starter project (Springer Nature `sn-jnl` template, Bioinformatics section layout).

- `main.tex`: preamble, authors, abstract, back matter; pulls in each section with `\input`
- `ref.bib`: references
- `sec/`: one file per section (introduction, methods, algorithm, implementation, discussion), so edits to different sections don't conflict
- `sn-jnl.cls`, `sn-nature.bst`: template files, don't edit

## build

On Perlmutter:

```
module load texlive/2024
cd report && latexmk -pdf main.tex
```

Locally: any TeX Live with `latexmk`. Clean up with `latexmk -C`.

Build outputs (including `main.pdf`) are gitignored, so pull and build rather than committing PDFs.
