# GenbankMiner

A Python tool for mining and summarising GenBank genome files. Point it at a folder of `.gb` / `.gbk` / `.gbff` files and it produces a genome summary, a flat feature table with protein sequences, and a multi-sheet copy number report — all with a single function call or CLI command.

---

## Features

- Handles single or multi-record GenBank files (chromosomes, plasmids, contigs)
- Extracts all CDS, tRNA, rRNA, and ncRNA features into a flat CSV table including translated protein sequences
- Reports gene copy numbers by gene name, product description, and normalised product
- Breaks down copy numbers by replicon (chromosome vs plasmid)
- Counts rRNA and tRNA copies by anticodon
- Outputs a plain-text genome summary, a copy number report text file, and a multi-sheet Excel workbook

---

## Requirements

- Python 3.8+
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [openpyxl](https://openpyxl.readthedocs.io/)

Install dependencies:

```bash
pip install biopython pandas openpyxl
```

---

## Usage

### Command line

```bash
python genbankminer.py /path/to/gb/folder
```

Optionally specify a custom output directory:

```bash
python genbankminer.py /path/to/gb/folder /path/to/output
```

### Python / Jupyter notebook

```python
from genbankminer import run_genbank_miner

out = run_genbank_miner("data/")
```

With a custom output directory:

```python
out = run_genbank_miner("data/", output_dir="results/my_genome")
```

The function returns a dict so you can continue working with the results in the same session:

```python
out["feature_df"]   # pandas DataFrame of all features
out["summary"]      # genome summary dict
out["results"]      # full copy number results dict
out["dfs"]          # individual copy number DataFrames
out["output_dir"]   # path where files were written
```

---

## Input

Place your GenBank files in a single folder. The following file extensions are recognised:

```
.gb  .gbk  .gbff  .genbank
```

Multiple files in the same folder are treated as a single genome (e.g. separate files for chromosome and plasmids). Multi-record files (e.g. draft assemblies with many contigs) are handled automatically.

---

## Output

All files are written to `genbankminer_output/` inside your input folder by default, or to the path you specify.

| File | Description |
|---|---|
| `summary.txt` | Organism name, genome size, replicon list, and total feature type counts |
| `all_features.csv` | Flat table of every CDS, tRNA, rRNA, and ncRNA feature with locus tag, gene name, product, protein ID, and translated protein sequence |
| `copy_number_report.txt` | Plain-text copy number report printed during the run |
| `copy_number.xlsx` | Excel workbook with one sheet per analysis (see below) |

### Excel sheets

| Sheet | Contents |
|---|---|
| `by_gene_name` | Copy count per named gene (e.g. `dnaA`, `rpoB`) |
| `by_product` | Copy count per product description string |
| `by_normalized_product` | Copy count after normalising product strings (collapses minor naming variants) with a column showing the raw variants that were merged |
| `by_replicon` | Copy count per gene broken down by chromosome vs plasmid |
| `trna` | tRNA copy count per anticodon |
| `rrna` | rRNA copy count per product (5S, 16S, 23S) |

---

## Individual functions

The helper functions are importable independently if you want to use parts of the pipeline in your own code:

```python
from genbankminer import (
    load_genome,
    genome_summary,
    build_feature_table,
    count_by_gene_name,
    count_by_product,
    count_by_normalized_product,
    copy_count_table,
    count_rna_copies,
    full_copy_number_report,
)

records = load_genome("data/")

# flat feature table
df = build_feature_table(records)

# search for genes by keyword
hits = df[df["product"].str.contains("transporter", case=False, na=False)]

# copy counts
gene_copies = count_by_gene_name(records)
trna, rrna  = count_rna_copies(records)
```

---

## Notes

- Gene copy counting by name (`by_gene_name`) only includes features that have a `gene` qualifier — hypothetical proteins and many poorly annotated genes will not appear here. Use `by_product` or `by_normalized_product` for broader coverage.
- Plasmid detection uses a combination of topology, sequence length (< 500 kbp circular), and the presence of the word "plasmid" in the record description or keywords. This heuristic works well for standard NCBI RefSeq records but may need adjustment for non-standard annotations.
- The `translation` column in `all_features.csv` contains the amino acid sequence as deposited in the GenBank record. It is absent for tRNA, rRNA, and ncRNA rows.


---
## Citations

Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, Volume 25, Issue 11, June 2009, Pages 1422–1423, https://doi.org/10.1093/bioinformatics/btp163

Anthropic. (2025). Claude (Claude Sonnet 4.6) [Large language model]. https://www.anthropic.com

---

## License

MIT
