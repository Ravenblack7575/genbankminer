"""
genbankminer.py

Run from command line:
    python genbankminer.py /path/to/gb/folder

Or import and call directly:
    from genbankminer import run_genbank_miner
    run_genbank_miner("data/", output_dir="results/")
"""

from Bio import SeqIO
from collections import Counter
import os
import re
import pandas as pd


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_genome(folder_path):
    """Load all GenBank files from a folder. Returns a list of SeqRecords."""
    paths = [
        os.path.join(folder_path, f) for f in os.listdir(folder_path)
        if f.endswith((".gb", ".gbk", ".gbff", ".genbank"))
    ]
    if not paths:
        raise FileNotFoundError(f"No GenBank files found in: {folder_path}")

    records = []
    for path in sorted(paths):
        recs = list(SeqIO.parse(path, "genbank"))
        print(f"  {os.path.basename(path)}: {len(recs)} record(s)")
        records.extend(recs)

    print(f"  --> {len(records)} total record(s) loaded\n")
    return records


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def genome_summary(records):
    """Return a summary dict describing the genome at a high level."""
    organism  = records[0].annotations.get("organism", "unknown")
    total_bp  = sum(len(r.seq) for r in records)
    all_feats = [f for r in records for f in r.features]
    feat_counts = Counter(f.type for f in all_feats)

    replicons = []
    for r in records:
        topology  = r.annotations.get("topology", "?")
        is_plasmid = (
            topology == "circular" and len(r.seq) < 500_000
            or "plasmid" in r.description.lower()
            or "plasmid" in r.annotations.get("keywords", [""])[0].lower()
        )
        replicons.append({
            "id":         r.id,
            "length_bp":  len(r.seq),
            "type":       "plasmid" if is_plasmid else "chromosome/contig",
            "topology":   topology,
            "features":   Counter(f.type for f in r.features),
        })

    return {
        "organism":       organism,
        "num_records":    len(records),
        "total_bp":       total_bp,
        "replicons":      replicons,
        "total_features": feat_counts,
    }


def format_summary(summary):
    """Return the genome summary as a formatted string."""
    lines = []
    lines.append(f"Organism   : {summary['organism']}")
    lines.append(f"Records    : {summary['num_records']}")
    lines.append(f"Total size : {summary['total_bp']:,} bp")
    lines.append("")
    lines.append(f"{'Replicon':<25} {'Length (bp)':>14}  {'Type':<22} {'Topology'}")
    lines.append("-" * 72)
    for r in summary["replicons"]:
        lines.append(
            f"  {r['id']:<23} {r['length_bp']:>14,}  {r['type']:<22} {r['topology']}"
        )
    lines.append("")
    lines.append(f"{'Feature type':<25} {'Count':>8}")
    lines.append("-" * 35)
    for ftype, count in sorted(summary["total_features"].items()):
        lines.append(f"  {ftype:<23} {count:>8}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Feature table
# ---------------------------------------------------------------------------

def build_feature_table(records):
    """Build a flat DataFrame of all CDS/tRNA/rRNA/ncRNA features."""
    rows = []
    for rec in records:
        organism = rec.annotations.get("organism", "unknown")
        source   = rec.annotations.get("source",   "unknown")
        for f in rec.features:
            if f.type in ("CDS", "tRNA", "rRNA", "ncRNA"):
                rows.append({
                    "organism":    organism,
                    "source":      source,
                    "record_id":   rec.id,
                    "type":        f.type,
                    "locus_tag":   f.qualifiers.get("locus_tag",   ["?"])[0],
                    "gene":        f.qualifiers.get("gene",         [""])[0],
                    "product":     f.qualifiers.get("product",      ["hypothetical protein"])[0],
                    "protein_id":  f.qualifiers.get("protein_id",   [""])[0],
                    "translation": f.qualifiers.get("translation",  [""])[0],
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Copy number functions
# ---------------------------------------------------------------------------

def count_by_gene_name(records, feature_type="CDS"):
    names = [
        f.qualifiers["gene"][0]
        for rec in records
        for f in rec.features
        if f.type == feature_type and "gene" in f.qualifiers
    ]
    return Counter(names)


def count_by_product(records, feature_type="CDS"):
    products = [
        f.qualifiers.get("product", ["hypothetical protein"])[0]
        for rec in records
        for f in rec.features
        if f.type == feature_type
    ]
    return Counter(products)


def normalize_product(s):
    s = s.lower().strip()
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r",?\s*(subunit|domain|protein|chain|component)\s*\w*$", "", s)
    s = re.sub(r"\d+$", "", s).strip()
    return s


def count_by_normalized_product(records):
    norm_counts   = Counter()
    norm_examples = {}
    for rec in records:
        for f in rec.features:
            if f.type == "CDS":
                raw  = f.qualifiers.get("product", ["hypothetical protein"])[0]
                norm = normalize_product(raw)
                norm_counts[norm] += 1
                norm_examples.setdefault(norm, set()).add(raw)
    return norm_counts, norm_examples


def copy_count_table(records):
    rows = []
    for rec in records:
        is_plasmid = (
            "plasmid" in rec.description.lower()
            or "plasmid" in rec.annotations.get("keywords", [""])[0].lower()
        )
        replicon = "plasmid" if is_plasmid else "chromosome"
        for f in rec.features:
            if f.type == "CDS":
                rows.append({
                    "record_id": rec.id,
                    "replicon":  replicon,
                    "gene":      f.qualifiers.get("gene",      [""])[0],
                    "product":   f.qualifiers.get("product",   ["hypothetical protein"])[0],
                    "locus_tag": f.qualifiers.get("locus_tag", ["?"])[0],
                })
    df = pd.DataFrame(rows)
    gene_counts = (
        df[df["gene"] != ""]
        .groupby(["gene", "replicon"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    gene_counts["total"] = gene_counts.drop(columns="gene").sum(axis=1)
    gene_counts = gene_counts.sort_values("total", ascending=False)
    return df, gene_counts


def count_rna_copies(records):
    trna = Counter()
    rrna = Counter()
    for rec in records:
        for f in rec.features:
            if f.type == "tRNA":
                trna[f.qualifiers.get("product", ["unknown"])[0]] += 1
            elif f.type == "rRNA":
                rrna[f.qualifiers.get("product", ["unknown"])[0]] += 1
    return trna, rrna


# ---------------------------------------------------------------------------
# Copy number report
# ---------------------------------------------------------------------------

def full_copy_number_report(records, top_n=20):
    """Run all copy number analyses and return a results dict."""
    by_name                    = count_by_gene_name(records)
    by_product                 = count_by_product(records)
    norm_counts, norm_examples = count_by_normalized_product(records)
    cct_df, gene_counts        = copy_count_table(records)
    trna, rrna                 = count_rna_copies(records)

    lines = []

    lines.append("=== BY GENE NAME ===")
    dups = {g: n for g, n in by_name.items() if n > 1}
    lines.append(f"  {len(by_name)} named genes, {len(dups)} with >1 copy")
    for g, n in sorted(dups.items(), key=lambda x: -x[1])[:top_n]:
        lines.append(f"  {n:4}x  {g}")

    lines.append("\n=== BY PRODUCT (top repeated) ===")
    for p, n in by_product.most_common(top_n):
        if n > 1:
            lines.append(f"  {n:4}x  {p}")

    lines.append("\n=== BY NORMALIZED PRODUCT (top repeated) ===")
    for norm, n in norm_counts.most_common(top_n):
        if n > 1:
            variants = norm_examples[norm]
            variant_str = f"  [{', '.join(list(variants)[:2])}]" if len(variants) > 1 else ""
            lines.append(f"  {n:4}x  {norm}{variant_str}")

    lines.append("\n=== COPIES BY REPLICON ===")
    multi_copy = gene_counts[gene_counts["total"] > 1]
    if not multi_copy.empty:
        lines.append(multi_copy.to_string(index=False))
    else:
        lines.append("  No multi-copy genes found by replicon")

    lines.append("\n=== rRNA COPIES ===")
    for product, n in sorted(rrna.items()):
        lines.append(f"  {n:3}x  {product}")

    lines.append(f"\n=== tRNA COPIES ({sum(trna.values())} total, {len(trna)} anticodons) ===")
    for product, n in sorted(trna.items()):
        lines.append(f"  {n:3}x  {product}")

    report_text = "\n".join(lines)
    print(report_text)

    return {
        "report_text":   report_text,
        "by_name":       by_name,
        "by_product":    by_product,
        "norm_counts":   norm_counts,
        "norm_examples": norm_examples,
        "cct_df":        cct_df,
        "gene_counts":   gene_counts,
        "trna":          trna,
        "rrna":          rrna,
    }


# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

def results_to_dataframes(results):
    """Convert the results dict into a dict of named DataFrames."""
    dfs = {}

    dfs["by_gene_name"] = pd.DataFrame(
        results["by_name"].items(), columns=["gene", "copy_count"]
    ).sort_values("copy_count", ascending=False).reset_index(drop=True)

    dfs["by_product"] = pd.DataFrame(
        results["by_product"].items(), columns=["product", "copy_count"]
    ).sort_values("copy_count", ascending=False).reset_index(drop=True)

    norm_rows = [
        {
            "normalized_product": norm,
            "copy_count":         count,
            "raw_variants":       "; ".join(sorted(results["norm_examples"][norm])),
        }
        for norm, count in results["norm_counts"].items()
    ]
    dfs["by_normalized_product"] = (
        pd.DataFrame(norm_rows)
        .sort_values("copy_count", ascending=False)
        .reset_index(drop=True)
    )

    dfs["by_replicon"] = results["gene_counts"].reset_index(drop=True)

    dfs["trna"] = pd.DataFrame(
        results["trna"].items(), columns=["product", "copy_count"]
    ).sort_values("product").reset_index(drop=True)

    dfs["rrna"] = pd.DataFrame(
        results["rrna"].items(), columns=["product", "copy_count"]
    ).sort_values("product").reset_index(drop=True)

    return dfs


def save_to_excel(dfs, path):
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in dfs.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f"  sheet '{sheet_name}'  ({len(df)} rows)")
    print(f"  --> saved {path}")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_genbank_miner(folder_path, output_dir=None):
    """
    Given a folder of GenBank files, produce:
      - summary.txt       : genome summary text
      - all_features.csv  : flat table of all CDS/RNA features with sequences
      - copy_number.xlsx  : copy number report across multiple sheets

    Parameters
    ----------
    folder_path : str
        Path to folder containing .gb / .gbk / .gbff files.
    output_dir : str, optional
        Where to write output files. Defaults to folder_path/genbankminer_output/
    """
    folder_path = os.path.abspath(folder_path)

    if output_dir is None:
        output_dir = os.path.join(folder_path, "genbankminer_output")
    os.makedirs(output_dir, exist_ok=True)

    sep = "-" * 50

    # --- Load ---
    print(sep)
    print("Loading GenBank files")
    print(sep)
    records = load_genome(folder_path)

    # --- Summary ---
    print(sep)
    print("Genome Summary")
    print(sep)
    summary     = genome_summary(records)
    summary_txt = format_summary(summary)
    print(summary_txt)

    summary_path = os.path.join(output_dir, "summary.txt")
    with open(summary_path, "w") as fh:
        fh.write(summary_txt + "\n")
    print(f"\n  --> saved {summary_path}")

    # --- Feature table ---
    print(f"\n{sep}")
    print("Building feature table")
    print(sep)
    feature_df   = build_feature_table(records)
    features_path = os.path.join(output_dir, "all_features.csv")
    feature_df.to_csv(features_path, index=False)
    print(f"  {len(feature_df)} features extracted")
    print(f"  --> saved {features_path}")

    # --- Copy number report ---
    print(f"\n{sep}")
    print("Copy Number Report")
    print(sep)
    results = full_copy_number_report(records)

    # save copy number report text alongside summary
    report_path = os.path.join(output_dir, "copy_number_report.txt")
    with open(report_path, "w") as fh:
        fh.write(results["report_text"] + "\n")
    print(f"\n  --> saved {report_path}")

    # save excel
    dfs          = results_to_dataframes(results)
    excel_path   = os.path.join(output_dir, "copy_number.xlsx")
    print("")
    save_to_excel(dfs, excel_path)

    print(f"\n{sep}")
    print(f"Done. All outputs written to: {output_dir}")
    print(sep)

    return {
        "records":     records,
        "summary":     summary,
        "feature_df":  feature_df,
        "results":     results,
        "dfs":         dfs,
        "output_dir":  output_dir,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python genbankminer.py <folder_path> [output_dir]")
        sys.exit(1)

    folder   = sys.argv[1]
    out_dir  = sys.argv[2] if len(sys.argv) > 2 else None
    run_genbank_miner(folder, output_dir=out_dir)
