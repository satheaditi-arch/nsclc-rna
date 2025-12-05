# ============================================================
#Differential Expression 
# ============================================================

import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


INPUT_FILE = "data/processed/adata_clustered.h5ad"
OUTPUT_DIR = "results/de_genes"
FIG_DIR = "results/figures"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


print("Loading clustered AnnData...")
adata = sc.read_h5ad(INPUT_FILE)
print("Data shape:", adata.shape)


#malignant vs non malignant 
POSSIBLE_COLUMNS = ["cell_type", "malignancy", "tumor_status", "annotation"]

label_col = None
for c in POSSIBLE_COLUMNS:
    if c in adata.obs.columns:
        label_col = c
        break

if label_col is None:
    raise ValueError(
        "No malignant/non-malignant column found.\n"
        "Expected one of: cell_type, malignancy, tumor_status, annotation\n"
        "Please add this manually before running DE."
    )

print(f"Using label column: {label_col}")

# Standardize labels
adata.obs["malignant_binary"] = (
    adata.obs[label_col]
    .astype(str)
    .str.lower()
    .str.contains("malignant|tumor|cancer")
)

print("Malignant counts:")
print(adata.obs["malignant_binary"].value_counts())
adata = adata[adata.obs["malignant_binary"].isin([True, False])].copy()

if "counts" in adata.layers:
    counts = adata.layers["counts"]
else:
    counts = adata.X.copy()

counts = pd.DataFrame(
    counts.toarray() if hasattr(counts, "toarray") else counts,
    index=adata.obs_names,
    columns=adata.var_names
)

metadata = adata.obs[["malignant_binary"]].copy()
metadata["malignant_binary"] = metadata["malignant_binary"].astype(int)

print("Running PyDESeq2...")

dds = DeseqDataSet(
    counts=counts,
    metadata=metadata,
    design_factors="malignant_binary",
    refit_cooks=True
)

dds.deseq2()

stat_res = DeseqStats(dds)
stat_res.summary()

res_df = stat_res.results_df.copy()
res_df["gene"] = res_df.index

res_df["significant"] = (
    (res_df["padj"] < 0.05) &
    (np.abs(res_df["log2FoldChange"]) > 1)
)

res_df.to_csv(f"{OUTPUT_DIR}/de_genes_malignant_vs_nonmalignant_FULL.csv")

sig_df = res_df[res_df["significant"]].sort_values("padj")
sig_df.to_csv(f"{OUTPUT_DIR}/de_genes_malignant_vs_nonmalignant_SIG.csv")

print("Top significant genes:")
print(sig_df.head(10))

plt.figure(figsize=(8, 7))

plt.scatter(
    res_df["log2FoldChange"],
    -np.log10(res_df["padj"]),
    c=res_df["significant"],
    alpha=0.6
)

plt.axvline(1, linestyle="--")
plt.axvline(-1, linestyle="--")
plt.axhline(-np.log10(0.05), linestyle="--")

plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10 Adjusted P-value")
plt.title("Malignant vs Non-Malignant Differential Expression")

plt.savefig(f"{FIG_DIR}/volcano_malignant_vs_nonmalignant.png", dpi=300)
plt.close()


top100 = sig_df.head(100)
top100.to_csv(f"{OUTPUT_DIR}/top100_biomarker_genes.csv", index=False)
