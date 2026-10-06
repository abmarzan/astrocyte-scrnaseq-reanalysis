# %% [markdown]
# # Astrocyte single-cell RNA-seq: a training re-analysis of GSE114000 with Scanpy
#
# **What this is.** A self-directed training exercise in single-cell RNA-seq analysis.
# It re-analyses the public Smart-seq2 dataset of adult mouse cortex and hippocampus
# astrocytes from Batiuk, Martirosyan et al., *Nature Communications* 2020
# (GEO accession GSE114000) and checks how far a standard Scanpy workflow recovers
# the cell types and astrocyte subtypes (AST1 to AST5) that the authors published.
#
# **What this is not.** It is not new biology and not a replication of the paper's
# own pipeline. The authors' labels are used only at the end, as a reference to
# compare against.
#
# Steps: load counts and metadata, quality control, normalisation, highly variable
# genes, PCA, neighbour graph, UMAP, Leiden clustering, marker genes, cell-type
# annotation, astrocyte sub-clustering, comparison with the published labels.

# %%
import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

SEED = 0
np.random.seed(SEED)
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, dpi_save=200, frameon=False, figsize=(4.5, 4))

ROOT = Path(".")
DATA = ROOT / "data"
FIG = ROOT / "figures"
RES = ROOT / "results"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

COUNTS = DATA / "GSE114000_Counts_Batiuk_Martirosyan_Supplementary_Data3.tsv.gz"
META = DATA / "GSE114000_Metadata_Batiuk_Martirosyan_Supplementary_Data1.xlsx"
print("scanpy", sc.__version__, "| anndata", ad.__version__)

# %% [markdown]
# ## 1. Load the count matrix and the published metadata
#
# The count file is genes by cells (HTSeq counts, one Smart-seq2 library per cell).
# Scanpy expects cells by genes, so the table is transposed. The metadata sheet
# gives the brain area, plate and the authors' cell-type label for every cell.

# %%
counts = pd.read_csv(COUNTS, sep="\t", skiprows=2, index_col=0)
meta = pd.read_excel(META, header=1)
meta["Cell ID"] = meta["Cell ID"].astype(str).str.replace(r"^X", "", regex=True)
meta = meta.set_index("Cell ID")
print("counts (genes x cells):", counts.shape)
print("metadata rows:", meta.shape[0])
print("cells in both:", len(set(counts.columns) & set(meta.index)))

adata = ad.AnnData(X=sparse.csr_matrix(counts.T.values.astype(np.float32)))
adata.obs_names = counts.columns.astype(str)
adata.var_names = counts.index.astype(str)
adata.var_names_make_unique()
adata.obs["region"] = meta.reindex(adata.obs_names)["Brain area"].astype(str).values
adata.obs["plate"] = meta.reindex(adata.obs_names)["Plate name"].astype(str).values
adata.obs["sample"] = meta.reindex(adata.obs_names)["Sample name"].astype(str).values
adata.obs["published_type"] = meta.reindex(adata.obs_names)["Type"].astype(str).values
adata

# %%
print(adata.obs["region"].value_counts().to_string())
print()
print(adata.obs["published_type"].value_counts().to_string())

# %% [markdown]
# ## 2. Quality control
#
# Three gene sets are flagged: mitochondrial genes (`mt-`), ERCC spike-ins and
# ribosomal protein genes. For each cell we compute the number of detected genes,
# the total counts and the share of counts from mitochondrial genes and spike-ins.
# A high mitochondrial or spike-in share usually means a damaged or empty well.

# %%
adata.var["mt"] = adata.var_names.str.startswith("mt-")
adata.var["ercc"] = adata.var_names.str.startswith("ERCC-")
adata.var["ribo"] = adata.var_names.str.match(r"^Rp[sl]\d")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ercc", "ribo"], percent_top=None, log1p=True, inplace=True)
qc_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ercc"]
print(adata.obs[qc_cols].describe().round(2).to_string())

# %%
fig, axes = plt.subplots(1, 4, figsize=(15, 3.6))
titles = ["Genes detected per cell", "Total counts per cell", "Mitochondrial counts (%)", "ERCC spike-in counts (%)"]
for ax, col, title in zip(axes, qc_cols, titles):
    parts = [adata.obs.loc[adata.obs["region"] == r, col].values for r in ["Cortex", "Hippocampus"]]
    ax.violinplot(parts, showmedians=True)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Cortex", "Hippocampus"])
    ax.set_title(title, fontsize=10)
    if col == "total_counts":
        ax.set_yscale("log")
    ax.grid(False)
fig.tight_layout()
fig.savefig(FIG / "01_qc_violin.png", dpi=200, bbox_inches="tight")
plt.show()

# %% [markdown]
# ### Filtering
#
# The authors already removed low-quality libraries before depositing the data,
# so the thresholds here are deliberately mild and data-driven: a cell is removed
# when it is an outlier by more than 5 median absolute deviations (MAD) on the log
# number of genes or log total counts, or when its mitochondrial share is more
# than 5 MADs above the median. Genes seen in fewer than 5 cells are dropped, and
# the spike-ins are removed after QC because they are not biological signal.

# %%
def mad_outlier(values, nmads, side="both"):
    med = np.median(values)
    mad = np.median(np.abs(values - med))
    lower = values < med - nmads * mad
    upper = values > med + nmads * mad
    return {"both": lower | upper, "upper": upper, "lower": lower}[side]

adata.obs["outlier_depth"] = mad_outlier(adata.obs["log1p_total_counts"].values, 5) | mad_outlier(
    adata.obs["log1p_n_genes_by_counts"].values, 5
)
adata.obs["outlier_mt"] = mad_outlier(adata.obs["pct_counts_mt"].values, 5, side="upper")
n_before = adata.n_obs
print("depth/complexity outliers:", int(adata.obs["outlier_depth"].sum()))
print("mitochondrial outliers:", int(adata.obs["outlier_mt"].sum()))
adata = adata[~(adata.obs["outlier_depth"] | adata.obs["outlier_mt"])].copy()
adata = adata[:, ~adata.var["ercc"]].copy()
sc.pp.filter_genes(adata, min_cells=5)
print(f"cells kept: {adata.n_obs} of {n_before}; genes kept: {adata.n_vars}")

# %% [markdown]
# ## 3. Normalisation and feature selection
#
# Counts are scaled to the median library size and log-transformed. The 2,000
# most variable genes are then selected for dimensionality reduction. The raw
# counts stay in a separate layer.

# %%
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
print("highly variable genes:", int(adata.var["highly_variable"].sum()))

# %% [markdown]
# ## 4. PCA, neighbour graph, UMAP and Leiden clustering

# %%
sc.tl.pca(adata, n_comps=50, mask_var="highly_variable", random_state=SEED)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
sc.tl.leiden(adata, resolution=0.6, key_added="leiden", flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
print(adata.obs["leiden"].value_counts().sort_index().to_string())

# %%
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
sc.pl.umap(adata, color="leiden", ax=axes[0], show=False, title="Leiden clusters (this analysis)", legend_loc="on data", legend_fontsize=8)
sc.pl.umap(adata, color="region", ax=axes[1], show=False, title="Brain region")
sc.pl.umap(adata, color="published_type", ax=axes[2], show=False, title="Published cell type (Batiuk et al. 2020)")
fig.tight_layout()
fig.savefig(FIG / "02_umap_all_cells.png", dpi=200, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 5. Which cell type is each cluster?
#
# Clusters are named from canonical marker genes, without looking at the
# published labels. Each cluster gets the cell type whose markers have the
# highest mean scaled expression in it.

# %%
MARKERS = {
    "Astrocyte": ["Aldoc", "Aqp4", "Gja1", "Slc1a3", "Slc1a2", "Aldh1l1"],
    "Neuron": ["Snap25", "Syt1", "Rbfox3", "Stmn2"],
    "Oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mobp"],
    "OPC": ["Pdgfra", "Cspg4", "Olig1"],
    "Microglia/macrophage": ["Cx3cr1", "C1qa", "P2ry12", "Csf1r"],
    "Endothelial": ["Cldn5", "Flt1", "Pecam1"],
    "Mural": ["Pdgfrb", "Rgs5", "Kcnj8", "Acta2"],
    "Ependymal": ["Foxj1", "Ccdc153", "Tmem212"],
}
MARKERS = {k: [g for g in v if g in adata.var_names] for k, v in MARKERS.items()}

for name, genes in MARKERS.items():
    sc.tl.score_genes(adata, genes, score_name=f"score_{name}", random_state=SEED)
score_cols = [f"score_{k}" for k in MARKERS]
cluster_scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
cluster_scores.columns = list(MARKERS)
cluster_call = cluster_scores.idxmax(axis=1)
adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_call).astype("category")
print(cluster_scores.round(2).assign(call=cluster_call).to_string())

# %%
sc.pl.dotplot(adata, MARKERS, groupby="leiden", standard_scale="var", show=False)
plt.savefig(FIG / "03_marker_dotplot.png", dpi=200, bbox_inches="tight")
plt.show()

# %% [markdown]
# ### Data-driven marker genes per cluster
#
# A Wilcoxon rank-sum test of each cluster against all other cells, as an
# independent check on the names given above.

# %%
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", pts=True)
top = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(10)
top.columns = [f"{c} ({cluster_call[c]})" for c in top.columns]
top.to_csv(RES / "top10_markers_per_cluster.csv", index=False)
print(top.T.to_string(header=False))

# %% [markdown]
# ### Agreement with the published broad cell types
#
# The authors label non-astrocyte cells by type and astrocytes as AST1 to AST5.
# For this comparison all AST labels are collapsed to "Astrocyte".

# %%
def broad(label):
    return "Astrocyte" if str(label).startswith("AST") else str(label)

adata.obs["published_broad"] = adata.obs["published_type"].map(broad)
broad_table = pd.crosstab(adata.obs["cell_type"], adata.obs["published_broad"])
broad_table.to_csv(RES / "celltype_vs_published_broad.csv")
print(broad_table.to_string())
is_astro_mine = adata.obs["cell_type"] == "Astrocyte"
is_astro_pub = adata.obs["published_broad"] == "Astrocyte"
agree = float((is_astro_mine == is_astro_pub).mean())
print(f"\nAstrocyte versus non-astrocyte agreement with the published labels: {agree:.1%}")

# %% [markdown]
# ## 6. Astrocyte subtypes
#
# The cells called astrocytes here are taken forward, highly variable genes are
# re-selected within them and they are clustered again. The published subtype
# labels (AST1 to AST5) are then used as the reference.

# %%
astro = adata[adata.obs["cell_type"] == "Astrocyte"].copy()
astro.X = astro.layers["counts"].copy()
astro.uns.pop("log1p", None)  # X holds raw counts again, so clear the log1p flag
sc.pp.filter_genes(astro, min_cells=5)
sc.pp.normalize_total(astro)
sc.pp.log1p(astro)
sc.pp.highly_variable_genes(astro, n_top_genes=2000, flavor="seurat")
sc.tl.pca(astro, n_comps=30, mask_var="highly_variable", random_state=SEED)
sc.pp.neighbors(astro, n_neighbors=15, n_pcs=15, random_state=SEED)
sc.tl.umap(astro, random_state=SEED)
print("astrocytes analysed:", astro.n_obs)

rows = []
has_pub = astro.obs["published_type"].str.startswith("AST").values
for res in [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]:
    key = f"leiden_{res}"
    sc.tl.leiden(astro, resolution=res, key_added=key, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
    rows.append(
        {
            "resolution": res,
            "clusters": astro.obs[key].nunique(),
            "ARI": adjusted_rand_score(astro.obs["published_type"][has_pub], astro.obs[key][has_pub]),
            "NMI": normalized_mutual_info_score(astro.obs["published_type"][has_pub], astro.obs[key][has_pub]),
        }
    )
sweep = pd.DataFrame(rows)
sweep.to_csv(RES / "astrocyte_resolution_sweep.csv", index=False)
print(sweep.round(3).to_string(index=False))

# %% [markdown]
# The resolution is not tuned to the published answer. It is fixed at the value
# that gives five clusters (the number of subtypes the paper reports); if several
# values give five, the lowest is used. The sweep above shows how much the
# agreement depends on that choice.

# %%
five = sweep[sweep["clusters"] == 5]
RES_CHOSEN = float(five["resolution"].min()) if len(five) else float(sweep.iloc[(sweep["clusters"] - 5).abs().argsort().iloc[0]]["resolution"])
astro.obs["subcluster"] = astro.obs[f"leiden_{RES_CHOSEN}"].astype("category")
ari = adjusted_rand_score(astro.obs["published_type"][has_pub], astro.obs["subcluster"][has_pub])
nmi = normalized_mutual_info_score(astro.obs["published_type"][has_pub], astro.obs["subcluster"][has_pub])
print(f"resolution used: {RES_CHOSEN} -> {astro.obs['subcluster'].nunique()} clusters; ARI = {ari:.3f}; NMI = {nmi:.3f}")

# %%
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
sc.pl.umap(astro, color="subcluster", ax=axes[0], show=False, title="Astrocyte sub-clusters (this analysis)", legend_loc="on data")
sc.pl.umap(astro, color="published_type", ax=axes[1], show=False, title="Published subtype (AST1 to AST5)")
sc.pl.umap(astro, color="region", ax=axes[2], show=False, title="Brain region")
fig.tight_layout()
fig.savefig(FIG / "04_umap_astrocytes.png", dpi=200, bbox_inches="tight")
plt.show()

# %%
conf = pd.crosstab(astro.obs["subcluster"], astro.obs["published_type"])
conf = conf[[c for c in conf.columns if c.startswith("AST")] + [c for c in conf.columns if not c.startswith("AST")]]
conf.to_csv(RES / "astrocyte_subclusters_vs_published.csv")
frac = conf.div(conf.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(1.1 * conf.shape[1] + 2, 0.6 * conf.shape[0] + 1.6))
im = ax.imshow(frac.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(conf.shape[1]))
ax.set_xticklabels(conf.columns, rotation=45, ha="right")
ax.set_yticks(range(conf.shape[0]))
ax.set_yticklabels([f"cluster {i}" for i in conf.index])
for i in range(conf.shape[0]):
    for j in range(conf.shape[1]):
        v = int(conf.values[i, j])
        if v:
            ax.text(j, i, str(v), ha="center", va="center", fontsize=8, color="white" if frac.values[i, j] > 0.55 else "black")
ax.set_xlabel("Published label")
ax.set_ylabel("This analysis")
ax.set_title(f"Cells per sub-cluster and published subtype (ARI = {ari:.2f})", fontsize=10)
ax.grid(False)
fig.colorbar(im, ax=ax, label="Share of the sub-cluster")
fig.tight_layout()
fig.savefig(FIG / "05_subcluster_vs_published.png", dpi=200, bbox_inches="tight")
plt.show()
print(conf.to_string())

# %% [markdown]
# ### Region composition and marker genes of the sub-clusters
#
# The paper's central finding is that astrocyte subtypes differ between cortex
# and hippocampus. The composition of each sub-cluster by region shows whether
# this simple re-analysis sees the same thing.

# %%
comp = pd.crosstab(astro.obs["subcluster"], astro.obs["region"])
comp_frac = comp.div(comp.sum(axis=1), axis=0)
comp.to_csv(RES / "astrocyte_subclusters_by_region.csv")
ax = comp_frac.plot(kind="barh", stacked=True, figsize=(6, 0.5 * comp.shape[0] + 1.4), color=["#4C78A8", "#F58518"], width=0.8)
ax.set_xlabel("Share of cells")
ax.set_ylabel("Astrocyte sub-cluster")
ax.set_xlim(0, 1)
ax.legend(title="Region", bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
ax.grid(False)
plt.tight_layout()
plt.savefig(FIG / "06_subcluster_region_composition.png", dpi=200, bbox_inches="tight")
plt.show()
print(comp.assign(hippocampus_share=comp_frac.get("Hippocampus", 0).round(2)).to_string())

# %%
sc.tl.rank_genes_groups(astro, "subcluster", method="wilcoxon", pts=True)
sub_top = pd.DataFrame(astro.uns["rank_genes_groups"]["names"]).head(15)
sub_top.columns = [f"cluster {c}" for c in sub_top.columns]
sub_top.to_csv(RES / "astrocyte_subcluster_top15_markers.csv", index=False)
print(sub_top.head(10).T.to_string(header=False))
sc.tl.dendrogram(astro, "subcluster")
sc.pl.rank_genes_groups_dotplot(astro, n_genes=5, standard_scale="var", show=False)
plt.savefig(FIG / "07_subcluster_marker_dotplot.png", dpi=200, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 7. A technical confound worth checking: dissociation-induced genes
#
# Enzymatic tissue dissociation switches on immediate-early genes (Fos, Jun,
# Egr1 and others) in some cells. That signal is a preparation artefact, not
# resting biology, and it can form its own cluster. Scoring it per cell shows
# where it sits. In work on stress-related disorders this matters twice over,
# because the same genes also respond to real stress in vivo.

# %%
IEG = [g for g in ["Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Ier2", "Ier3", "Dusp1", "Atf3", "Hspa1a", "Hspa1b", "Nr4a1", "Klf2", "Klf4", "Btg2"] if g in adata.var_names]
sc.tl.score_genes(adata, IEG, score_name="ieg_score", random_state=SEED)
ieg_by_cluster = adata.obs.groupby("leiden", observed=True)["ieg_score"].agg(["mean", "median", "size"]).round(2)
ieg_by_cluster["cell_type"] = cluster_call
ieg_by_cluster.to_csv(RES / "ieg_score_by_cluster.csv")
print("genes in the score:", ", ".join(IEG))
print(ieg_by_cluster.to_string())
astro.obs["ieg_score"] = adata.obs["ieg_score"].reindex(astro.obs_names).values
print()
print(astro.obs.groupby("subcluster", observed=True)["ieg_score"].agg(["mean", "size"]).round(2).to_string())

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
sc.pl.umap(adata, color="ieg_score", ax=axes[0], show=False, title="Immediate-early gene score, all cells", cmap="viridis")
sc.pl.umap(astro, color="ieg_score", ax=axes[1], show=False, title="Immediate-early gene score, astrocytes", cmap="viridis")
fig.tight_layout()
fig.savefig(FIG / "08_dissociation_ieg_score.png", dpi=200, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 8. Summary of what the numbers say

# %%
summary = {
    "cells_in_deposited_matrix": int(n_before),
    "cells_after_qc": int(adata.n_obs),
    "genes_after_filtering": int(adata.n_vars),
    "leiden_clusters_all_cells": int(adata.obs["leiden"].nunique()),
    "cells_called_astrocyte": int(astro.n_obs),
    "astrocyte_vs_other_agreement_with_published": round(agree, 4),
    "astrocyte_subcluster_resolution": RES_CHOSEN,
    "astrocyte_subclusters": int(astro.obs["subcluster"].nunique()),
    "ARI_subclusters_vs_published_AST": round(float(ari), 4),
    "NMI_subclusters_vs_published_AST": round(float(nmi), 4),
    "ieg_high_cluster": str(ieg_by_cluster["mean"].idxmax()),
    "ieg_high_cluster_cells": int(ieg_by_cluster.loc[ieg_by_cluster["mean"].idxmax(), "size"]),
}
pd.Series(summary).to_csv(RES / "summary.csv", header=["value"])
for k, v in summary.items():
    print(f"{k}: {v}")
adata.obs.to_csv(RES / "cell_annotations_all_cells.csv")
astro.obs[["region", "plate", "published_type", "subcluster"]].to_csv(RES / "cell_annotations_astrocytes.csv")
