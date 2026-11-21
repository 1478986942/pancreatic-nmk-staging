# -*- coding: utf-8 -*-
"""
GSE202051 Practice â€?Nature-style figures from human PDAC snRNA-seq
Data: 80,000 nuclei, 44 patients, 18 naive / 25 treated
Output: painting_libraries/practice_figs_gse202051/

Figure plan:
  A) Cell type composition bar (wide + tall)          â€?Fig 5a style
  B) UMAP colored by broad_celltypes                   â€?Fig 3g style  
  C) UMAP colored by response within immune subset      â€?immune-focused
  D) Violin + points: SERPIN fraction by response       â€?bar+dot style
  E) Top marker genes per cell type (dotplot-like)      â€?Fig 5b style
  F) Heatmap of top 10 markers per major type            â€?Fig 3c style
"""
import os, warnings, random
warnings.filterwarnings('ignore')
random.seed(42)

import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
import seaborn as sns

sns.set_style("whitegrid")
rcParams['font.family'] = 'Arial'
rcParams['axes.linewidth'] = 0.4
rcParams['lines.linewidth'] = 0.5

ROOT  = r"E:\pancreatic_cancer\submit"
ADATA = os.path.join(ROOT, "gse202051_hum", "processed", "GSE202051_qc_filtered.h5ad")
OUT   = os.path.join(ROOT, "painting_libraries", "practice_figs_gse202051")
os.makedirs(OUT, exist_ok=True)

# â”€â”€ Load â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
print("Loading adata...")
adata = ad.read_h5ad(ADATA)
print(f"  Shape: {adata.shape}")
print(f"  broad_celltypes: {dict(adata.obs['broad_celltypes'].value_counts())}")
print(f"  response:        {dict(adata.obs['response'].value_counts())}")
print(f"  status:          {dict(adata.obs['status'].value_counts())}")

# â”€â”€ Palette (Nature style â€?colorblind friendly) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
PAL_CT    = {"Epithelial":   "#d62728",
             "Endothelial":  "#1f77b4",
             "Immune":       "#2ca02c",
             "Fibroblast":   "#ff7f0e",
             "Endocrine":    "#9467bd",
             "Schwann":      "#8c564b",
             "unknown":      "#bcbd22"}
PAL_RESP  = {"Moderate response": "#2ca02c",
             "Poor response":     "#d62728",
             "Minimal response":  "#ff7f0e",
             "Untreated":         "#7f7f7f"}

def nature_theme():
    fig = plt.gcf()
    for ax in fig.axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.35)
        ax.spines['bottom'].set_linewidth(0.35)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontname('Arial')
            label.set_fontsize(9)
        ax.xaxis.grid(True, color="#ededed", linewidth=0.5)
        ax.yaxis.grid(True, color="#ededed", linewidth=0.5)
        ax.set_axisbelow(True)

# â”€â”€ Helper: annotate panel letter â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
def add_letter(ax, letter, y_offset=0.0):
    ax.text(-0.08, 1.0 + y_offset, letter, transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='bottom', ha='left',
            fontname='Arial')

# =====================================================================
# A. Cell type composition â€?two complementary stacked bars
#    A1: within each cell type, naive vs treated proportion
#    A2: within each status, cell type composition (sums to 100%)
# =====================================================================
print("\n[Fig A] Cell type composition ...")
ct_df = adata.obs.groupby(['broad_celltypes', 'status']).size().reset_index(name='count')

order = ['Epithelial','Immune','Fibroblast','Endothelial','Endocrine','Schwann']
PAL_STATUS = {"naive": "#4d4d4d", "treated": "#1f77b4"}

# --- Pivot 1: rows = celltype, cols = status, normalised per row ---
piv_row = ct_df.pivot(index='broad_celltypes', columns='status', values='count').fillna(0)
piv_row = piv_row.loc[[o for o in order if o in piv_row.index]]
piv_row = piv_row.div(piv_row.sum(axis=1), axis=0) * 100

# --- Pivot 2: rows = celltype, cols = status, normalised per column ---
piv_col = ct_df.pivot(index='broad_celltypes', columns='status', values='count').fillna(0)
piv_col = piv_col.loc[[o for o in order if o in piv_col.index]]
piv_col = piv_col.div(piv_col.sum(axis=0), axis=1) * 100

fig_a, axes = plt.subplots(1, 2, figsize=(8, 3.2), gridspec_kw={'width_ratios': [1.25, 1]})
ax_a1, ax_a2 = axes

# --- A1: horizontal, within celltype status split ---
status_cols = [c for c in ["naive", "treated"] if c in piv_row.columns]
left = np.zeros(len(piv_row))
for st in status_cols:
    vals = piv_row[st].values
    ax_a1.barh(piv_row.index, vals, left=left,
               color=PAL_STATUS[st], edgecolor='white', linewidth=0.6,
               label=st)
    left += vals
ax_a1.set_xlabel('% of cells (naive vs treated)', fontsize=9)
ax_a1.set_xlim(0, 100)
ax_a1.set_title('A', fontsize=13, fontweight='bold', loc='left', pad=6)
ax_a1.legend(title=None, loc='lower center', ncol=2, fontsize=7,
             frameon=False, bbox_to_anchor=(0.5, -0.42), handlelength=1.2)
nature_theme()
ax_a1.yaxis.grid(False)

# --- A2: vertical, within status celltype composition ---
ct_order_col = [o for o in order if o in piv_col.index]
bottom = np.zeros(len(piv_col.columns))
for ct in ct_order_col:
    vals = piv_col.loc[ct].values
    ax_a2.bar(piv_col.columns, vals, bottom=bottom,
              color=PAL_CT[ct], edgecolor='white', linewidth=0.6,
              label=ct)
    bottom += vals
ax_a2.set_ylabel('% of cells', fontsize=9)
ax_a2.set_xlabel('Treatment status', fontsize=9)
ax_a2.set_ylim(0, 100)
ax_a2.set_title('B', fontsize=13, fontweight='bold', loc='left', pad=6)
ax_a2.legend(title=None, loc='center left', bbox_to_anchor=(1.02, 0.5),
             fontsize=7, frameon=False, handlelength=1.0)
nature_theme()
ax_a2.xaxis.grid(False)

plt.tight_layout(pad=1.2)
plt.savefig(os.path.join(OUT, "gse202051_A_celltype_composition.png"), dpi=250, bbox_inches='tight')
plt.close()
print("  Saved A")


# =====================================================================
# B. UMAP coloured by broad_celltypes (Fig 3g style, sampled)
# =====================================================================
print("[Fig B] UMAP broad_celltypes ...")
X_umap = adata.obsm['X_umap']
n_sample_B = min(25000, adata.n_obs)
rng = np.random.default_rng(42)
idx_B = rng.choice(adata.n_obs, n_sample_B, replace=False)
sub_B = adata[idx_B].copy()

fig_b, ax_b = plt.subplots(figsize=(5.8, 4.8))
for ct in sub_B.obs['broad_celltypes'].cat.categories:
    mask = sub_B.obs['broad_celltypes'] == ct
    ax_b.scatter(sub_B.obsm['X_umap'][mask, 0], sub_B.obsm['X_umap'][mask, 1],
                 s=1.5, alpha=0.6, c=PAL_CT[ct], label=ct, edgecolor='none')
ax_b.set_xlabel('UMAP 1', fontsize=9)
ax_b.set_ylabel('UMAP 2', fontsize=9)
ax_b.set_title('B', fontsize=13, fontweight='bold', loc='left', pad=4)
ax_b.legend(loc='upper right', fontsize=7, framealpha=0.9, ncol=2, 
            handletextpad=0.3, borderaxespad=0.3)
nature_theme()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "gse202051_B_umap_broadcelltype.png"), dpi=250, bbox_inches='tight')
plt.close()
print("  Saved B")


# =====================================================================
# C. UMAP immune cells only, colored by response  
# =====================================================================
print("[Fig C] UMAP immune x response ...")
immune_mask = adata.obs['broad_celltypes'] == 'Immune'
imm_adata = adata[immune_mask].copy()
resp_order = ['Moderate response', 'Minimal response', 'Poor response']
imm_adata.obs['response_clean'] = imm_adata.obs['response'].cat.remove_unused_categories()

n_sample_C = min(15000, len(imm_adata))
idx_C = rng.choice(len(imm_adata), n_sample_C, replace=False)
sub_C = imm_adata[idx_C]

fig_c, ax_c = plt.subplots(figsize=(5.8, 4.8))
for resp in resp_order:
    mask = sub_C.obs['response'] == resp
    ax_c.scatter(sub_C.obsm['X_umap'][mask, 0], sub_C.obsm['X_umap'][mask, 1],
                 s=1.2, alpha=0.5, c=PAL_RESP[resp], label=resp, edgecolor='none')
ax_c.set_xlabel('UMAP 1', fontsize=9)
ax_c.set_ylabel('UMAP 2', fontsize=9)
ax_c.set_title('C', fontsize=13, fontweight='bold', loc='left', pad=4)
ax_c.legend(loc='upper right', fontsize=7, framealpha=0.9)
nature_theme()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "gse202051_C_umap_immune_response.png"), dpi=250, bbox_inches='tight')
plt.close()
print("  Saved C")


# =====================================================================
# D. SERPIN fraction per patient â€?box+violon+jitter (Fig 2 style)
# =====================================================================
print("[Fig D] SERPIN fraction by response ...")
# Compute per-patient SERPIN fraction from SPP1/MARCO expression
malignant_mask = adata.obs['broad_celltypes'] == 'Epithelial'
mal_adata = adata[malignant_mask].copy()
patient_resp = mal_adata.obs[['pid', 'response']].drop_duplicates()
# Get patient-level SERPINE1+SPP1+ fraction
gene_expr = mal_adata.X.toarray() if hasattr(mal_adata.X, 'toarray') else mal_adata.X
serpin_genes = ['SERPINE1', 'SERPINB2', 'SPP1', 'MARCO']
gene_idx = [i for i, g in enumerate(adata.var_names) if g in serpin_genes]
if len(gene_idx) < len(serpin_genes):
    print(f"  WARNING: not all SERPIN genes found, using {len(gene_idx)}")
gene_means = gene_expr[:, gene_idx].mean(axis=1)  # per-cell mean SERPIN signal

mal_adata.obs['serpin_mean'] = gene_means
pat_stats = mal_adata.obs.groupby('pid')['serpin_mean'].agg(['mean','count'])
pat_resp  = mal_adata.obs[['pid','response']].drop_duplicates().set_index('pid')['response']
pat_df    = pat_stats.join(pat_resp).rename(columns={'mean':'serpin_frac'}).dropna()
# Only analyse treated patients with response labels
pat_df = pat_df[pat_df['response'].isin(['Moderate response','Minimal response','Poor response'])]
# Also keep only patients with enough cells
pat_df = pat_df[pat_df['count'] >= 15]

fig_d, ax_d = plt.subplots(figsize=(4.8, 4.0))
responses = pat_df['response'].cat.remove_unused_categories()
groups = []
values = []
for resp in responses:
    sub = pat_df[pat_df['response']==resp]
    groups.extend([resp]*len(sub))
    values.extend(sub['serpin_frac'].tolist())
df_plot = pd.DataFrame({'Response': groups, 'Serpin Fraction': values})

# Violin + jitter
sns.violinplot(x='Response', y='Serpin Fraction', data=df_plot, 
               inner=None, palette=PAL_RESP, scale='count', ax=ax_d, linewidth=0.8)
sns.stripplot(x='Response', y='Serpin Fraction', data=df_plot, 
              color='white', edgecolor='black', linewidth=0.3,
              size=4, alpha=0.7, ax=ax_d, jitter=True)
ax_d.set_xlabel('Response category', fontsize=9)
ax_d.set_ylabel('Mean SERPIN signal (malignant)', fontsize=9)
ax_d.set_title('D', fontsize=13, fontweight='bold', loc='left', pad=4)
ax_d.tick_params(axis='x', rotation=15, labelsize=8)
nature_theme()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "gse202051_D_serpin_fraction_by_response.png"), dpi=250, bbox_inches='tight')
plt.close()
print("  Saved D")


# =====================================================================
# E. Dotplot: top markers per cell type (enrichplot dotplot style)
# =====================================================================
print("[Fig E] Marker dotplot ...")
# The source adata.X is already log-normalised (verified: values are non-integer,
# log1p scale). Do NOT re-normalise â€?just rank genes on the existing matrix.
counts = adata.layers['counts']
print("  Using source log-normalised X for ranking")

markers_out = []
try:
    # Use scanpy's built-in ranking (Wilcoxon)
    sc.tl.rank_genes_groups(adata, groupby='broad_celltypes', method='wilcoxon')
    for ct in adata.obs['broad_celltypes'].cat.categories:
        grp = adata.uns['rank_genes_groups']['names'][ct][:5]
        pvals = adata.uns['rank_genes_groups']['pvals_adj'][ct][:5]
        logfc = adata.uns['rank_genes_groups']['logfoldchanges'][ct][:5]
        for g, lfc, p in zip(grp, logfc, pvals):
            if pd.notna(g) and pd.notna(lfc):
                markers_out.append({
                    'celltype': ct, 'gene': g,
                    'log2FC': float(lfc) if np.isfinite(float(lfc)) else 0,
                    'p_adj': max(float(p), 1e-300) if pd.notna(p) else 1e-300})
    
    df_markers = pd.DataFrame(markers_out)
    print(f"  Markers collected: {len(df_markers)}")
    
    fig_e, ax_e = plt.subplots(figsize=(6.5, 4.2))
    keep_ct = [ct for ct in ['Epithelial','Fibroblast','Immune','Endothelial']
               if ct in df_markers['celltype'].unique()]

    # --- True marker dotplot: rows = cell types, cols = top genes ---
    # fraction expressing from counts layer (sparse-safe), mean from log-normalised X
    var_list = list(adata.var_names)

    # Top 4 genes per cell type (keep order)
    top_genes = []
    for ct in keep_ct:
        gs = adata.uns['rank_genes_groups']['names'][ct][:4]
        for g in gs:
            if g in adata.var_names:
                gidx = var_list.index(g)
                top_genes.append((ct, g, gidx))
    seen = set(); gene_order = []
    for _, g, gi in top_genes:
        if g not in seen:
            seen.add(g); gene_order.append((g, gi))

    gene_labels = [g for g, _ in gene_order]
    gene_idxs = [gi for _, gi in gene_order]
    ct_order_plot = keep_ct

    frac_mat = np.zeros((len(ct_order_plot), len(gene_idxs)))
    mean_mat = np.zeros((len(ct_order_plot), len(gene_idxs)))
    for r, ct in enumerate(ct_order_plot):
        mask = (adata.obs['broad_celltypes'] == ct).values
        for c, gi in enumerate(gene_idxs):
            cnt_col = counts[mask, gi]
            cnt_col = cnt_col.toarray().flatten() if hasattr(cnt_col, 'toarray') else np.asarray(cnt_col).flatten()
            frac_mat[r, c] = (cnt_col > 0).mean() * 100
            expr_col = adata.X[mask, gi]
            expr_col = expr_col.toarray().flatten() if hasattr(expr_col, 'toarray') else np.asarray(expr_col).flatten()
            mean_mat[r, c] = expr_col.mean()
    # Scale mean per gene (z-score across cell types)
    mu = mean_mat.mean(axis=0, keepdims=True)
    sd = mean_mat.std(axis=0, keepdims=True) + 1e-9
    z_mat = (mean_mat - mu) / sd

    xs, ys, ss, cc = [], [], [], []
    for r in range(len(ct_order_plot)):
        for c in range(len(gene_idxs)):
            xs.append(c); ys.append(r)
            ss.append(max(frac_mat[r, c], 1) * 4 + 3)
            cc.append(z_mat[r, c])
    scatter = ax_e.scatter(xs, ys, s=ss, c=cc, cmap='RdBu_r',
                           vmin=-2, vmax=2, edgecolor='black', linewidth=0.3)
    ax_e.set_xticks(range(len(gene_labels)))
    ax_e.set_xticklabels(gene_labels, rotation=45, ha='right', fontsize=7, style='italic')
    ax_e.set_yticks(range(len(ct_order_plot)))
    ax_e.set_yticklabels(ct_order_plot, fontsize=8)
    for tick_label, ct in zip(ax_e.get_yticklabels(), ct_order_plot):
        tick_label.set_color(PAL_CT.get(ct, '#333'))
        tick_label.set_fontweight('bold')
    ax_e.set_xlabel('Top marker genes', fontsize=9)
    ax_e.set_ylabel('Cell type', fontsize=9)
    ax_e.set_title('E', fontsize=13, fontweight='bold', loc='left', pad=4)
    cbar = plt.colorbar(scatter, ax=ax_e, shrink=0.7, pad=0.02)
    cbar.set_label('Scaled mean expr', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    # Size legend
    for frac in [25, 50, 75]:
        ax_e.scatter([], [], s=frac*4+3, c='lightgrey', edgecolor='black',
                     linewidth=0.3, label=f'{frac}%')
    ax_e.legend(title='% expressing', loc='upper left', bbox_to_anchor=(1.12, 1.0),
                fontsize=6, title_fontsize=7, frameon=False,
                labelspacing=1.2, borderpad=0.8)
    nature_theme()
    ax_e.grid(False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "gse202051_E_marker_dotplot.png"), dpi=250, bbox_inches='tight')
    plt.close()
    print("  Saved E")
except Exception as ex:
    import traceback; traceback.print_exc()
    print(f"  Fig E failed: {ex}")


# =====================================================================
# F. Heatmap: top 8 markers Ã— cell types (ComplexHeatmap style)
# =====================================================================
print("[Fig F] Marker heatmap ...")
if markers_out:
    df_all = pd.DataFrame(markers_out)
    df_all['celltype'] = df_all['celltype'].astype('category')
    df_all = df_all.rename(columns={'pvals_adj':'p_adj','log1pFC':'log2FC'})
    # Top 8 per celltype
    top8_per_ct = []
    for ct in df_all['celltype'].cat.categories:
        top8 = df_all[df_all['celltype']==ct].nlargest(8, 'log2FC')
        top8_per_ct.append(top8)
    df_top8 = pd.concat(top8_per_ct, ignore_index=True)
    
    # Build expression matrix: rows=top8 markers, cols=cell types
    all_genes = df_top8['gene'].unique()
    celltypes_used = sorted(df_top8['celltype'].unique())
    
    mat_data = np.zeros((len(all_genes), len(celltypes_used)))
    gene_to_row = {g: i for i, g in enumerate(all_genes)}
    ct_to_col = {ct: i for i, ct in enumerate(celltypes_used)}
    
    for _, row in df_top8.iterrows():
        gene, ct = row['gene'], row['celltype']
        idx = adata.var_names.get_loc(gene)
        ct_mask = adata.obs['broad_celltypes'] == ct
        vals = np.asarray(adata.X[ct_mask, idx].toarray()).flatten().astype(float)
        mat_data[gene_to_row[gene], ct_to_col[ct]] = vals.mean() if len(vals)>0 else 0
    
    mat = pd.DataFrame(mat_data, index=all_genes, columns=celltypes_used)
    # Z-score across genes for each cell type
    mat_z = (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-9)
    
    fig_f, ax_f = plt.subplots(figsize=(5.5, 4.5))
    sns.heatmap(mat_z, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                cbar_kws={'shrink':0.7, 'label':'Z-score'},
                linewidths=0.3, linecolor='white', ax=ax_f, square=False)
    ax_f.set_xlabel('Cell type', fontsize=9)
    ax_f.set_ylabel('Marker gene', fontsize=9)
    ax_f.set_title('F', fontsize=13, fontweight='bold', loc='left', pad=4)
    ax_f.tick_params(axis='x', rotation=25, labelsize=8)
    ax_f.tick_params(axis='y', labelsize=8)
    nature_theme()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "gse202051_F_marker_heatmap.png"), dpi=250, bbox_inches='tight')
    plt.close()
    print("  Saved F")

print(f"\nAll practice figures saved to: {OUT}")
print("Files:", sorted(os.listdir(OUT)))
