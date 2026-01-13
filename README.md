# Single-Cell RNA-Seq Analysis to Identify Potential Biomarkers of Non-Small Cell Lung Cancer (NSCLC)

**Authors:** Aditi S, Jenna K, Liv R, Jennifer G, Montserrat P  
**Goal:** Build a reproducible discovery + validation pipeline to identify NSCLC-associated gene signatures and pathways using **single-cell RNA-seq** (training) and **bulk RNA-seq / microarray** cohorts (validation), with clinical relevance assessed via **TCGA survival**.

---

## Project Overview

Non-small cell lung cancer (NSCLC) comprises **>85%** of lung cancers and is highly heterogeneous, making robust biomarker discovery challenging. This project replicates and extends a comprehensive bioinformatics workflow (inspired by Sultana et al., 2023) using **Python**, emphasizing:

- Tumor vs non-malignant cell-state differences at single-cell resolution
- Pseudobulk differential expression to reduce scRNA noise
- Pathway-level reproducibility across datasets
- Clinical survival validation using TCGA LUAD/LUSC
- A lightweight diagnostic model via ROC/AUC evaluation

**Key training genes (single-cell):** `ARHGAP9`, `C16orf54`, `CHST11`, `AOAH`  
**Most reproducible bulk marker:** `C16orf54` (consistent tumor-suppressed signal)

**Consistent pathway themes across datasets:**  
- **p53 / DNA damage response**
- **Hypoxia**
- **Reactive Oxygen Species (ROS)**
- immune signaling programs (context-dependent)

---

## Data Sources

### Training (Single-cell)
- **Single-Cell Lung Cancer Atlas (CELLxGENE collection)**  
  https://cellxgene.cziscience.com/collections/edb893ee-4066-4128-9aec-5eb2b03f8287

### Validation (Bulk)
- **GSE40419** (paired tumor–normal LUAD; bulk RNA-seq)  
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40419
- **GSE19188** (91 tumor, 65 adjacent normal lung samples; microarray/bulk expression)  
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188
- **GSE30219** (307 NSCLC tumors only; Affymetrix; tumor-only stratification + GSEA)  
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30219

### Normal baseline reference
- **Tabula Sapiens Lung** (normal lung baseline expression via CELLxGENE)  
  https://tabula-sapiens.sf.czbiohub.org/about

### Clinical survival cohorts
- **TCGA LUAD** and **TCGA LUSC** (GDC portal)  
  https://portal.gdc.cancer.gov/projects/TCGA-LUAD  
  https://portal.gdc.cancer.gov/projects/TCGA-LUSC  

---

## Methods Summary (Pipeline)

1. **QC & Preprocessing (Scanpy)**
   - Subsample ~30k cells for compute efficiency
   - QC filters (counts / genes detected / mito%)
   - Normalize to 10k counts + log1p
   - Select HVGs (Seurat v3 method)

2. **Clustering & Visualization**
   - PCA → kNN graph → UMAP
   - Leiden clustering
   - Cluster marker detection (Wilcoxon)

3. **Pseudobulk Differential Expression**
   - Aggregate cells by sample/donor and condition (malignant vs non-malignant)
   - DE testing on pseudobulk profiles (Wilcoxon + FDR)
   - Export full DE tables + volcano plots + top genes

4. **Trajectory Inference**
   - PAGA on Leiden clusters (topology/connectivity)
   - Diffusion maps + diffusion pseudotime (DPT) rooted in a major non-malignant cluster

5. **Pathway Enrichment**
   - Enrichr via GSEApy on top biomarker genes (GO BP, KEGG)
   - Also run relaxed thresholds when gene lists are small

6. **Gene Set Enrichment Analysis (GSEA)**
   - Preranked GSEA on pseudobulk DE stats
   - MSigDB Hallmark + KEGG
   - Summarize enriched biological programs via NES and FDR

7. **Regulatory / Co-expression Network**
   - Pearson correlation across malignant biomarker genes
   - Build network edges for |r| ≥ threshold
   - Identify hubs (degree centrality), visualize modules

8. **ROC Curve Analysis**
   - Logistic regression classifier on biomarker gene expression
   - Train/test split (stratified), standardize features
   - Evaluate accuracy/precision/recall/AUC + plot ROC

9. **Baseline Expression in Normal Lung**
   - Use Tabula Sapiens Lung to deprioritize genes broadly expressed in healthy lung

10. **TCGA Survival Analysis**
   - Kaplan–Meier (log-rank) for high vs low expression
   - Cox proportional hazards model with clinical covariates where available

11. **External Validation in Independent Bulk Cohorts**
   - Confirm tumor vs normal directionality in GSE40419 & GSE19188
   - Tumor-only stratification in GSE30219 (high vs low C16orf54) + Hallmark GSEA

---

## Key Findings

### Training (single-cell)
- Clear malignant vs non-malignant transcriptional differences
- **C16orf54** highlighted as a strong candidate marker (among top gene candidates)

### Bulk Validation (reproducibility)
- In **GSE40419** and **GSE19188**:
  - All four candidates trend lower in tumors
  - **C16orf54** shows **strong, consistent tumor-suppressed** behavior and best reproducibility

### Tumor-only validation (GSE30219)
- **Low C16orf54** tumors enrich proliferative programs:
  - **E2F Targets**, **G2M Checkpoint**, **MYC Targets**
- **High C16orf54** tumors enrich immune/inflammatory programs:
  - Interferon responses, TNF-α/NF-κB, IL-6/JAK/STAT3, complement

## Citation

If you use or reference this work, please cite the original study that inspired and informed this pipeline:

Sultana A, Alam MS, Liu X, Sharma R, Singla RK, Gundamaraju R, Shen B.  
**Single-cell RNA-seq analysis to identify potential biomarkers for diagnosis and prognosis of non-small cell lung cancer by using comprehensive bioinformatics approaches.**  
*Translational Oncology.* 2023 Jan;27:101571.  
https://doi.org/10.1016/j.tranon.2022.101571

This project replicates and extends the above study using a fully Python-based workflow (Scanpy, GSEApy, lifelines, scikit-learn), incorporating additional validation datasets and pathway-level analyses.
## How This Project Differs from Sultana et al. (2023)

This project was inspired by the study of Sultana et al. (2023) but introduces several important methodological, technical, and analytical extensions that distinguish it as an independent and expanded investigation.

**1. Python-Based Reimplementation**  
While Sultana et al. conducted their analysis primarily in R using Seurat and Monocle, this project fully reimplements the pipeline in **Python**, leveraging Scanpy, GSEApy, lifelines, NetworkX, and scikit-learn. This improves accessibility for Python-based workflows and enables seamless integration with modern machine learning tools.

**2. Pseudobulk Differential Expression Framework**  
Rather than relying solely on single-cell–level differential expression, this project adopts a **pseudobulk strategy** that aggregates cells by sample and condition. This reduces single-cell noise, incorporates biological replication, and improves the robustness and reproducibility of differential expression results.

**3. Explicit Training–Validation Separation**  
This study explicitly separates **training** (Single-Cell Lung Cancer Atlas) and **validation** (GSE40419, GSE19188) datasets. Additional orthogonal validation is performed using a **tumor-only cohort (GSE30219)**, enabling intra-tumoral stratification independent of normal tissue comparisons—an extension not emphasized in the original study.

**4. Baseline Normal Lung Filtering Using Tabula Sapiens**  
Candidate biomarkers are filtered using **Tabula Sapiens Lung** single-cell reference data to deprioritize genes with high baseline expression in healthy lung tissue. This step increases cancer specificity and was not a core component of the original pipeline.

**5. Integrated Pathway-Level Reproducibility Analysis**  
Given the heterogeneity of NSCLC, this project emphasizes **pathway-level consistency** across datasets rather than relying solely on gene-level overlap. Hallmark and KEGG gene set enrichment analyses are used systematically to identify conserved biological programs such as hypoxia, p53 signaling, ROS response, and immune activation.

**6. Regulatory Network and Systems-Level Interpretation**  
Beyond differential expression, this project constructs **gene–gene co-expression networks** to identify hub genes and modular regulatory programs. The resulting networks support a systems-level model of NSCLC progression rather than a single-driver gene hypothesis.

**7. Quantitative Diagnostic Evaluation via ROC Analysis**  
A supervised classification step using **logistic regression and ROC/AUC evaluation** is incorporated to quantify the diagnostic potential of identified biomarkers. This provides a translational performance metric not explicitly reported in Sultana et al.

**8. Extended Clinical Validation and Interpretation**  
Clinical relevance is assessed using both **Kaplan–Meier survival analysis and Cox proportional hazards modeling** in TCGA LUAD and LUSC cohorts, with additional interpretation of subtype-specific effects and pathway associations.

Together, these extensions transform the original descriptive bioinformatics framework into a **scalable, reproducible discovery and validation engine** suitable for biomarker prioritization, translational research, and future integration with machine learning–driven diagnostic or therapeutic pipelines.

