"""
STREAMLINED NETWORK DIFFUSION ANALYSIS
======================================

This script contains ONLY the essential code needed to generate the publication
Sankey plots from your notebooks. All exploratory and alternative methods have
been removed.

Input Requirements:
-------------------
1. Proteome/phosphoproteome data with log2FC values
2. GCC network from igraph (g_gcc)
3. OmniPath edges with activation/inhibition annotations

Output:
-------
1. Annotated g_gcc with all diffusion metrics
2. final_driver_table (comprehensive summary)
3. Top driver/sink tables (top_early, top_mid, top_late, top_sinks)
4. Sankey plots for t=10, 600, 1800
"""

import numpy as np
import pandas as pd
import igraph as ig
from scipy import sparse
from scipy.sparse.linalg import inv

# ============================================================================
# PART 1: LOAD DATA
# ============================================================================

def load_and_prepare_data():
    """
    Load node data and GCC network
    Attach phosphosite metrics and metadata to g_gcc
    """
    
    # Load your data (adjust paths as needed)
    # proteome_df = pd.read_csv("proteome_data.csv")
    # phospho_df = pd.read_csv("phosphoproteome_data.csv")
    # g_gcc = ig.Graph.Read_Pickle("g_gcc.pkl")
    
    # Ensure g_gcc has the following node attributes:
    # - name (UniProt ID)
    # - gene_symbol
    # - mean_logFC_10, mean_logFC_600, mean_logFC_1800
    
    return g_gcc


# ============================================================================
# PART 2: PERSONALIZED PAGERANK DIFFUSION
# ============================================================================

def compute_ppr_diffusion(g_gcc, timepoint, alpha=0.85):
    """
    Compute Personalized PageRank diffusion for a given timepoint
    
    Parameters:
    -----------
    g_gcc : igraph.Graph
        Network graph
    timepoint : str
        "10", "600", or "1800"
    alpha : float
        Restart probability (default: 0.85)
    
    Returns:
    --------
    ppr_scores : np.array
        PPR diffusion scores for each node
    """
    
    logfc_col = f"mean_logFC_{timepoint}"
    
    # Get adjacency matrix (directed)
    adj = np.array(g_gcc.get_adjacency(attribute=None).data)
    
    # Row-normalize
    row_sums = adj.sum(axis=1)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    A = adj / row_sums[:, np.newaxis]
    
    # Personalization vector (weighted by |log2FC|)
    logfc_values = np.array([v.get(logfc_col, 0) for v in g_gcc.vs])
    p = np.abs(logfc_values)
    p = p / p.sum() if p.sum() > 0 else np.ones(len(p)) / len(p)
    
    # Compute PPR: α(I - (1-α)A^T)^(-1) * p
    n = len(g_gcc.vs)
    I = np.eye(n)
    
    # Use sparse matrices for efficiency
    A_sparse = sparse.csr_matrix(A.T)
    M = I - (1 - alpha) * A_sparse.toarray()
    M_inv = np.linalg.inv(M)
    
    ppr_scores = alpha * M_inv @ p
    
    return ppr_scores


def run_diffusion_all_timepoints(g_gcc):
    """
    Run PPR diffusion for all timepoints and attach to graph
    """
    
    print("Running PPR diffusion...")
    
    for tp in ["10", "600", "1800"]:
        print(f"  Computing PPR for t={tp}...")
        ppr = compute_ppr_diffusion(g_gcc, tp)
        
        # Attach to graph
        g_gcc.vs[f"ppr_{tp}"] = ppr
    
    print("✓ Diffusion complete")
    return g_gcc


# ============================================================================
# PART 3: PAIRWISE INFLUENCE MATRICES
# ============================================================================

def compute_pairwise_influence(g_gcc, timepoint, alpha=0.85):
    """
    Compute pairwise influence matrix: who influences whom
    
    Returns:
    --------
    influence_matrix : np.array
        influence[i,j] = influence of node i on node j
    """
    
    logfc_col = f"mean_logFC_{timepoint}"
    
    # Get adjacency matrix
    adj = np.array(g_gcc.get_adjacency(attribute=None).data)
    row_sums = adj.sum(axis=1)
    row_sums[row_sums == 0] = 1
    A = adj / row_sums[:, np.newaxis]
    
    # PPR transition matrix
    n = len(g_gcc.vs)
    I = np.eye(n)
    A_sparse = sparse.csr_matrix(A.T)
    M = I - (1 - alpha) * A_sparse.toarray()
    M_inv = np.linalg.inv(M)
    ppr_matrix = alpha * M_inv
    
    # Effective influence: PPR[i,j] * |log2FC[i]|
    logfc_values = np.array([abs(v.get(logfc_col, 0)) for v in g_gcc.vs])
    influence_matrix = ppr_matrix * logfc_values[:, np.newaxis]
    
    return influence_matrix


# ============================================================================
# PART 4: POLARITY CALCULATION
# ============================================================================

def compute_polarity(g_gcc, timepoint):
    """
    Compute polarity: (outgoing - incoming) / (outgoing + incoming)
    
    Polarity > 0: driver (net outgoing influence)
    Polarity < 0: sink (net incoming influence)
    """
    
    influence_matrix = compute_pairwise_influence(g_gcc, timepoint)
    
    # Outgoing = sum of row (influence TO others)
    outgoing = influence_matrix.sum(axis=1)
    
    # Incoming = sum of column (influence FROM others)
    incoming = influence_matrix.sum(axis=0)
    
    # Polarity
    total = outgoing + incoming
    total[total == 0] = 1  # Avoid division by zero
    polarity = (outgoing - incoming) / total
    
    return polarity, influence_matrix


def compute_polarity_all_timepoints(g_gcc):
    """
    Compute polarity for all timepoints and attach to graph
    """
    
    print("Computing polarity...")
    
    influence_matrices = {}
    
    for tp in ["10", "600", "1800"]:
        print(f"  Computing polarity for t={tp}...")
        polarity, influence_matrix = compute_polarity(g_gcc, tp)
        
        # Attach to graph
        g_gcc.vs[f"polarity_{tp}"] = polarity
        influence_matrices[tp] = influence_matrix
    
    print("✓ Polarity complete")
    return g_gcc, influence_matrices


# ============================================================================
# PART 5: DRIVER EXPLANATION METRICS
# ============================================================================

def compute_driver_metrics(g_gcc, timepoint, influence_matrix, de_threshold=1.5):
    """
    Compute driver explanation metrics for each node
    
    Metrics:
    --------
    - N_explained: number of DE nodes influenced
    - explained_strength: sum of influence × |log2FC| for DE targets
    - fraction_explained: proportion of total network change explained
    """
    
    logfc_col = f"mean_logFC_{timepoint}"
    
    # Get log2FC values
    logfc_values = np.array([v.get(logfc_col, 0) for v in g_gcc.vs])
    
    # DE mask (|log2FC| > threshold)
    de_mask = np.abs(logfc_values) > de_threshold
    
    n_nodes = len(g_gcc.vs)
    
    n_explained = []
    explained_strength = []
    
    for i in range(n_nodes):
        # Influence to DE nodes
        influence_to_de = influence_matrix[i, de_mask]
        
        # N_explained: count of DE nodes influenced
        n_explained.append(np.sum(influence_to_de > 0))
        
        # explained_strength: sum of (influence × |log2FC|) for DE targets
        de_logfc = np.abs(logfc_values[de_mask])
        explained_strength.append(np.sum(influence_to_de * de_logfc))
    
    # fraction_explained
    total_change = np.sum(explained_strength)
    fraction_explained = np.array(explained_strength) / total_change if total_change > 0 else np.zeros(n_nodes)
    
    return {
        "N_explained": np.array(n_explained),
        "explained_strength": np.array(explained_strength),
        "fraction_explained": fraction_explained
    }


def compute_driver_metrics_all_timepoints(g_gcc, influence_matrices):
    """
    Compute driver metrics for all timepoints and attach to graph
    """
    
    print("Computing driver explanation metrics...")
    
    for tp in ["10", "600", "1800"]:
        print(f"  Computing metrics for t={tp}...")
        
        metrics = compute_driver_metrics(g_gcc, tp, influence_matrices[tp])
        
        # Attach to graph
        g_gcc.vs[f"N_explained_{tp}"] = metrics["N_explained"]
        g_gcc.vs[f"explained_strength_{tp}"] = metrics["explained_strength"]
        g_gcc.vs[f"fraction_explained_{tp}"] = metrics["fraction_explained"]
    
    print("✓ Driver metrics complete")
    return g_gcc


# ============================================================================
# PART 6: DRIVER CLASSIFICATION & SCORING
# ============================================================================

def classify_and_score_drivers(g_gcc):
    """
    Classify nodes as early/mid/late drivers and compute scores
    
    Driver Score = explained_strength + (N_explained × polarity)
    """
    
    print("Classifying drivers and computing scores...")
    
    for tp, score_name in [("10", "early_score"), 
                           ("600", "mid_score"), 
                           ("1800", "late_score")]:
        
        explained_strength = np.array(g_gcc.vs[f"explained_strength_{tp}"])
        n_explained = np.array(g_gcc.vs[f"N_explained_{tp}"])
        polarity = np.array(g_gcc.vs[f"polarity_{tp}"])
        
        # Compute score
        score = explained_strength + (n_explained * polarity)
        
        # Attach to graph
        g_gcc.vs[score_name] = score
    
    print("✓ Driver classification complete")
    return g_gcc


# ============================================================================
# PART 7: BUILD FINAL DRIVER TABLE
# ============================================================================

def build_final_driver_table(g_gcc):
    """
    Build comprehensive driver summary table with all metrics
    """
    
    print("Building final driver table...")
    
    data = []
    
    for v in g_gcc.vs:
        row = {
            "UniProt": v["name"],
            "gene_symbol": v.get("gene_symbol", ""),
        }
        
        # Add all timepoint-specific metrics
        for tp in ["10", "600", "1800"]:
            row[f"mean_logFC_{tp}"] = v.get(f"mean_logFC_{tp}", np.nan)
            row[f"ppr_{tp}"] = v.get(f"ppr_{tp}", np.nan)
            row[f"N_explained_{tp}"] = v.get(f"N_explained_{tp}", np.nan)
            row[f"explained_strength_{tp}"] = v.get(f"explained_strength_{tp}", np.nan)
            row[f"polarity_{tp}"] = v.get(f"polarity_{tp}"], np.nan)
        
        # Add scores
        row["early_score"] = v.get("early_score", np.nan)
        row["mid_score"] = v.get("mid_score", np.nan)
        row["late_score"] = v.get("late_score", np.nan)
        
        data.append(row)
    
    final_table = pd.DataFrame(data)
    
    print(f"✓ Final table created: {len(final_table)} rows, {len(final_table.columns)} columns")
    
    return final_table


# ============================================================================
# PART 8: SELECT TOP DRIVERS & SINKS
# ============================================================================

def select_top_drivers_and_sinks(final_table, n_top=30):
    """
    Select top drivers for each timepoint and top sinks
    """
    
    print(f"Selecting top {n_top} drivers/sinks...")
    
    # Top early drivers
    top_early = final_table.nlargest(n_top, "early_score")[
        ["UniProt", "gene_symbol", "early_score", "mean_logFC_10", "ppr_10", "polarity_10"]
    ].copy()
    
    # Top mid drivers
    top_mid = final_table.nlargest(n_top, "mid_score")[
        ["UniProt", "gene_symbol", "mid_score", "mean_logFC_600", "ppr_600", "polarity_600"]
    ].copy()
    
    # Top late drivers
    top_late = final_table.nlargest(n_top, "late_score")[
        ["UniProt", "gene_symbol", "late_score", "mean_logFC_1800", "ppr_1800", "polarity_1800"]
    ].copy()
    
    # Top sinks (lowest polarity across timepoints)
    final_table["min_polarity"] = final_table[["polarity_10", "polarity_600", "polarity_1800"]].min(axis=1)
    top_sinks = final_table.nsmallest(n_top, "min_polarity")[
        ["UniProt", "gene_symbol", "min_polarity", "polarity_10", "polarity_600", "polarity_1800"]
    ].copy()
    
    print(f"✓ Top drivers/sinks selected")
    print(f"  Early drivers: {len(top_early)}")
    print(f"  Mid drivers: {len(top_mid)}")
    print(f"  Late drivers: {len(top_late)}")
    print(f"  Sinks: {len(top_sinks)}")
    
    return top_early, top_mid, top_late, top_sinks


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def run_complete_analysis(g_gcc):
    """
    Run the complete network diffusion analysis pipeline
    
    Returns all data needed for publication Sankey plots
    """
    
    print("\n" + "="*80)
    print("STARTING NETWORK DIFFUSION ANALYSIS")
    print("="*80 + "\n")
    
    # Step 1: Run diffusion
    g_gcc = run_diffusion_all_timepoints(g_gcc)
    
    # Step 2: Compute polarity
    g_gcc, influence_matrices = compute_polarity_all_timepoints(g_gcc)
    
    # Step 3: Compute driver metrics
    g_gcc = compute_driver_metrics_all_timepoints(g_gcc, influence_matrices)
    
    # Step 4: Classify and score
    g_gcc = classify_and_score_drivers(g_gcc)
    
    # Step 5: Build final table
    final_driver_table = build_final_driver_table(g_gcc)
    
    # Step 6: Select top drivers/sinks
    top_early, top_mid, top_late, top_sinks = select_top_drivers_and_sinks(final_driver_table)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80 + "\n")
    
    return {
        "g_gcc": g_gcc,
        "final_driver_table": final_driver_table,
        "top_early": top_early,
        "top_mid": top_mid,
        "top_late": top_late,
        "top_sinks": top_sinks,
        "influence_matrices": influence_matrices
    }


# ============================================================================
# USAGE
# ============================================================================

if __name__ == "__main__":
    """
    To use this script:
    
    1. Load your g_gcc network with node attributes attached
    2. Run the analysis:
    
        results = run_complete_analysis(g_gcc)
        
    3. Export publication data:
    
        from export_publication_data import export_all_publication_data
        
        export_all_publication_data(
            g_gcc=results["g_gcc"],
            final_driver_table=results["final_driver_table"],
            top_early=results["top_early"],
            top_mid=results["top_mid"],
            top_late=results["top_late"],
            top_sinks=results["top_sinks"]
        )
    """
    print(__doc__)
