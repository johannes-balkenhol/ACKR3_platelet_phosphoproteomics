#!/usr/bin/env python3
"""
Publication Export Script: Network Diffusion Analysis
=====================================================

This script exports all tables, figures, and data needed for publication.
It should be run after completing the full network diffusion analysis.

Requirements:
- Completed diffusion analysis with annotated g_gcc graph
- Driver classification tables (top_early, top_mid, top_late, top_sinks)
- Final driver summary table

Outputs:
- CSV tables for supplementary materials
- High-resolution Sankey plots (PNG and HTML)
- Network statistics report
- Annotated network file (GraphML)
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import pickle
from pathlib import Path
from datetime import datetime

# ============================================================================
# CONFIGURATION
# ============================================================================

OUTPUT_DIR = Path("publication_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

# Create subdirectories
(OUTPUT_DIR / "tables").mkdir(exist_ok=True)
(OUTPUT_DIR / "figures").mkdir(exist_ok=True)
(OUTPUT_DIR / "networks").mkdir(exist_ok=True)
(OUTPUT_DIR / "reports").mkdir(exist_ok=True)

# Timepoints to analyze
TIMEPOINTS = ["10", "600", "1800"]

# Manual proteins of interest to include
PROTEINS_OF_INTEREST = ["ARRB1", "FLNA", "TRIO", "PDE5A", "OPHN1"]

# ============================================================================
# FUNCTION: CREATE SANKEY PLOT
# ============================================================================

def create_sankey_plot(g_gcc, tp, top_early, top_mid, top_late, top_sinks, 
                       proteins_of_interest=None):
    """
    Create publication-quality Sankey diagram for network diffusion analysis
    
    Parameters
    ----------
    g_gcc : igraph.Graph
        Giant connected component with attached metrics
    tp : str
        Timepoint ("10", "600", or "1800")
    top_early, top_mid, top_late, top_sinks : pd.DataFrame
        Top ranked driver/sink tables
    proteins_of_interest : list, optional
        Additional protein symbols to include
    
    Returns
    -------
    fig : plotly.graph_objects.Figure
        Sankey diagram
    selected_nodes : list
        List of protein IDs included in the plot
    """
    
    # Settings
    logfc_col = f"mean_logFC_{tp}"
    infl_col  = f"ppr_{tp}"
    score_col = {"10": "early_score", 
                 "600": "mid_score", 
                 "1800": "late_score"}[tp]
    
    # Select nodes
    selected = (
        set(top_early["UniProt"].dropna()) |
        set(top_mid["UniProt"].dropna()) |
        set(top_late["UniProt"].dropna()) |
        set(top_sinks["UniProt"].dropna())
    )
    
    # Add proteins of interest
    if proteins_of_interest:
        for sym in proteins_of_interest:
            hits = [v["name"] for v in g_gcc.vs if v.get("gene_symbol") == sym]
            selected.update(hits)
    
    # Filter to nodes actually in the graph
    selected_nodes = [n for n in selected if n in g_gcc.vs["name"]]
    node_index = {n: i for i, n in enumerate(selected_nodes)}
    
    print(f"  Timepoint {tp}: {len(selected_nodes)} nodes selected")
    
    # Build node attributes
    node_labels = []
    node_colors = []
    node_sizes  = []
    
    for n in selected_nodes:
        v = g_gcc.vs.find(name=n)
        
        # Label
        node_labels.append(v.get("gene_symbol", n))
        
        # Color from signed mean logFC
        fc = v.get(logfc_col)
        if fc is None or np.isnan(fc):
            node_colors.append("rgba(200,200,200,0.4)")
        elif fc > 0:
            alpha = min(0.85, 0.3 + abs(fc)/3)
            node_colors.append(f"rgba(255,0,0,{alpha})")
        else:
            alpha = min(0.85, 0.3 + abs(fc)/3)
            node_colors.append(f"rgba(0,0,255,{alpha})")
        
        # Node size from driver score
        score = v.get(score_col)
        if score is None or np.isnan(score):
            size = 0.1
        else:
            size = abs(score) + 0.2
        node_sizes.append(size)
    
    # Build edge attributes
    sources = []
    targets = []
    values = []
    edge_colors = []
    
    for e in g_gcc.es:
        src = g_gcc.vs[e.source]["name"]
        tgt = g_gcc.vs[e.target]["name"]
        
        if src not in node_index or tgt not in node_index:
            continue
        
        # Edge weight from influence
        infl = g_gcc.vs[e.source].get(infl_col)
        w = 0.05 if infl is None or np.isnan(infl) else max(0.05, abs(infl))
        
        # Color from edge type
        edge_type = e.get("type", "other")
        if edge_type == "activation":
            col = "rgba(255,0,0,0.35)"
        elif edge_type == "inhibition":
            col = "rgba(0,0,255,0.35)"
        else:
            col = "rgba(150,150,150,0.25)"
        
        sources.append(node_index[src])
        targets.append(node_index[tgt])
        values.append(w)
        edge_colors.append(col)
    
    print(f"  Timepoint {tp}: {len(sources)} edges included")
    
    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=20,
            thickness=20,
            label=node_labels,
            color=node_colors,
            customdata=node_sizes,
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=edge_colors,
        )
    )])
    
    fig.update_layout(
        title=f"Network Diffusion Analysis – t={tp} min",
        width=1100,
        height=900,
        font=dict(size=14)
    )
    
    return fig, selected_nodes

# ============================================================================
# FUNCTION: EXPORT TABLES
# ============================================================================

def export_tables(final_driver_table, top_early, top_mid, top_late, top_sinks):
    """Export all driver analysis tables as CSV files"""
    
    print("\n" + "="*80)
    print("EXPORTING TABLES")
    print("="*80)
    
    # Main driver summary table
    output_file = OUTPUT_DIR / "tables" / "supplementary_table_driver_summary.csv"
    final_driver_table.to_csv(output_file, index=False)
    print(f"✓ Saved: {output_file}")
    print(f"  Rows: {len(final_driver_table)}, Columns: {len(final_driver_table.columns)}")
    
    # Top drivers by timepoint
    tables = {
        "top_early_drivers.csv": top_early,
        "top_mid_drivers.csv": top_mid,
        "top_late_drivers.csv": top_late,
        "top_sink_nodes.csv": top_sinks
    }
    
    for filename, df in tables.items():
        output_file = OUTPUT_DIR / "tables" / filename
        df.to_csv(output_file, index=False)
        print(f"✓ Saved: {output_file}")
        print(f"  Rows: {len(df)}")
    
    # Create a summary statistics table
    stats_data = []
    for tp in TIMEPOINTS:
        n_upregulated = len(final_driver_table[
            final_driver_table[f"mean_logFC_{tp}"] > 0
        ])
        n_downregulated = len(final_driver_table[
            final_driver_table[f"mean_logFC_{tp}"] < 0
        ])
        
        stats_data.append({
            "Timepoint_min": tp,
            "Total_Proteins": len(final_driver_table),
            "Upregulated": n_upregulated,
            "Downregulated": n_downregulated,
            "Early_Drivers": len(top_early),
            "Mid_Drivers": len(top_mid),
            "Late_Drivers": len(top_late),
            "Sink_Nodes": len(top_sinks)
        })
    
    stats_df = pd.DataFrame(stats_data)
    output_file = OUTPUT_DIR / "tables" / "summary_statistics.csv"
    stats_df.to_csv(output_file, index=False)
    print(f"✓ Saved: {output_file}")
    
    return stats_df

# ============================================================================
# FUNCTION: EXPORT FIGURES
# ============================================================================

def export_figures(g_gcc, top_early, top_mid, top_late, top_sinks):
    """Generate and export all Sankey plots"""
    
    print("\n" + "="*80)
    print("GENERATING SANKEY PLOTS")
    print("="*80)
    
    all_nodes = {}
    
    for tp in TIMEPOINTS:
        print(f"\nGenerating plot for t={tp} min...")
        
        fig, selected_nodes = create_sankey_plot(
            g_gcc, tp, top_early, top_mid, top_late, top_sinks,
            proteins_of_interest=PROTEINS_OF_INTEREST
        )
        
        all_nodes[tp] = selected_nodes
        
        # Save as interactive HTML
        html_file = OUTPUT_DIR / "figures" / f"sankey_t{tp}_interactive.html"
        fig.write_html(str(html_file))
        print(f"✓ Saved HTML: {html_file}")
        
        # Save as high-resolution PNG
        try:
            png_file = OUTPUT_DIR / "figures" / f"sankey_t{tp}_publication.png"
            fig.write_image(str(png_file), width=1100, height=900, scale=3)
            print(f"✓ Saved PNG: {png_file}")
        except Exception as e:
            print(f"✗ Could not save PNG (install kaleido): {e}")
        
        # Save as SVG for vector graphics
        try:
            svg_file = OUTPUT_DIR / "figures" / f"sankey_t{tp}_publication.svg"
            fig.write_image(str(svg_file), width=1100, height=900)
            print(f"✓ Saved SVG: {svg_file}")
        except Exception as e:
            print(f"✗ Could not save SVG: {e}")
    
    # Save list of nodes included in each plot
    nodes_df = pd.DataFrame({
        f"t{tp}_nodes": pd.Series(all_nodes[tp]) for tp in TIMEPOINTS
    })
    nodes_file = OUTPUT_DIR / "tables" / "sankey_nodes_included.csv"
    nodes_df.to_csv(nodes_file, index=False)
    print(f"\n✓ Saved node list: {nodes_file}")
    
    return all_nodes

# ============================================================================
# FUNCTION: EXPORT NETWORK
# ============================================================================

def export_network(g_gcc):
    """Export annotated network in multiple formats"""
    
    print("\n" + "="*80)
    print("EXPORTING NETWORK FILES")
    print("="*80)
    
    # Save as pickle (Python)
    pickle_file = OUTPUT_DIR / "networks" / "g_gcc_annotated.pkl"
    with open(pickle_file, "wb") as f:
        pickle.dump(g_gcc, f)
    print(f"✓ Saved pickle: {pickle_file}")
    
    # Save as GraphML (Cytoscape compatible)
    graphml_file = OUTPUT_DIR / "networks" / "g_gcc_annotated.graphml"
    g_gcc.write_graphml(str(graphml_file))
    print(f"✓ Saved GraphML: {graphml_file}")
    
    # Export edge list
    edge_list = []
    for e in g_gcc.es:
        edge_list.append({
            "source": g_gcc.vs[e.source]["name"],
            "target": g_gcc.vs[e.target]["name"],
            "source_symbol": g_gcc.vs[e.source].get("gene_symbol", ""),
            "target_symbol": g_gcc.vs[e.target].get("gene_symbol", ""),
            "type": e.get("type", "")
        })
    
    edge_df = pd.DataFrame(edge_list)
    edge_file = OUTPUT_DIR / "networks" / "edge_list.csv"
    edge_df.to_csv(edge_file, index=False)
    print(f"✓ Saved edge list: {edge_file}")
    print(f"  Total edges: {len(edge_df)}")
    
    return edge_df

# ============================================================================
# FUNCTION: GENERATE REPORT
# ============================================================================

def generate_report(g_gcc, final_driver_table, stats_df, edge_df):
    """Generate a comprehensive analysis report"""
    
    print("\n" + "="*80)
    print("GENERATING ANALYSIS REPORT")
    print("="*80)
    
    report = []
    report.append("="*80)
    report.append("NETWORK DIFFUSION ANALYSIS - PUBLICATION REPORT")
    report.append("="*80)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")
    
    # Network statistics
    report.append("NETWORK STATISTICS")
    report.append("-"*80)
    report.append(f"Total nodes in GCC: {g_gcc.vcount()}")
    report.append(f"Total edges in GCC: {g_gcc.ecount()}")
    report.append(f"Network density: {g_gcc.density():.4f}")
    report.append("")
    
    # Edge type distribution
    edge_types = edge_df["type"].value_counts()
    report.append("Edge Type Distribution:")
    for edge_type, count in edge_types.items():
        report.append(f"  {edge_type}: {count} ({count/len(edge_df)*100:.1f}%)")
    report.append("")
    
    # Driver statistics
    report.append("DRIVER CLASSIFICATION SUMMARY")
    report.append("-"*80)
    for _, row in stats_df.iterrows():
        report.append(f"Timepoint: {row['Timepoint_min']} min")
        report.append(f"  Total proteins analyzed: {row['Total_Proteins']}")
        report.append(f"  Upregulated: {row['Upregulated']}")
        report.append(f"  Downregulated: {row['Downregulated']}")
        report.append(f"  Early drivers: {row['Early_Drivers']}")
        report.append(f"  Mid drivers: {row['Mid_Drivers']}")
        report.append(f"  Late drivers: {row['Late_Drivers']}")
        report.append(f"  Sink nodes: {row['Sink_Nodes']}")
        report.append("")
    
    # Top driver genes
    report.append("TOP DRIVER GENES BY CATEGORY")
    report.append("-"*80)
    
    for category, score_col in [
        ("Early Drivers", "early_score"),
        ("Mid Drivers", "mid_score"),
        ("Late Drivers", "late_score")
    ]:
        top_genes = final_driver_table.nlargest(10, score_col)
        report.append(f"\n{category} (Top 10 by {score_col}):")
        for i, (_, row) in enumerate(top_genes.iterrows(), 1):
            symbol = row.get("gene_symbol", row["UniProt"])
            score = row[score_col]
            report.append(f"  {i}. {symbol} (score: {score:.3f})")
    
    # Save report
    report_text = "\n".join(report)
    report_file = OUTPUT_DIR / "reports" / "analysis_report.txt"
    with open(report_file, "w") as f:
        f.write(report_text)
    
    print(f"✓ Saved report: {report_file}")
    print("\n" + report_text)
    
    return report_text

# ============================================================================
# MAIN EXPORT FUNCTION
# ============================================================================

def export_all_publication_data(g_gcc, final_driver_table, top_early, top_mid, 
                                top_late, top_sinks):
    """
    Main function to export all publication data
    
    This is the function you call after completing your diffusion analysis
    """
    
    print("\n" + "="*80)
    print("STARTING PUBLICATION DATA EXPORT")
    print("="*80)
    print(f"Output directory: {OUTPUT_DIR.absolute()}")
    
    # Export tables
    stats_df = export_tables(final_driver_table, top_early, top_mid, 
                            top_late, top_sinks)
    
    # Generate and export figures
    all_nodes = export_figures(g_gcc, top_early, top_mid, top_late, top_sinks)
    
    # Export network
    edge_df = export_network(g_gcc)
    
    # Generate report
    report = generate_report(g_gcc, final_driver_table, stats_df, edge_df)
    
    print("\n" + "="*80)
    print("EXPORT COMPLETE!")
    print("="*80)
    print(f"\nAll files saved to: {OUTPUT_DIR.absolute()}")
    print("\nDirectory structure:")
    print(f"  {OUTPUT_DIR}/")
    print(f"    ├── tables/")
    print(f"    │   ├── supplementary_table_driver_summary.csv")
    print(f"    │   ├── top_early_drivers.csv")
    print(f"    │   ├── top_mid_drivers.csv")
    print(f"    │   ├── top_late_drivers.csv")
    print(f"    │   ├── top_sink_nodes.csv")
    print(f"    │   ├── summary_statistics.csv")
    print(f"    │   └── sankey_nodes_included.csv")
    print(f"    ├── figures/")
    print(f"    │   ├── sankey_t10_interactive.html")
    print(f"    │   ├── sankey_t10_publication.png")
    print(f"    │   ├── sankey_t600_interactive.html")
    print(f"    │   ├── sankey_t600_publication.png")
    print(f"    │   ├── sankey_t1800_interactive.html")
    print(f"    │   └── sankey_t1800_publication.png")
    print(f"    ├── networks/")
    print(f"    │   ├── g_gcc_annotated.pkl")
    print(f"    │   ├── g_gcc_annotated.graphml")
    print(f"    │   └── edge_list.csv")
    print(f"    └── reports/")
    print(f"        └── analysis_report.txt")
    
    return OUTPUT_DIR

# ============================================================================
# USAGE EXAMPLE
# ============================================================================

if __name__ == "__main__":
    """
    To use this script, ensure you have:
    1. Loaded your g_gcc graph with all attributes attached
    2. Created final_driver_table
    3. Created top_early, top_mid, top_late, top_sinks DataFrames
    
    Then simply run:
    
    export_all_publication_data(
        g_gcc=g_gcc,
        final_driver_table=final_driver_table,
        top_early=top_early,
        top_mid=top_mid,
        top_late=top_late,
        top_sinks=top_sinks
    )
    """
    
    print(__doc__)
    print("\nThis script is ready to use. Import it in your Jupyter notebook:")
    print("\n  from export_publication_data import export_all_publication_data")
    print("\n  export_all_publication_data(")
    print("      g_gcc=g_gcc,")
    print("      final_driver_table=final_driver_table,")
    print("      top_early=top_early,")
    print("      top_mid=top_mid,")
    print("      top_late=top_late,")
    print("      top_sinks=top_sinks")
    print("  )")
