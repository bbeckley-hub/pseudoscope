#!/usr/bin/env python3
"""
PseudoScope Comprehensive Visualizer
================================================================================
Generates a rich, vertically-scrolling HTML dashboard with:
- ST distribution (bar chart)
- O‑type distribution (bar chart)
- Top 20-50 AMR genes (stacked by database) – can change to top 100 by editing variable
- Top 20-50 virulence genes (stacked by database)
- Co‑occurrence network (top 100 edges)
- FASTA QC box plots (GC%, AT%, N50, Total Length) with points & outliers (using TSV)
- PCA plot of gene presence/absence (colored by ST) – gracefully handles single sample
- Co‑occurrence heatmap (top 30 AMR genes)

All plots are interactive, include scientific explanations, and use varied colors.

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.0.0
Date: 2026-04-22
================================================================================
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Helper functions
# =============================================================================

def parse_qc_from_tsv(tsv_path):
    """Extract QC metrics from TSV file."""
    df = pd.read_csv(tsv_path, sep='\t')
    # Rename columns to simpler names
    df.rename(columns={'Filename': 'Sample'}, inplace=True)
    df['Sample'] = df['Sample'].apply(lambda x: Path(x).stem)
    # Keep only needed metrics
    keep = ['Sample', 'GC Content (%)', 'AT Content (%)', 'N50', 'Total Length']
    existing = [c for c in keep if c in df.columns]
    df = df[existing]
    for col in ['GC Content (%)', 'AT Content (%)', 'N50', 'Total Length']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def load_sample_overview(csv_path):
    """Load sample overview (ST, O‑type)."""
    df = pd.read_csv(csv_path)
    if 'Sample' not in df.columns:
        df.rename(columns={df.columns[0]: 'Sample'}, inplace=True)
    return df

def plot_st_distribution(df):
    """Bar chart of ST counts."""
    st_counts = df['ST'].value_counts().reset_index()
    st_counts.columns = ['ST', 'Count']
    fig = px.bar(st_counts, x='ST', y='Count', title='Sequence Type (ST) Distribution',
                 labels={'Count': 'Number of genomes', 'ST': 'ST'},
                 color='Count', color_continuous_scale='Blues')
    fig.update_layout(xaxis_tickangle=-45, height=500)
    return fig

def plot_otype_distribution(df):
    """Bar chart of O‑type counts."""
    o_counts = df['O_Type'].value_counts().reset_index()
    o_counts.columns = ['O_Type', 'Count']
    fig = px.bar(o_counts, x='O_Type', y='Count', title='O‑Serotype Distribution',
                 labels={'Count': 'Number of genomes', 'O_Type': 'O‑type'},
                 color='Count', color_continuous_scale='Greens')
    fig.update_layout(xaxis_tickangle=-45, height=500)
    return fig

def plot_top_genes_stacked(csv_path, title, top_n=50, color_seq=None):
    """Stacked bar chart of top genes by database."""
    if not csv_path.exists():
        return None
    df = pd.read_csv(csv_path)
    # Sum counts per gene to get top N
    gene_total = df.groupby('Gene')['Count'].sum().nlargest(top_n).index
    df_top = df[df['Gene'].isin(gene_total)]
    # Use a visible qualitative palette
    if color_seq is None:
        color_seq = px.colors.qualitative.Set1 + px.colors.qualitative.Dark2
    fig = px.bar(df_top, x='Gene', y='Count', color='Database', title=title,
                 labels={'Count': 'Number of genomes', 'Gene': 'Gene'},
                 barmode='stack', color_discrete_sequence=color_seq)
    fig.update_layout(xaxis_tickangle=-45, height=550, legend_title='Database')
    return fig

def build_network_from_cooccurrence(csv_path, top_edges=100):
    """Build network graph from gene_cooccurrence.csv."""
    df = pd.read_csv(csv_path)
    if 'Gene1' not in df.columns or 'Gene2' not in df.columns:
        if 'Gene 1' in df.columns:
            df.rename(columns={'Gene 1': 'Gene1', 'Gene 2': 'Gene2'}, inplace=True)
    df_top = df.nlargest(top_edges, 'Count')
    G = nx.Graph()
    for _, row in df_top.iterrows():
        G.add_edge(row['Gene1'], row['Gene2'], weight=row['Count'])
    return G

def plot_network(G):
    """Interactive network plot."""
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    node_x = [pos[n][0] for n in G.nodes()]
    node_y = [pos[n][1] for n in G.nodes()]
    node_degree = [G.degree(n) for n in G.nodes()]
    hover_text = [f"Gene: {n}<br>Degree: {d}" for n, d in zip(G.nodes(), node_degree)]
    edge_trace = go.Scatter(x=edge_x, y=edge_y, mode='lines',
                            line=dict(width=0.8, color='lightgray'), hoverinfo='none')
    node_trace = go.Scatter(x=node_x, y=node_y, mode='markers+text',
                            marker=dict(size=[8 + d*1.5 for d in node_degree],
                                        color='#8B4513', line=dict(width=1, color='white')),
                            text=[n[:12] for n in G.nodes()], textposition='top center',
                            hovertext=hover_text, hoverinfo='text')
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(title=f'Gene Co‑occurrence Network (top {len(G.edges())} edges)',
                      showlegend=False, xaxis_visible=False, yaxis_visible=False,
                      plot_bgcolor='white', height=700)
    return fig

def plot_qc_boxplots_with_points(qc_df):
    """Four separate box plots with all points (jitter) and outliers."""
    if qc_df.empty:
        return None
    metrics = ['GC Content (%)', 'AT Content (%)', 'N50', 'Total Length']
    available = [m for m in metrics if m in qc_df.columns]
    if not available:
        return None
    fig = make_subplots(rows=2, cols=2, subplot_titles=available,
                        shared_yaxes=False)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    for i, metric in enumerate(available):
        row = i // 2 + 1
        col = i % 2 + 1
        # Remove any NaN values
        values = qc_df[metric].dropna()
        if len(values) == 0:
            continue
        box = go.Box(y=values, name=metric, boxpoints='all', jitter=0.3,
                     pointpos=-1.8, marker=dict(color=colors[i], size=4),
                     line=dict(color=colors[i]), fillcolor='rgba(0,0,0,0)')
        fig.add_trace(box, row=row, col=col)
        if metric in ['N50', 'Total Length']:
            fig.update_yaxes(type='log', row=row, col=col)
    fig.update_layout(height=800, title_text="FASTA QC Metrics Distribution (with all points)",
                      showlegend=False)
    return fig

def plot_pca(amr_csv, vir_csv, sample_df):
    """PCA of gene presence/absence (AMR + virulence) colored by ST."""
    # Check for at least 2 samples
    if len(sample_df) < 2:
        print("Skipping PCA: need at least 2 samples.")
        return None
    if not amr_csv.exists() or not vir_csv.exists():
        return None
    amr = pd.read_csv(amr_csv)
    vir = pd.read_csv(vir_csv)
    all_genes = set(amr['Gene'].unique()) | set(vir['Gene'].unique())
    samples = sample_df['Sample'].tolist()
    sample_genes = {s: set() for s in samples}
    for _, row in amr.iterrows():
        for s in row['Genomes'].split(';'):
            if s in sample_genes:
                sample_genes[s].add(row['Gene'])
    for _, row in vir.iterrows():
        for s in row['Genomes'].split(';'):
            if s in sample_genes:
                sample_genes[s].add(row['Gene'])
    gene_list = list(all_genes)
    matrix = np.zeros((len(samples), len(gene_list)), dtype=int)
    for i, s in enumerate(samples):
        for j, g in enumerate(gene_list):
            if g in sample_genes[s]:
                matrix[i, j] = 1
    if matrix.shape[1] < 2:
        return None
    pca = PCA(n_components=2, random_state=42)
    scaled = StandardScaler(with_mean=False).fit_transform(matrix)
    pca_result = pca.fit_transform(scaled)
    pca_df = pd.DataFrame({'PC1': pca_result[:,0], 'PC2': pca_result[:,1],
                           'Sample': samples})
    pca_df = pca_df.merge(sample_df[['Sample', 'ST']], on='Sample', how='left')
    fig = px.scatter(pca_df, x='PC1', y='PC2', color='ST',
                     title=f'PCA of Gene Presence/Absence (AMR+Virulence)',
                     labels={'PC1': f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)',
                             'PC2': f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)'},
                     hover_data=['Sample'],
                     color_discrete_sequence=px.colors.qualitative.Set1)
    fig.update_traces(marker=dict(size=10, line=dict(width=1, color='DarkSlateGrey')))
    fig.update_layout(height=600)
    return fig

def plot_cooccurrence_heatmap(amr_csv, top_n=30):
    """Co‑occurrence heatmap of top AMR genes."""
    if not amr_csv.exists():
        return None
    amr = pd.read_csv(amr_csv)
    top_genes = amr.groupby('Gene')['Count'].sum().nlargest(top_n).index.tolist()
    if len(top_genes) < 2:
        return None
    n = len(top_genes)
    cooc = np.zeros((n, n), dtype=int)
    gene_to_idx = {g:i for i,g in enumerate(top_genes)}
    genome_genes = {}
    for _, row in amr.iterrows():
        if row['Gene'] not in top_genes:
            continue
        for genome in row['Genomes'].split(';'):
            if genome not in genome_genes:
                genome_genes[genome] = set()
            genome_genes[genome].add(row['Gene'])
    for genes in genome_genes.values():
        gene_list = [g for g in genes if g in top_genes]
        for i, g1 in enumerate(gene_list):
            for g2 in gene_list[i+1:]:
                cooc[gene_to_idx[g1], gene_to_idx[g2]] += 1
                cooc[gene_to_idx[g2], gene_to_idx[g1]] += 1
    fig = px.imshow(cooc, x=top_genes, y=top_genes, text_auto=True,
                    title=f'AMR Gene Co‑occurrence Heatmap (Top {top_n})',
                    color_continuous_scale='RdBu', aspect='auto')
    fig.update_layout(height=700, xaxis_tickangle=-45)
    return fig

# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Generate PseudoScope visual dashboard')
    parser.add_argument('-i', '--input-dir', default='.',
                        help='Directory containing input files (CSVs and QC TSV/HTML)')
    args = parser.parse_args()

    current_dir = Path(args.input_dir)
    if not current_dir.exists():
        print(f"ERROR: Directory {current_dir} not found.")
        sys.exit(1)

    # Input files
    sample_csv = current_dir / "sample_overview.csv"
    amr_csv = current_dir / "amr_genes.csv"
    vir_csv = current_dir / "virulence_genes.csv"
    cooc_csv = current_dir / "gene_cooccurrence.csv"
    qc_tsv = current_dir / "PseudoScope_FASTA_QC_summary.tsv"
    qc_html = current_dir / "PseudoScope_FASTA_QC_summary.html"

    if not sample_csv.exists():
        print("ERROR: sample_overview.csv not found.")
        sys.exit(1)

    print("Loading sample data...")
    sample_df = load_sample_overview(sample_csv)

    print("Creating ST distribution plot...")
    fig_st = plot_st_distribution(sample_df)

    print("Creating O‑type distribution plot...")
    fig_otype = plot_otype_distribution(sample_df)

    print("Creating AMR stacked bar chart (top 50)...")
    fig_amr_stack = plot_top_genes_stacked(amr_csv, "Top 50 AMR Genes (by database)",
                                           color_seq=px.colors.qualitative.Set1 + px.colors.qualitative.Dark2)

    print("Creating virulence stacked bar chart (top 50)...")
    fig_vir_stack = plot_top_genes_stacked(vir_csv, "Top 50 Virulence Genes (by database)",
                                           color_seq=px.colors.qualitative.Set2 + px.colors.qualitative.Pastel1)

    print("Building co‑occurrence network (top 100 edges)...")
    if cooc_csv.exists():
        G = build_network_from_cooccurrence(cooc_csv, top_edges=100)
        fig_network = plot_network(G)
    else:
        fig_network = None

    print("Parsing QC data...")
    if qc_tsv.exists():
        qc_df = parse_qc_from_tsv(qc_tsv)
    elif qc_html.exists():
        # Fallback to HTML parser (simpler to use TSV, but keeping for compatibility)
        from bs4 import BeautifulSoup
        def parse_qc_from_html_fallback(html_path):
            with open(html_path, 'r', encoding='utf-8') as f:
                soup = BeautifulSoup(f.read(), 'html.parser')
            tables = soup.find_all('table')
            table = tables[1] if len(tables) >= 2 else tables[0]
            if not table:
                return pd.DataFrame()
            rows = table.find_all('tr')
            if len(rows) < 2:
                return pd.DataFrame()
            headers = [th.get_text().strip() for th in rows[0].find_all('th')]
            data = []
            for row in rows[1:]:
                cols = row.find_all('td')
                if cols:
                    data.append([col.get_text().strip() for col in cols])
            if not data:
                return pd.DataFrame()
            df = pd.DataFrame(data, columns=headers)
            if 'Filename' in df.columns:
                df['Sample'] = df['Filename'].apply(lambda x: Path(x).stem)
            else:
                df['Sample'] = df.iloc[:,0].apply(lambda x: Path(x).stem)
            keep = ['Sample', 'GC Content (%)', 'AT Content (%)', 'N50', 'Total Length']
            existing = [c for c in keep if c in df.columns]
            df = df[existing]
            for col in ['GC Content (%)', 'AT Content (%)', 'N50', 'Total Length']:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            return df
        qc_df = parse_qc_from_html_fallback(qc_html)
    else:
        qc_df = pd.DataFrame()

    if not qc_df.empty:
        fig_qc = plot_qc_boxplots_with_points(qc_df)
    else:
        fig_qc = None

    print("Computing PCA...")
    fig_pca = plot_pca(amr_csv, vir_csv, sample_df)

    print("Creating co‑occurrence heatmap (top 30 AMR genes)...")
    fig_heatmap = plot_cooccurrence_heatmap(amr_csv, top_n=30)

    # Count total samples
    total_samples = len(sample_df)
    st_count = sample_df['ST'].nunique()
    o_count = sample_df['O_Type'].nunique()
    amr_total = len(pd.read_csv(amr_csv)['Gene'].unique()) if amr_csv.exists() else 0
    vir_total = len(pd.read_csv(vir_csv)['Gene'].unique()) if vir_csv.exists() else 0

    # Assemble HTML (vertical scrolling, clean header/footer)
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    html_parts = []
    html_parts.append(f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <title>PseudoScope Comprehensive Dashboard</title>
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; background: #f0f2f5; }}
            .container {{ max-width: 1400px; margin: 20px auto; background: white; padding: 30px; border-radius: 20px; box-shadow: 0 0 30px rgba(0,0,0,0.1); }}
            .main-header {{ background: linear-gradient(135deg, #5C2E0B 0%, #8B4513 100%); color: white; padding: 20px 30px; border-radius: 15px; margin-bottom: 30px; }}
            .main-header h1 {{ margin: 0 0 10px 0; color: white; }}
            .metadata-bar {{ background: rgba(255,255,255,0.1); padding: 12px; border-radius: 10px; margin-top: 15px; display: flex; justify-content: space-around; flex-wrap: wrap; gap: 15px; font-size: 0.9em; }}
            h2 {{ color: #5C2E0B; margin-top: 40px; border-left: 5px solid #8B4513; padding-left: 15px; }}
            .explanation {{ background: #f9f9e0; border-left: 5px solid #8B4513; padding: 15px; margin: 20px 0; border-radius: 8px; font-size: 0.95em; }}
            .footer {{ text-align: center; padding: 15px; margin-top: 50px; background: linear-gradient(135deg, #5C2E0B, #8B4513); color: white; border-radius: 12px; font-size: 0.85em; }}
            .footer a {{ color: #ffd700; text-decoration: none; }}
            .footer a:hover {{ text-decoration: underline; }}
        </style>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>
    <body>
    <div class="container">
        <div class="main-header">
            <h1>🧬 PseudoScope Comprehensive Dashboard</h1>
            <p>Interactive visualisation of <em>Pseudomonas aeruginosa</em> genomic data – ST and O‑type distributions, resistance/virulence gene prevalence, co‑occurrence networks, assembly quality, and multivariate analysis.</p>
            <div class="metadata-bar">
                <div><strong>📅 Analysis Date:</strong> {now}</div>
                <div><strong>🧬 Samples:</strong> {total_samples}</div>
                <div><strong>🔬 STs:</strong> {st_count}</div>
                <div><strong>🏷️ O‑types:</strong> {o_count}</div>
                <div><strong>💊 AMR Genes:</strong> {amr_total}</div>
                <div><strong>🦠 Virulence Genes:</strong> {vir_total}</div>
            </div>
        </div>
    """)

    # Section: Typing
    html_parts.append("<h2>📊 MLST & Serotype Distribution</h2>")
    html_parts.append('<div class="explanation"><strong>ST and O‑type prevalence</strong> – Bars show the number of genomes per Sequence Type (Pasteur scheme) and O‑serotype. Hover for exact counts.</div>')
    html_parts.append(fig_st.to_html(full_html=False, include_plotlyjs='cdn'))
    html_parts.append(fig_otype.to_html(full_html=False, include_plotlyjs='cdn'))

    # Section: Top Genes
    html_parts.append("<h2>💊 Top AMR & Virulence Genes</h2>")
    html_parts.append('<div class="explanation"><strong>Stacked bar charts</strong> – Each bar represents a gene; segments show contributions from different databases (AMRfinder, NCBI, CARD, ResFinder, etc.). Only the top 50 most prevalent genes are shown (adjustable in script).</div>')
    if fig_amr_stack:
        html_parts.append(fig_amr_stack.to_html(full_html=False, include_plotlyjs='cdn'))
    if fig_vir_stack:
        html_parts.append(fig_vir_stack.to_html(full_html=False, include_plotlyjs='cdn'))

    # Section: Co‑occurrence Network
    html_parts.append("<h2>🔗 Gene Co‑occurrence Network</h2>")
    html_parts.append('<div class="explanation"><strong>Network of gene co‑occurrence</strong> – Nodes are genes; edges connect genes that appear together in at least one genome. Edge count is based on the top 100 co‑occurrence pairs. Node size reflects degree (number of connections). Hover for details.</div>')
    if fig_network:
        html_parts.append(fig_network.to_html(full_html=False, include_plotlyjs='cdn'))
    else:
        html_parts.append("<p>Co‑occurrence network data not available.</p>")

    # Section: QC Box Plots
    html_parts.append("<h2>📈 FASTA Assembly Quality</h2>")
    html_parts.append('<div class="explanation"><strong>QC metrics distribution</strong> – Each box plot shows the spread of GC%, AT%, N50 (contiguity), and total assembly length. All individual data points are shown as dots (jittered). N50 and Total Length are on a log scale.</div>')
    if fig_qc:
        html_parts.append(fig_qc.to_html(full_html=False, include_plotlyjs='cdn'))
    else:
        html_parts.append("<p>QC data not found. Please ensure PseudoScope_FASTA_QC_summary.tsv is present.</p>")

    # Section: PCA
    html_parts.append("<h2>📐 Principal Component Analysis (PCA)</h2>")
    html_parts.append('<div class="explanation"><strong>Gene presence/absence PCA</strong> – Each point is a genome, colored by ST. The plot reduces the high‑dimensional gene profile to two principal components. Clusters indicate similar resistance/virulence gene complements.</div>')
    if fig_pca:
        html_parts.append(fig_pca.to_html(full_html=False, include_plotlyjs='cdn'))
    else:
        html_parts.append("<p>⚠️ PCA plot not available (need at least 2 samples or insufficient gene data).</p>")

    # Section: Co‑occurrence Heatmap
    html_parts.append("<h2>🔥 AMR Gene Co‑occurrence Heatmap</h2>")
    html_parts.append('<div class="explanation"><strong>Pairwise co‑occurrence counts</strong> – Darker red indicates that two AMR genes appear together more often. Helps identify linked resistance determinants. Shown for top 30 AMR genes.</div>')
    if fig_heatmap:
        html_parts.append(fig_heatmap.to_html(full_html=False, include_plotlyjs='cdn'))
    else:
        html_parts.append("<p>Heatmap not available (insufficient AMR genes).</p>")

    # Footer (brown gradient)
    html_parts.append(f"""
        <div class="footer">
            <strong>PseudoScope Comprehensive Dashboard v1.0.0</strong><br>
            University of Ghana Medical School | Department of Medical Biochemistry<br>
            Author: Brown Beckley &lt;brownbeckley94@gmail.com&gt; | GitHub: <a href="https://github.com/bbeckley-hub/pseudoscope" target="_blank">bbeckley-hub/pseudoscope</a><br>
            Generated on {now}
        </div>
    </div>
    </body>
    </html>
    """)

    output_path = current_dir / "genius_pseudomonas_visual_dashboard.html"
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("\n".join(html_parts))

    print(f"✅ Dashboard saved: {output_path}")
    print("Open the file in your browser to explore all interactive plots.")

if __name__ == "__main__":
    main()