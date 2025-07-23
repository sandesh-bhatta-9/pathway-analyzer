import streamlit as st
import networkx as nx
import pandas as pd
from collections import Counter
from bioservices import KEGG
import plotly.graph_objects as go
import base64
from pathlib import Path

# --- Page Configuration ---
st.set_page_config(layout="wide")
st.title('Interactive Pathway Overlap Analyzer')
st.write("This tool analyzes pathways from KEGG, provides detailed, linked gene tables, and generates a dynamic network graph.")

# --- Constants ---
ICON_NAMES = ["Gene Cards", "NCBI", "ENSEMBL", "KEGG", "GEO"]
ICONS_DIR = Path(__file__).parent / "Icons"  # folder alongside this script

# --- Cached: load all icons as base64 strings once ---
@st.cache_data
def load_icons_b64():
    icon_map = {}
    for name in ICON_NAMES:
        fp = ICONS_DIR / f"{name}.png"
        if fp.is_file():
            with open(fp, "rb") as f:
                b64 = base64.b64encode(f.read()).decode("utf-8")
            icon_map[name] = f"data:image/png;base64,{b64}"
        else:
            icon_map[name] = None  # fallback to text if missing
    return icon_map

ICON_B64 = load_icons_b64()

# --- KEGG helpers (cached) ---
@st.cache_data
def get_all_kegg_pathways():
    k = KEGG()
    pathways_raw = k.list("pathway/hsa")
    pathway_dict = {}
    for line in pathways_raw.strip().split('\n'):
        path_id, name_desc = line.split('\t')
        name = name_desc.split(' - ')[0]
        pathway_dict[name] = path_id.replace('path:', '')
    return pathway_dict

@st.cache_data
def get_genes(pathway_id: str):
    k = KEGG()
    k.organism = "hsa"
    raw_data = k.get(pathway_id)
    if not raw_data:
        return set()
    genes = set()
    record = False
    for line in raw_data.strip().split('\n'):
        if line.startswith('GENE'):
            record = True
            line_content = line[12:].strip()
        elif record and line.startswith(' '):
            line_content = line.strip()
        else:
            record = False
            continue
        if line_content:
            parts = line_content.split()
            if len(parts) > 1 and parts[0].isdigit():
                genes.add(parts[1].replace(';', ''))
    return genes

# --- Build icon HTML ---
def generate_icon_links(gene_name: str) -> dict:
    urls = {
        "Gene Cards": f"https://www.genecards.org/Search/Keyword?queryString={gene_name}",
        "NCBI":      f"https://www.ncbi.nlm.nih.gov/gene/?term={gene_name}",
        "ENSEMBL":   f"https://useast.ensembl.org/Search/Results?q={gene_name}",
        "KEGG":      f"https://www.kegg.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&dbkey=genes&keywords={gene_name}",
        "GEO":       f"https://www.ncbi.nlm.nih.gov/gds/?term={gene_name}",
    }

    html = {}
    for db, url in urls.items():
        b64src = ICON_B64.get(db)
        if b64src:
            html[db] = (
                f'<a href="{url}" target="_blank">'
                f'<img src="{b64src}" width="24" height="24" '
                f'style="margin:2px; border-radius:4px; border:1px solid #ccc;">'
                f'</a>'
            )
        else:
            # Fallback text link if icon file missing
            html[db] = f'<a href="{url}" target="_blank">{db}</a>'
    return html

# --- Sidebar pathway selection ---
all_pathways = get_all_kegg_pathways()
st.sidebar.header('Select Pathways for Analysis')

show_all = st.sidebar.checkbox("Show all KEGG pathways (~350)")
if show_all:
    pathway_options = all_pathways
    default_selection = list(pathway_options.keys())[:5]
else:
    keywords = ['cancer', 'ascorbate', 'vitamin d', 'diabetes', 'drug metabolism', 'xenobiotics', 'fluoropyrimidine']
    pathway_options = {
        name: pid for name, pid in all_pathways.items()
        if any(keyword in name.lower() for keyword in keywords)
    }
    default_selection = [
        "Pathways in cancer", "Colorectal cancer", "Pancreatic cancer",
        "Type II diabetes mellitus", "Drug metabolism - cytochrome P450"
    ]
    default_selection = [p for p in default_selection if p in pathway_options]

selected_pathways = st.sidebar.multiselect(
    'Choose pathways:',
    options=list(pathway_options.keys()),
    default=default_selection
)

# --- Main content ---
if selected_pathways:
    with st.spinner('Fetching gene data...'):
        pathway_genes_map = {
            name: get_genes(pathway_options[name])
            for name in selected_pathways
        }

    st.header('Shared Gene Analysis')

    # Flatten + count
    all_genes_flat = [gene for genes in pathway_genes_map.values() for gene in genes]
    gene_counts = Counter(all_genes_flat)

    # Rows: only genes in >=2 pathways
    rows = []
    for gene, freq in gene_counts.items():
        if freq >= 2:
            member_of_pathways = [name for name, genes in pathway_genes_map.items() if gene in genes]
            links = generate_icon_links(gene)
            row = {
                "Gene": gene,
                "Frequency": freq,
                "Pathways": ", ".join(member_of_pathways),
                "Gene Cards": links["Gene Cards"],
                "NCBI": links["NCBI"],
                "ENSEMBL": links["ENSEMBL"],
                "KEGG": links["KEGG"],
                "GEO": links["GEO"],
            }
            rows.append(row)

    if rows:
        df_shared = pd.DataFrame(rows)
        # sort: frequency desc then gene
        df_shared = df_shared.sort_values(by=["Frequency", "Gene"], ascending=[False, True])

        st.subheader("🔗 Shared Genes Table with Links")
        st.write(
            df_shared.to_html(escape=False, index=False),
            unsafe_allow_html=True
        )
    else:
        st.write("No shared genes found for the selected pathways.")

    # --- Network plot ---
    st.header('Interactive Network Visualization')
    st.info('You can zoom, pan, and hover over nodes to see their names.')

    with st.spinner('Generating interactive graph...'):
        G = nx.Graph()
        # Pathway nodes
        for name in pathway_genes_map.keys():
            G.add_node(name, type='pathway')
        # Gene nodes
        for gene, count in gene_counts.items():
            G.add_node(gene, type='gene', count=count)
        # Edges
        for name, genes in pathway_genes_map.items():
            for gene in genes:
                G.add_edge(name, gene)

        pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)

        edge_x, edge_y = [], []
        for n1, n2 in G.edges():
            x0, y0 = pos[n1]
            x1, y1 = pos[n2]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

        edge_trace = go.Scatter(
            x=edge_x,
            y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines'
        )

        node_x, node_y, node_text, node_color = [], [], [], []
        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            data = G.nodes[node]
            if data['type'] == 'pathway':
                node_text.append(f"Pathway: {node}")
                node_color.append('skyblue')
            else:
                count = data['count']
                hover_info = f"Gene: {node}<br>Shared in {count} pathways" if count > 1 else f"Gene: {node}"
                node_text.append(hover_info)
                node_color.append('red' if count > 1 else 'lightgray')

        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            mode='markers',
            hoverinfo='text',
            text=node_text,
            marker=dict(
                showscale=False,
                color=node_color,
                size=10,
                line_width=2
            )
        )

        fig = go.Figure(
            data=[edge_trace, node_trace],
            layout=go.Layout(
                title=dict(
                    text='<br>Interactive Network of Pathways and Genes',
                    font=dict(size=16)
                ),
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20, l=5, r=5, t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
            )
        )
        st.plotly_chart(fig, use_container_width=True)

else:
    st.info("Please select pathways from the sidebar to begin the analysis.")
