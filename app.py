import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from collections import Counter
from bioservices import KEGG
import base64
from pathlib import Path

# --- Page Configuration ---
st.set_page_config(page_title="Pathway Overlap Analyzer", layout="wide")
st.title("🧬 Pathway Overlap Analyzer (Genes & Metabolites)")

# --- Icon and KEGG Utilities ---
ICON_NAMES = ["Gene Cards", "NCBI", "ENSEMBL", "KEGG", "GEO"]
ICONS_DIR = Path(__file__).parent / "Icons"

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
            icon_map[name] = None
    return icon_map

ICON_B64 = load_icons_b64()

@st.cache_data
def get_all_kegg_pathways():
    k = KEGG(); pathways_raw = k.list("pathway/hsa")
    d = {};
    for line in pathways_raw.strip().split('\n'):
        pid, name_desc = line.split('\t'); name = name_desc.split(' - ')[0]
        d[name] = pid.replace('path:', '')
    return d

@st.cache_data
def get_pathway_data(pathway_id: str):
    k = KEGG(); k.organism = "hsa"; data = k.get(pathway_id) or ""
    genes, metabolites = set(), set()
    gene_record, metabolite_record = False, False
    for line in data.split('\n'):
        if line.startswith('GENE'): gene_record = True; metabolite_record = False; line_cont = line[12:].strip()
        elif line.startswith('COMPOUND'): metabolite_record = True; gene_record = False; line_cont = line[12:].strip()
        elif (gene_record or metabolite_record) and line.startswith(' '): line_cont = line.strip()
        else: gene_record, metabolite_record = False, False; continue
        if line_cont:
            parts = line_cont.split()
            if parts:
                if gene_record and parts[0].isdigit(): genes.add(parts[1])
                elif metabolite_record: metabolites.add(parts[0])
    return genes, metabolites

def generate_icon_links(name: str, id_type: str) -> dict:
    html = {}
    if id_type == 'gene':
        urls = {
            "Gene Cards": f"https://www.genecards.org/Search/Keyword?queryString={name}",
            "NCBI": f"https://www.ncbi.nlm.nih.gov/gene/?term={name}",
            "ENSEMBL": f"https://useast.ensembl.org/Search/Results?q={name}",
            "GEO": f"https://www.ncbi.nlm.nih.gov/gds/?term={name}",
        }
        for db, url in urls.items():
            b64src = ICON_B64.get(db)
            if b64src: html[db] = (f'<a href="{url}" target="_blank"><img src="{b64src}" width="24" height="24" style="margin:2px; border-radius:4px; border:1px solid #ccc;"></a>')
            else: html[db] = f'<a href="{url}" target="_blank">{db}</a>'
    elif id_type == 'metabolite':
        url = f"https://www.kegg.jp/entry/{name}"
        b64src = ICON_B64.get("KEGG")
        if b64src: html["KEGG"] = (f'<a href="{url}" target="_blank"><img src="{b64src}" width="24" height="24" style="margin:2px; border-radius:4px; border:1px solid #ccc;"></a>')
        else: html["KEGG"] = f'<a href="{url}" target="_blank">KEGG</a>'
    return html

# --- Sidebar ---
all_paths = get_all_kegg_pathways()
st.sidebar.header("Select KEGG Pathways")
show_all = st.sidebar.checkbox("Show all (~350)", value=False)
if show_all:
    options = all_paths
    default = list(options.keys())[:5]
else:
    keywords = ['cancer','vitamin d','diabetes','drug']
    options = {n:pid for n,pid in all_paths.items() if any(k in n.lower() for k in keywords)}
    default = list(options.keys())[:8]
sel = st.sidebar.multiselect("Pathways:", list(options.keys()), default=default, max_selections=15)

# --- Main App Logic ---
if sel:
    with st.spinner("Fetching and analyzing pathway data..."):
        p2data = {n: get_pathway_data(options[n]) for n in sel}
        p2g = {n: data[0] for n, data in p2data.items()}
        p2m = {n: data[1] for n, data in p2data.items()}

    # --- Gene Analysis Section ---
    st.header("📊 Shared Gene Analysis")
    flat_genes = [g for genes in p2g.values() for g in genes]
    gene_counts = Counter(flat_genes)
    genes_by_freq = {cnt: [g for g, c in gene_counts.items() if c == cnt] for cnt in range(2, len(sel) + 1)}
    if any(genes_by_freq.values()):
        for freq in sorted(genes_by_freq.keys(), reverse=True):
            gene_list = genes_by_freq[freq]
            if not gene_list: continue
            st.subheader(f"Genes in {freq} pathways")
            page_data = [{'Gene': g, 'Pathways': ', '.join([n for n, gs in p2g.items() if g in gs]), **generate_icon_links(g, 'gene')} for g in gene_list]
            df_page = pd.DataFrame(page_data).sort_values(by="Gene")
            if len(df_page) <= 20: st.write(f"Showing all {len(df_page)} shared genes:"); display_df = df_page
            else: st.write(f"Showing a preview of the top 20 of {len(df_page)} total genes:"); display_df = df_page.head(20)
            cols = ['Gene','Pathways','Gene Cards','NCBI','ENSEMBL','GEO']
            st.write(display_df[cols].to_html(escape=False,index=False), unsafe_allow_html=True)
            csv_string = df_page.to_csv(index=False).encode('utf-8')
            st.download_button(label=f"Download Full List ({len(df_page)} genes)", data=csv_string, file_name=f"shared_genes_{freq}_pathways.csv", mime="text/csv", key=f"dl_gene_{freq}")
            st.markdown("---")
    else: st.write("⚠️ No shared genes found.")

    # --- Metabolite Analysis Section ---
    st.header("🧪 Shared Metabolite Analysis")
    flat_metabolites = [m for metabolites in p2m.values() for m in metabolites]
    metabolite_counts = Counter(flat_metabolites)
    metabolites_by_freq = {cnt: [m for m, c in metabolite_counts.items() if c == cnt] for cnt in range(2, len(sel) + 1)}
    if any(metabolites_by_freq.values()):
        for freq in sorted(metabolites_by_freq.keys(), reverse=True):
            metabolite_list = metabolites_by_freq[freq]
            if not metabolite_list: continue
            st.subheader(f"Metabolites in {freq} pathways")
            page_data = [{'Metabolite': m, 'Pathways': ', '.join([n for n, ms in p2m.items() if m in ms]), **generate_icon_links(m, 'metabolite')} for m in metabolite_list]
            df_page = pd.DataFrame(page_data).sort_values(by="Metabolite")
            if len(df_page) <= 20: st.write(f"Showing all {len(df_page)} shared metabolites:"); display_df = df_page
            else: st.write(f"Showing a preview of the top 20 of {len(df_page)} total metabolites:"); display_df = df_page.head(20)
            cols = ['Metabolite','Pathways','KEGG']
            st.write(display_df[cols].to_html(escape=False,index=False), unsafe_allow_html=True)
            csv_string = df_page.to_csv(index=False).encode('utf-8')
            st.download_button(label=f"Download Full List ({len(df_page)} metabolites)", data=csv_string, file_name=f"shared_metabolites_{freq}_pathways.csv", mime="text/csv", key=f"dl_metabolite_{freq}")
            st.markdown("---")
    else: st.write("⚠️ No shared metabolites found.")

    # --- Network Visualization ---
    st.header("🕸️ Pathway Network")
    # NEW: Radio buttons to select the network view
    view_choice = st.radio("Select Network View:",
                           ('Genes & Metabolites (Combined)', 'Genes Only', 'Metabolites Only'),
                           horizontal=True)

    G = nx.Graph()
    # Populate graph based on choice
    if view_choice == 'Genes & Metabolites (Combined)':
        for p in p2g: G.add_node(p, type="pathway")
        for g, c in gene_counts.items(): G.add_node(g, type="gene", count=c)
        for m, c in metabolite_counts.items(): G.add_node(m, type="metabolite", count=c)
        for p, gs in p2g.items():
            for g in gs: G.add_edge(p, g)
        for p, ms in p2m.items():
            for m in ms: G.add_edge(p, m)
    elif view_choice == 'Genes Only':
        for p in p2g: G.add_node(p, type="pathway")
        for g, c in gene_counts.items(): G.add_node(g, type="gene", count=c)
        for p, gs in p2g.items():
            for g in gs: G.add_edge(p, g)
    elif view_choice == 'Metabolites Only':
        for p in p2m: G.add_node(p, type="pathway")
        for m, c in metabolite_counts.items(): G.add_node(m, type="metabolite", count=c)
        for p, ms in p2m.items():
            for m in ms: G.add_edge(p, m)
    
    if G.nodes:
        pos = nx.spring_layout(G, k=0.5, seed=42)
        edge_x, edge_y = [], []
        for u, v in G.edges():
            x0, y0 = pos[u]; x1, y1 = pos[v]
            edge_x += [x0, x1, None]; edge_y += [y0, y1, None]
        edge_trace = go.Scatter(x=edge_x, y=edge_y, mode='lines', line=dict(color='#888', width=0.5), hoverinfo='none')

        node_x, node_y, text, color = [], [], [], []
        for node, attr in G.nodes(data=True):
            x, y = pos[node]
            node_x.append(x); node_y.append(y)
            if attr["type"] == "pathway":
                text.append(node); color.append("skyblue")
            elif attr["type"] == "gene":
                count = attr.get("count", 1); text.append(f"{node} (in {count} pathways)"); color.append("red" if count > 1 else "lightcoral")
            elif attr["type"] == "metabolite":
                count = attr.get("count", 1); text.append(f"{node} (in {count} pathways)"); color.append("gold" if count > 1 else "lightyellow")
        
        node_trace = go.Scatter(x=node_x, y=node_y, mode='markers', marker=dict(size=10, color=color, line_width=2), text=text, hoverinfo='text')
        fig = go.Figure(data=[edge_trace, node_trace], layout=go.Layout(title=f"{view_choice} Network", xaxis=dict(showgrid=False, zeroline=False, showticklabels=False), yaxis=dict(showgrid=False, zeroline=False, showticklabels=False), margin=dict(b=20, l=5, r=5, t=40), hovermode='closest'))
        st.plotly_chart(fig, use_container_width=True)
else:
    st.info("☝️ Please select at least one pathway from the sidebar to begin.")