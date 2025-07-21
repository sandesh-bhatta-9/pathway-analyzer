import streamlit as st
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from bioservices import KEGG
from adjustText import adjust_text

st.set_page_config(layout="wide")
st.title('Polished Pathway Overlap Analyzer')
st.write("This tool analyzes pathways from KEGG, provides detailed, linked gene tables, and generates a dynamic network graph.")



@st.cache_data
def get_all_kegg_pathways():
    """Fetches a list of all human pathways from KEGG."""
    k = KEGG()
    pathways_raw = k.list("pathway/hsa")
    pathway_dict = {}
    for line in pathways_raw.strip().split('\n'):
        path_id, name_desc = line.split('\t')
        name = name_desc.split(' - ')[0]
        clean_id = path_id.replace('path:', '')
        pathway_dict[name] = clean_id
    return pathway_dict

@st.cache_data
def get_genes(pathway_id):
    """Fetches the gene set for a specific KEGG pathway."""
    k = KEGG()
    k.organism = "hsa"
    raw_data = k.get(pathway_id)
    if not raw_data: return set()
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
    


def generate_gene_links(gene_name):
    """Generates markdown hyperlinks for a given gene name."""
    base_urls = {
        "GeneCards": f"https://www.genecards.org/Search/Keyword?queryString={gene_name}",
        "NCBI": f"https://www.ncbi.nlm.nih.gov/gene/?term={gene_name}",
        "ENSEMBL": f"https://useast.ensembl.org/Search/Results?q={gene_name}",
        "GEO": f"https://www.ncbi.nlm.nih.gov/gds/?term={gene_name}",
        "UCSC": f"https://genome.ucsc.edu/cgi-bin/hgSearch?search={gene_name}"
    }
    return {db: f"[{db}]({url})" for db, url in base_urls.items()}

# --- Main App Logic ---

all_pathways = get_all_kegg_pathways()
keywords = ['cancer', 'ascorbate', 'vitamin d', 'diabetes', 'drug metabolism', 'xenobiotics', 'fluoropyrimidine']
filtered_pathways = {
    name: pid for name, pid in all_pathways.items()
    if any(keyword in name.lower() for keyword in keywords)
}

st.sidebar.header('Select Pathways for Analysis')
selected_pathways = st.sidebar.multiselect(
    'Choose from pathways related to cancer, vitamins, diabetes, and drugs:',
    options=list(filtered_pathways.keys()),
    default=list(filtered_pathways.keys())[:10]
)

if selected_pathways:
    with st.spinner('Fetching gene data...'):
        pathway_genes_map = {name: get_genes(filtered_pathways[name]) for name in selected_pathways}

    st.header('Shared Gene Analysis')
    all_genes_flat = [gene for sublist in pathway_genes_map.values() for gene in sublist]
    gene_counts = Counter(all_genes_flat)
    
    
    shared_gene_data = []
    for gene, frequency in gene_counts.items():
        if frequency >= 2:
            member_of_pathways = [name for name, genes in pathway_genes_map.items() if gene in genes]
            links = generate_gene_links(gene)
            gene_info = {'Gene': gene, 'Frequency': frequency, 'Pathways': ', '.join(member_of_pathways)}
            gene_info.update(links) 
            shared_gene_data.append(gene_info)

    if shared_gene_data:
        df_shared_enhanced = pd.DataFrame(shared_gene_data)
        df_shared_enhanced = df_shared_enhanced.sort_values(by='Frequency', ascending=False)

       
        grouped = df_shared_enhanced.groupby('Frequency')
        for freq, group in sorted(grouped, reverse=True):
            st.subheader(f"Genes Shared in {freq} Pathways")
            display_cols = ['Gene', 'Pathways', 'GeneCards', 'NCBI', 'ENSEMBL', 'GEO', 'UCSC']
            st.write(group[display_cols].to_html(escape=False, index=False), unsafe_allow_html=True)
    else:
        st.write("No shared genes found for the selected pathways.")
    

    st.header('Network Visualization')
    with st.spinner('Generating network graph...'):
        G = nx.Graph()
        
        def get_pathway_category(name):
            name_lower = name.lower()
            if 'cancer' in name_lower: return 'Cancer', 'forestgreen'
            if 'vitamin' in name_lower or 'ascorbate' in name_lower: return 'Vitamin', 'orchid'
            if 'diabetes' in name_lower: return 'Diabetes', 'darkorange'
            if any(k in name_lower for k in ['drug', 'xenobiotics', 'fluoropyrimidine']): return 'Drug Metabolism', 'royalblue'
            return 'Other', 'gray'

        node_colors, pathway_categories = {}, set()
        for name in pathway_genes_map.keys():
            category, color = get_pathway_category(name)
            G.add_node(name, type='pathway', category=category)
            node_colors[name] = color
            pathway_categories.add((category, color))

        for gene, count in gene_counts.items():
            G.add_node(gene, type='gene')
            node_colors[gene] = 'red' if count > 1 else 'lightgray'

        for name, genes in pathway_genes_map.items():
            for gene in genes:
                G.add_edge(name, gene)
        
        fig, ax = plt.subplots(figsize=(25, 20))
        pos = nx.spring_layout(G, k=0.6, iterations=50, seed=42)
        nx.draw(G, pos, ax=ax, with_labels=False, node_color=[node_colors[n] for n in G.nodes()], node_size=90, edge_color='#cccccc')

        pathway_labels = {n: n for n in G.nodes if G.nodes[n]['type'] == 'pathway'}
        texts = [ax.text(pos[n][0], pos[n][1], text, fontsize=12, fontweight='bold') for n, text in pathway_labels.items()]
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

        
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=f'{cat} Pathway', markerfacecolor=color, markersize=15) for cat, color in sorted(list(pathway_categories))]
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', label='Shared Gene', markerfacecolor='red', markersize=15))
        ax.legend(handles=legend_elements, loc='upper right', fontsize=16)

        st.pyplot(fig)
else:
    st.info("Please select pathways from the sidebar to begin the analysis.")