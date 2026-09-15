import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from collections import Counter
from bioservices import KEGG
import concurrent.futures
import base64
from pathlib import Path


# ---------------------------------------------------------
# 1. PAGE CONFIGURATION
# ---------------------------------------------------------

st.set_page_config(
    page_title="Pathway Analyzer",
    layout="wide"
)

st.title("🧬 Pathway Analyzer")
st.write(
    "Explore a Super Pathway connecting cancer, natural product, "
    "and pharmaceutical pathways through shared genes."
)


# ---------------------------------------------------------
# 2. ICONS
# ---------------------------------------------------------

ICON_NAMES = [
    "Gene Cards",
    "NCBI",
    "ENSEMBL",
    "KEGG",
    "GEO"
]

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


# ---------------------------------------------------------
# 3. KEGG FUNCTIONS
# ---------------------------------------------------------

@st.cache_data
def get_all_kegg_pathways():

    k = KEGG()

    pathways_raw = k.list("pathway/hsa")

    pathways = {}

    for line in pathways_raw.strip().split("\n"):

        try:

            pid, name_desc = line.split("\t")

            name = name_desc.split(" - ")[0]

            pathways[name] = pid.replace("path:", "")

        except ValueError:
            continue

    return pathways


@st.cache_data
def get_genes_from_pathway(pathway_id):

    k = KEGG()

    genes = set()

    try:

        data = k.get(pathway_id)

        if not data:
            return set()

        parsed = k.parse(data)

        if "GENE" not in parsed:
            return set()

        gene_data = parsed["GENE"]

        if isinstance(gene_data, dict):

            for description in gene_data.values():

                symbol = description.split(";")[0].strip()

                if symbol:
                    genes.add(symbol)

        elif isinstance(gene_data, list):

            for entry in gene_data:

                parts = entry.split()

                if len(parts) > 1:

                    symbol = parts[1].replace(";", "").strip()

                    if symbol:
                        genes.add(symbol)

        return genes

    except Exception as e:

        st.error(
            f"Error retrieving pathway {pathway_id}: {e}"
        )

        return set()


# ---------------------------------------------------------
# 4. GENE LINKS
# ---------------------------------------------------------

def generate_url_links(gene_name):

    return {

        "GeneCards_URL":
            f"https://www.genecards.org/Search/Keyword?queryString={gene_name}",

        "NCBI_URL":
            f"https://www.ncbi.nlm.nih.gov/gene/?term={gene_name}",

        "ENSEMBL_URL":
            f"https://useast.ensembl.org/Search/Results?q={gene_name}",

        "GEO_URL":
            f"https://www.ncbi.nlm.nih.gov/gds/?term={gene_name}"

    }


def generate_icon_links(gene_name):

    urls = {

        "Gene Cards":
            f"https://www.genecards.org/Search/Keyword?queryString={gene_name}",

        "NCBI":
            f"https://www.ncbi.nlm.nih.gov/gene/?term={gene_name}",

        "ENSEMBL":
            f"https://useast.ensembl.org/Search/Results?q={gene_name}",

        "GEO":
            f"https://www.ncbi.nlm.nih.gov/gds/?term={gene_name}"

    }

    html = {}

    for db, url in urls.items():

        b64src = ICON_B64.get(db)

        if b64src:

            html[db] = (
                f'<a href="{url}" target="_blank">'
                f'<img src="{b64src}" width="24" height="24" '
                f'style="margin:2px; border-radius:4px; '
                f'border:1px solid #ccc;">'
                f'</a>'
            )

        else:

            html[db] = (
                f'<a href="{url}" target="_blank">{db}</a>'
            )

    return html


# ---------------------------------------------------------
# 5. PATHWAY CATEGORIES
# ---------------------------------------------------------

all_paths = get_all_kegg_pathways()


CANCER_KEYWORDS = [
    "cancer",
    "glioma",
    "leukemia",
    "lymphoma",
    "melanoma",
    "gastric",
    "colorectal",
    "prostate",
    "breast",
    "renal cell"
]


CHEMO_KEYWORDS = [
    "fluoropyrimidine",
    "folate",
    "platinum"
]


NATURAL_PRODUCT_KEYWORDS = [
    "metabolism of xenobiotics by cytochrome p450",
    "drug metabolism - cytochrome p450",
    "drug metabolism - other enzymes",
    "steroid hormone biosynthesis",
    "retinol metabolism",
    "metabolism"
]


cancer_options = {
    name: pid
    for name, pid in all_paths.items()
    if any(
        keyword in name.lower()
        for keyword in CANCER_KEYWORDS
    )
}


chemo_options = {
    name: pid
    for name, pid in all_paths.items()
    if any(
        keyword in name.lower()
        for keyword in CHEMO_KEYWORDS
    )
}


natural_options = {
    name: pid
    for name, pid in all_paths.items()
    if any(
        keyword in name.lower()
        for keyword in NATURAL_PRODUCT_KEYWORDS
    )
}


# ---------------------------------------------------------
# 6. SIDEBAR
# ---------------------------------------------------------

st.sidebar.header("Select Pathways")


st.sidebar.subheader("Cancer Pathways")

sel_cancer = st.sidebar.multiselect(
    "Select cancer pathways:",
    list(cancer_options.keys())
)


st.sidebar.subheader("Natural Product Pathways")

sel_natural = st.sidebar.multiselect(
    "Select natural product pathways:",
    list(natural_options.keys())
)


st.sidebar.subheader("Pharmaceutical Pathways")

sel_chemo = st.sidebar.multiselect(
    "Select pharmaceutical pathways:",
    list(chemo_options.keys())
)


selected_pathways = (
    sel_cancer +
    sel_natural +
    sel_chemo
)


# ---------------------------------------------------------
# 7. MAIN ANALYSIS
# ---------------------------------------------------------

if selected_pathways:

    # -----------------------------------------------------
    # Fetch pathway genes
    # -----------------------------------------------------

    pathway_ids = {}

    pathway_type = {}

    for pathway in sel_cancer:

        pathway_ids[pathway] = cancer_options[pathway]
        pathway_type[pathway] = "cancer"


    for pathway in sel_natural:

        pathway_ids[pathway] = natural_options[pathway]
        pathway_type[pathway] = "natural"


    for pathway in sel_chemo:

        pathway_ids[pathway] = chemo_options[pathway]
        pathway_type[pathway] = "pharma"


    p2g = {}

    with st.spinner(
        f"Fetching {len(selected_pathways)} pathways..."
    ):

        with concurrent.futures.ThreadPoolExecutor() as executor:

            future_to_pathway = {

                executor.submit(
                    get_genes_from_pathway,
                    pathway_ids[pathway]
                ): pathway

                for pathway in selected_pathways
            }


            for future in concurrent.futures.as_completed(
                future_to_pathway
            ):

                pathway = future_to_pathway[future]

                try:
                    p2g[pathway] = future.result()

                except Exception:
                    p2g[pathway] = set()


    # -----------------------------------------------------
    # Gene mapping
    # -----------------------------------------------------

    flat_genes = [
        gene
        for genes in p2g.values()
        for gene in genes
    ]

    gene_counts = Counter(flat_genes)


    # Gene -> Cancer pathways

    gene_to_cancers = {}

    for gene in gene_counts:

        gene_to_cancers[gene] = set()

        for cancer in sel_cancer:

            if gene in p2g.get(cancer, set()):

                gene_to_cancers[gene].add(cancer)


    # Gene -> Natural product pathways

    gene_to_natural = {}

    for gene in gene_counts:

        gene_to_natural[gene] = set()

        for pathway in sel_natural:

            if gene in p2g.get(pathway, set()):

                gene_to_natural[gene].add(pathway)


    # Gene -> Pharmaceutical pathways

    gene_to_pharma = {}

    for gene in gene_counts:

        gene_to_pharma[gene] = set()

        for pathway in sel_chemo:

            if gene in p2g.get(pathway, set()):

                gene_to_pharma[gene].add(pathway)


    # -----------------------------------------------------
    # 8. SUPER PATHWAY
    # -----------------------------------------------------

    st.header("🕸️ Super Pathway")

    st.write(
        "The Super Pathway combines the selected cancer, "
        "natural product, and pharmaceutical pathways through "
        "their shared genes."
    )


    G = nx.Graph()


    # Add pathway nodes

    for pathway in selected_pathways:

        G.add_node(
            pathway,
            type="pathway",
            pathway_type=pathway_type[pathway]
        )


    # Add gene nodes

    for gene, count in gene_counts.items():

        G.add_node(
            gene,
            type="gene",
            count=count
        )


    # Add pathway-gene edges

    for pathway, genes in p2g.items():

        for gene in genes:

            G.add_edge(
                pathway,
                gene
            )


    # -----------------------------------------------------
    # Cancer colors
    # -----------------------------------------------------

    cancer_colors = [

        "crimson",
        "royalblue",
        "forestgreen",
        "darkorange",
        "purple",
        "deeppink",
        "brown",
        "teal",
        "goldenrod",
        "darkcyan"

    ]


    cancer_color_map = {}

    for i, cancer in enumerate(sel_cancer):

        cancer_color_map[cancer] = (
            cancer_colors[
                i % len(cancer_colors)
            ]
        )


    # -----------------------------------------------------
    # Network layout
    # -----------------------------------------------------

    pos = nx.spring_layout(
        G,
        k=0.7,
        seed=42
    )


    # -----------------------------------------------------
    # Edges
    # -----------------------------------------------------

    edge_x = []
    edge_y = []

    for u, v in G.edges():

        x0, y0 = pos[u]
        x1, y1 = pos[v]

        edge_x.extend(
            [x0, x1, None]
        )

        edge_y.extend(
            [y0, y1, None]
        )


    edge_trace = go.Scatter(

        x=edge_x,
        y=edge_y,

        mode="lines",

        line=dict(
            color="#999",
            width=0.6
        ),

        hoverinfo="none"

    )


    # -----------------------------------------------------
    # Pathway nodes and gene nodes
    # -----------------------------------------------------

    pathway_x = []
    pathway_y = []
    pathway_text = []
    pathway_color = []
    pathway_symbol = []


    gene_x = []
    gene_y = []
    gene_text = []
    gene_color = []


    # -----------------------------------------------------
    # Pathway nodes
    # -----------------------------------------------------

    for pathway in selected_pathways:

        x, y = pos[pathway]

        pathway_x.append(x)
        pathway_y.append(y)

        ptype = pathway_type[pathway]


        if ptype == "cancer":

            pathway_color.append(
                cancer_color_map[pathway]
            )

            pathway_symbol.append(
                "circle"
            )

            pathway_text.append(
                f"{pathway}<br>"
                f"Type: Cancer pathway"
            )


        elif ptype == "natural":

            pathway_color.append(
                "forestgreen"
            )

            pathway_symbol.append(
                "square"
            )

            pathway_text.append(
                f"{pathway}<br>"
                f"Type: Natural product pathway"
            )


        elif ptype == "pharma":

            pathway_color.append(
                "mediumblue"
            )

            pathway_symbol.append(
                "square"
            )

            pathway_text.append(
                f"{pathway}<br>"
                f"Type: Pharmaceutical pathway"
            )


    # -----------------------------------------------------
    # Gene nodes
    # -----------------------------------------------------

    total_cancers = len(sel_cancer)


    for gene in gene_counts:

        x, y = pos[gene]

        gene_x.append(x)
        gene_y.append(y)


        cancers = gene_to_cancers.get(
            gene,
            set()
        )

        natural_paths = gene_to_natural.get(
            gene,
            set()
        )

        pharma_paths = gene_to_pharma.get(
            gene,
            set()
        )


        cancer_names = ", ".join(
            sorted(cancers)
        )

        natural_names = ", ".join(
            sorted(natural_paths)
        )

        pharma_names = ", ".join(
            sorted(pharma_paths)
        )


        # Gene coloring based on cancer membership

        if (
            total_cancers > 1
            and len(cancers) == total_cancers
        ):

            gene_color.append(
                "hotpink"
            )

            cancer_status = (
                f"Shared by ALL {total_cancers} cancers"
            )


        elif len(cancers) > 1:

            gene_color.append(
                "cyan"
            )

            cancer_status = (
                f"Shared by {len(cancers)} cancers"
            )


        elif len(cancers) == 1:

            cancer_color = cancer_color_map[
                list(cancers)[0]
            ]

            gene_color.append(
                cancer_color
            )

            cancer_status = (
                f"Unique to {list(cancers)[0]}"
            )


        else:

            gene_color.append(
                "lightgray"
            )

            cancer_status = (
                "Not present in selected cancer pathways"
            )


        # Hover information

        hover_text = (
            f"<b>{gene}</b><br>"
            f"{cancer_status}<br><br>"
            f"<b>Cancer pathways:</b><br>"
            f"{cancer_names if cancer_names else 'None'}<br><br>"
            f"<b>Natural product pathways:</b><br>"
            f"{natural_names if natural_names else 'None'}<br><br>"
            f"<b>Pharmaceutical pathways:</b><br>"
            f"{pharma_names if pharma_names else 'None'}"
        )


        gene_text.append(
            hover_text
        )


    # -----------------------------------------------------
    # Pathway trace
    # -----------------------------------------------------

    pathway_trace = go.Scatter(

        x=pathway_x,
        y=pathway_y,

        mode="markers+text",

        marker=dict(

            size=24,

            color=pathway_color,

            symbol=pathway_symbol,

            line=dict(
                width=2,
                color="black"
            )

        ),

        text=[
            pathway
            for pathway in selected_pathways
        ],

        textposition="top center",

        hovertext=pathway_text,

        hoverinfo="text"

    )


    # -----------------------------------------------------
    # Gene trace
    # -----------------------------------------------------

    gene_trace = go.Scatter(

        x=gene_x,
        y=gene_y,

        mode="markers",

        marker=dict(

            size=11,

            color=gene_color,

            line=dict(
                width=1,
                color="black"
            )

        ),

        text=gene_text,

        hoverinfo="text"

    )


    # -----------------------------------------------------
    # Figure
    # -----------------------------------------------------

    fig = go.Figure(

        data=[
            edge_trace,
            pathway_trace,
            gene_trace
        ]

    )


    fig.update_layout(

        title=(
            "Super Pathway: Cancer, Natural Product "
            "and Pharmaceutical Pathways"
        ),

        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False
        ),

        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False
        ),

        margin=dict(
            b=20,
            l=5,
            r=5,
            t=60
        ),

        hovermode="closest",

        height=750,

        showlegend=False

    )


    st.plotly_chart(
        fig,
        use_container_width=True
    )


    # -----------------------------------------------------
    # 9. LEGEND
    # -----------------------------------------------------

    st.subheader("Legend")

    st.markdown(
        """
        **Pathway nodes**

        ◯ Cancer pathways — each cancer has its own color

        ■ Natural product pathways — green square

        ■ Pharmaceutical pathways — blue square

        **Gene nodes**

        🩷 Shared by all selected cancers

        🔵 Shared by some selected cancers

        Colored gene — associated with one selected cancer

        ⚪ Gene not present in selected cancer pathways
        """
    )


    # -----------------------------------------------------
    # 10. GENE ANALYSIS
    # -----------------------------------------------------

    st.markdown("---")

    st.header("📊 Gene Analysis")


    gene_filter = st.text_input(
        "Search for a gene (for example: TP53)"
    ).strip().upper()


    genes_all_cancers = sorted(

        [
            gene
            for gene, cancers
            in gene_to_cancers.items()

            if total_cancers > 1
            and len(cancers) == total_cancers
        ]

    )


    genes_some_cancers = sorted(

        [
            gene
            for gene, cancers
            in gene_to_cancers.items()

            if len(cancers) > 1
            and len(cancers) < total_cancers
        ]

    )


    genes_one_cancer = sorted(

        [
            gene
            for gene, cancers
            in gene_to_cancers.items()

            if len(cancers) == 1
        ]

    )


    tabs = st.tabs(

        [
            "Shared by ALL Cancers",
            "Shared by SOME Cancers",
            "Unique to ONE Cancer"
        ]

    )


    # -----------------------------------------------------
    # ALL CANCERS
    # -----------------------------------------------------

    with tabs[0]:

        if total_cancers < 2:

            st.info(
                "Select at least two cancer pathways "
                "to compare shared cancer genes."
            )

        elif not genes_all_cancers:

            st.warning(
                "No genes were found to be common "
                "to all selected cancer pathways."
            )

        else:

            st.write(
                f"Found {len(genes_all_cancers)} genes "
                f"common to all selected cancers."
            )


            display_data = []


            for gene in genes_all_cancers:

                row = {

                    "Gene": gene,

                    "Cancer Pathways":
                        ", ".join(
                            sorted(
                                gene_to_cancers[gene]
                            )
                        ),

                    "Natural Product Pathways":
                        ", ".join(
                            sorted(
                                gene_to_natural.get(
                                    gene,
                                    set()
                                )
                            )
                        ),

                    "Pharmaceutical Pathways":
                        ", ".join(
                            sorted(
                                gene_to_pharma.get(
                                    gene,
                                    set()
                                )
                            )
                        ),

                    **generate_icon_links(gene)

                }

                display_data.append(row)


            df = pd.DataFrame(
                display_data
            )


            if gene_filter:

                df = df[
                    df["Gene"]
                    .str.upper()
                    .str.contains(gene_filter)
                ]


            columns = [

                "Gene",
                "Cancer Pathways",
                "Natural Product Pathways",
                "Pharmaceutical Pathways",
                "Gene Cards",
                "NCBI",
                "ENSEMBL",
                "GEO"

            ]


            st.write(
                df[columns]
                .to_html(
                    escape=False,
                    index=False
                ),

                unsafe_allow_html=True
            )


    # -----------------------------------------------------
    # SOME CANCERS
    # -----------------------------------------------------

    with tabs[1]:

        if not genes_some_cancers:

            st.info(
                "No genes were found to be shared "
                "by some but not all selected cancers."
            )

        else:

            display_data = []


            for gene in genes_some_cancers:

                row = {

                    "Gene": gene,

                    "Cancer Pathways":
                        ", ".join(
                            sorted(
                                gene_to_cancers[gene]
                            )
                        ),

                    "Natural Product Pathways":
                        ", ".join(
                            sorted(
                                gene_to_natural.get(
                                    gene,
                                    set()
                                )
                            )
                        ),

                    "Pharmaceutical Pathways":
                        ", ".join(
                            sorted(
                                gene_to_pharma.get(
                                    gene,
                                    set()
                                )
                            )
                        ),

                    **generate_icon_links(gene)

                }

                display_data.append(row)


            df = pd.DataFrame(
                display_data
            )


            if gene_filter:

                df = df[
                    df["Gene"]
                    .str.upper()
                    .str.contains(gene_filter)
                ]


            columns = [

                "Gene",
                "Cancer Pathways",
                "Natural Product Pathways",
                "Pharmaceutical Pathways",
                "Gene Cards",
                "NCBI",
                "ENSEMBL",
                "GEO"

            ]


            st.write(
                df[columns]
                .to_html(
                    escape=False,
                    index=False
                ),

                unsafe_allow_html=True
            )


    # -----------------------------------------------------
    # ONE CANCER
    # -----------------------------------------------------

    with tabs[2]:

        if not genes_one_cancer:

            st.info(
                "No genes were found to be unique "
                "to a single selected cancer."
            )

        else:

            display_data = []


            for gene in genes_one_cancer:

                row = {

                    "Gene": gene,

                    "Cancer Pathway":
                        list(
                            gene_to_cancers[gene]
                        )[0],

                    "Natural Product Pathways":
                        ", ".join(
                            sorted(
                                gene_to_natural.get(
                                    gene,
                                    set()
                                )
                            )
                        ),

                    "Pharmaceutical Pathways":
                        ", ".join(
                            sorted(
                                gene_to_pharma.get(
                                    gene,
                                    set()
                                )
                            )
                        ),

                    **generate_icon_links(gene)

                }

                display_data.append(row)


            df = pd.DataFrame(
                display_data
            )


            if gene_filter:

                df = df[
                    df["Gene"]
                    .str.upper()
                    .str.contains(gene_filter)
                ]


            columns = [

                "Gene",
                "Cancer Pathway",
                "Natural Product Pathways",
                "Pharmaceutical Pathways",
                "Gene Cards",
                "NCBI",
                "ENSEMBL",
                "GEO"

            ]


            st.write(
                df[columns]
                .to_html(
                    escape=False,
                    index=False
                ),

                unsafe_allow_html=True
            )


    # -----------------------------------------------------
    # 11. DOWNLOAD DATA
    # -----------------------------------------------------

    st.markdown("---")

    st.subheader("📥 Download Gene Data")


    download_data = []


    for gene in sorted(gene_counts):

        row = {

            "Gene Name": gene,

            "Cancer Pathways":
                ", ".join(
                    sorted(
                        gene_to_cancers.get(
                            gene,
                            set()
                        )
                    )
                ),

            "Natural Product Pathways":
                ", ".join(
                    sorted(
                        gene_to_natural.get(
                            gene,
                            set()
                        )
                    )
                ),

            "Pharmaceutical Pathways":
                ", ".join(
                    sorted(
                        gene_to_pharma.get(
                            gene,
                            set()
                        )
                    )
                ),

            "Frequency":
                gene_counts[gene],

            **generate_url_links(gene)

        }

        download_data.append(row)


    df_download = pd.DataFrame(
        download_data
    )


    csv_string = df_download.to_csv(
        index=False
    ).encode("utf-8")


    st.download_button(

        label="Download Gene Data",

        data=csv_string,

        file_name="super_pathway_gene_data.csv",

        mime="text/csv"

    )


else:

    st.info(
        "☝️ Select at least one pathway from the sidebar "
        "to build the Super Pathway."
    )
