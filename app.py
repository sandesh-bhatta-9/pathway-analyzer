import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from collections import Counter
from bioservices import KEGG
import concurrent.futures
import base64
from pathlib import Path


# -----------------------------------------------------
# 1. PAGE CONFIGURATION
# -----------------------------------------------------

st.set_page_config(
    page_title="OncoNet Explorer",
    layout="wide"
)

st.title("OncoNet Explorer")
st.write(
    "Integrated visualization of cancer, natural product, "
    "and pharmaceutical pathway-gene associations."
)


# -----------------------------------------------------
# 2. ICON & KEGG UTILITIES
# -----------------------------------------------------

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


@st.cache_data
def get_all_kegg_pathways():
    k = KEGG()

    pathways_raw = k.list("pathway/hsa")

    pathways = {}

    for line in pathways_raw.strip().split("\n"):
        parts = line.split("\t")

        if len(parts) != 2:
            continue

        pid, name_desc = parts

        name = name_desc.split(" - ")[0]

        pathways[name] = pid.replace("path:", "")

    return pathways


@st.cache_data
def get_genes_from_pathway(pathway_id: str):

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
            f"Error parsing pathway {pathway_id}: {e}"
        )

        return set()


def generate_url_links(gene_name: str):

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


def generate_icon_links(gene_name: str):

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


# -----------------------------------------------------
# 3. SIDEBAR
# -----------------------------------------------------

all_paths = get_all_kegg_pathways()

st.sidebar.header("Select Pathways for Comparison")


CANCER_KEYWORDS = [
    "cancer",
    "glioma",
    "leukemia",
    "lymphoma",
    "melanoma",
    "gastric",
    "colorectal",
    "prostate",
    "breast"
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


# -----------------------------------------------------
# 4. EXAMPLE SELECTION
# -----------------------------------------------------

st.sidebar.markdown("---")

st.sidebar.subheader("Load Examples")


def set_example_state(simple=True):

    if simple:

        st.session_state.cancer_key = [
            k
            for k in cancer_options.keys()
            if "Melanoma" in k
        ]

        st.session_state.chemo_key = [
            k
            for k in chemo_options.keys()
            if "Platinum" in k
        ]

        st.session_state.natural_key = []

    else:

        st.session_state.cancer_key = [
            k
            for k in cancer_options.keys()
            if k in [
                "Melanoma",
                "Renal cell carcinoma",
                "Gastric cancer"
            ]
        ]

        st.session_state.chemo_key = [
            k
            for k in chemo_options.keys()
            if "Platinum" in k
        ]

        st.session_state.natural_key = [
            k
            for k in natural_options.keys()
            if "cytochrome P450" in k
        ]


if "cancer_key" not in st.session_state:

    st.session_state.cancer_key = [
        k
        for k in cancer_options.keys()
        if k in [
            "Melanoma",
            "Renal cell carcinoma"
        ]
    ]


if "chemo_key" not in st.session_state:

    st.session_state.chemo_key = []


if "natural_key" not in st.session_state:

    st.session_state.natural_key = []


st.sidebar.button(
    "Load Simple Example",
    on_click=set_example_state,
    args=(True,)
)


st.sidebar.button(
    "Load Complex Example",
    on_click=set_example_state,
    args=(False,)
)


st.sidebar.markdown("---")


# -----------------------------------------------------
# 5. PATHWAY SELECTION
# -----------------------------------------------------

st.sidebar.subheader("Cancer Pathways")

sel_cancer = st.sidebar.multiselect(
    "Select cancer pathways:",
    list(cancer_options.keys()),
    key="cancer_key"
)


st.sidebar.subheader("Therapeutic Pathways")

sel_chemo = st.sidebar.multiselect(
    "Select pharmaceutical pathways:",
    list(chemo_options.keys()),
    key="chemo_key"
)


sel_natural = st.sidebar.multiselect(
    "Select natural product pathways:",
    list(natural_options.keys()),
    key="natural_key"
)


combined_options = {
    **cancer_options,
    **chemo_options,
    **natural_options
}


selected_pathways = (
    sel_cancer +
    sel_chemo +
    sel_natural
)


# -----------------------------------------------------
# 6. MAIN APP LOGIC
# -----------------------------------------------------

if selected_pathways:

    # -------------------------------------------------
    # 6A. FETCH PATHWAY DATA
    # -------------------------------------------------

    p2g = {}

    with st.spinner(
        f"Fetching {len(selected_pathways)} pathways..."
    ):

        with concurrent.futures.ThreadPoolExecutor() as executor:

            future_to_pathway = {

                executor.submit(
                    get_genes_from_pathway,
                    combined_options[pathway]
                ): pathway

                for pathway in selected_pathways
            }

            for future in concurrent.futures.as_completed(
                future_to_pathway
            ):

                pathway = future_to_pathway[future]

                p2g[pathway] = future.result()


    # -------------------------------------------------
    # 6B. GENE MAPPING
    # -------------------------------------------------

    flat_genes = [
        gene
        for genes in p2g.values()
        for gene in genes
    ]

    gene_counts = Counter(flat_genes)


    g2c = {}

    for gene in gene_counts:

        g2c[gene] = set()

        for cancer_name in sel_cancer:

            if gene in p2g.get(cancer_name, set()):

                g2c[gene].add(cancer_name)


    # -------------------------------------------------
    # 6C. PATHWAY NETWORK
    # -------------------------------------------------

    st.header("Pathway-Gene Network")

    st.info(
        "The network represents pathway-gene associations among the "
        "selected cancer, natural product, and pharmaceutical pathways. "
        "Node position is determined by the network layout and should "
        "not be interpreted as a quantitative measure of biological similarity."
    )


    G = nx.Graph()


    for pathway in p2g:

        if pathway in sel_cancer:

            pathway_type = "cancer"

        elif pathway in sel_natural:

            pathway_type = "natural"

        elif pathway in sel_chemo:

            pathway_type = "pharmaceutical"

        else:

            pathway_type = "pathway"


        G.add_node(
            pathway,
            type=pathway_type
        )


    for gene, count in gene_counts.items():

        G.add_node(
            gene,
            type="gene",
            count=count
        )


    for pathway, genes in p2g.items():

        for gene in genes:

            G.add_edge(
                pathway,
                gene
            )


    pos = nx.spring_layout(
        G,
        k=0.5,
        seed=42
    )


    # -------------------------------------------------
    # 6D. NETWORK EDGES
    # -------------------------------------------------

    edge_x = []
    edge_y = []


    for u, v in G.edges():

        x0, y0 = pos[u]
        x1, y1 = pos[v]

        edge_x.extend([
            x0,
            x1,
            None
        ])

        edge_y.extend([
            y0,
            y1,
            None
        ])


    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(
            color="#888",
            width=0.6
        ),
        hoverinfo="none"
    )


    # -------------------------------------------------
    # 6E. CANCER COLORS
    # -------------------------------------------------

    cancer_colors = [

        "#D62728",
        "#1F77B4",
        "#2CA02C",
        "#9467BD",
        "#FF7F0E",
        "#17BECF",
        "#8C564B",
        "#E377C2"
    ]


    cancer_color_map = {

        cancer_name:
            cancer_colors[index % len(cancer_colors)]

        for index, cancer_name
        in enumerate(sel_cancer)
    }


    # -------------------------------------------------
    # 6F. PATHWAY NODE TRACES
    # -------------------------------------------------

    pathway_traces = []


    for pathway in p2g:

        node_type = G.nodes[pathway]["type"]

        x = pos[pathway][0]
        y = pos[pathway][1]


        if node_type == "cancer":

            node_color = cancer_color_map[pathway]
            node_symbol = "circle"

        elif node_type == "natural":

            node_color = "#228B22"
            node_symbol = "square"

        elif node_type == "pharmaceutical":

            node_color = "#0000CD"
            node_symbol = "square"

        else:

            node_color = "gray"
            node_symbol = "circle"


        pathway_trace = go.Scatter(

            x=[x],
            y=[y],

            mode="markers+text",

            text=[pathway],

            textposition="top center",

            marker=dict(

                size=18,

                color=node_color,

                symbol=node_symbol,

                line=dict(
                    width=1.5,
                    color="black"
                )
            ),

            hovertemplate=(
                "<b>%{text}</b>"
                "<br>Pathway type: "
                + node_type.capitalize()
                + "<extra></extra>"
            ),

            showlegend=False
        )


        pathway_traces.append(
            pathway_trace
        )


    # -------------------------------------------------
    # 6G. GENE NODE CATEGORIES
    # -------------------------------------------------

    gene_x = []
    gene_y = []

    gene_text = []
    gene_color = []


    total_selected_cancers = len(
        sel_cancer
    )


    for gene in gene_counts:

        x = pos[gene][0]
        y = pos[gene][1]

        gene_x.append(x)
        gene_y.append(y)


        cancers_for_gene = g2c.get(
            gene,
            set()
        )

        num_cancers = len(
            cancers_for_gene
        )


        if (
            total_selected_cancers > 1
            and num_cancers == total_selected_cancers
        ):

            color = "hotpink"

            category = (
                "Present in all selected cancer pathways"
            )


        elif (
            1 < num_cancers
            < total_selected_cancers
        ):

            color = "cyan"

            category = (
                "Present in multiple selected cancer pathways"
            )


        elif num_cancers == 1:

            cancer_name = list(
                cancers_for_gene
            )[0]

            color = cancer_color_map[
                cancer_name
            ]

            category = (
                "Present in one selected cancer pathway: "
                + cancer_name
            )


        else:

            color = "lightgray"

            category = (
                "Not present in the selected cancer pathways"
            )


        gene_color.append(
            color
        )


        gene_text.append(
            f"{gene}<br>{category}"
        )


    gene_trace = go.Scatter(

        x=gene_x,
        y=gene_y,

        mode="markers",

        marker=dict(

            size=10,

            color=gene_color,

            line=dict(
                width=1,
                color="black"
            )
        ),

        text=gene_text,

        hoverinfo="text",

        showlegend=False
    )


    # -------------------------------------------------
    # 6H. FIGURE
    # -------------------------------------------------

    fig = go.Figure()


    fig.add_trace(
        edge_trace
    )


    for trace in pathway_traces:

        fig.add_trace(
            trace
        )


    fig.add_trace(
        gene_trace
    )


    fig.update_layout(

        title=(
            "Integrated Cancer, Natural Product, "
            "and Pharmaceutical Pathway-Gene Network"
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

        plot_bgcolor="white",

        height=750
    )


    st.plotly_chart(
        fig,
        use_container_width=True
    )


    st.markdown("---")


    # -------------------------------------------------
    # 7. LEGEND
    # -------------------------------------------------

    st.subheader("Legend")


    st.markdown(
        """
        ### Pathway Nodes

        **Cancer pathways**  
        Circular nodes represent selected cancer pathways. Each cancer
        pathway is assigned a distinct color, which is also used to
        identify genes specific to that cancer pathway.

        **Natural product pathways**  
        Green square nodes represent pathways associated with natural products.

        **Pharmaceutical pathways**  
        Blue square nodes represent pharmaceutical or chemotherapy-related pathways.


        ### Gene Nodes

        **Pink**  
        Gene is present in **all selected cancer pathways**.

        **Cyan**  
        Gene is present in **more than one, but not all, selected cancer pathways**.

        **Cancer-specific color**  
        Gene is present in **only one selected cancer pathway**.
        The gene color corresponds to the color assigned to that
        cancer pathway.

        **Light gray**  
        Gene is **not present in any of the selected cancer pathways**.


        ### Network Interpretation

        **Pathway–gene edge**  
        An edge connecting a pathway and a gene indicates that the gene
        is associated with the corresponding pathway.

        **Shared genes**  
        Genes connected to multiple selected cancer pathways represent
        common molecular components across those pathways.

        **Cancer-specific genes**  
        Genes connected to only one selected cancer pathway represent
        pathway-specific gene associations within the selected analysis.

        **Natural product and pharmaceutical connections**  
        Connections between these pathway nodes and genes provide a
        visual representation of their relationships within the
        integrated pathway-gene network.
        """
    )


    st.markdown("---")


    # -------------------------------------------------
    # 8. SHARED GENE ANALYSIS
    # -------------------------------------------------

    st.header("Shared Gene Analysis")


    gene_filter = st.text_input(
        "Search for a specific gene in the tables below "
        "(e.g., TP53)"
    ).strip().upper()


    genes_all_cancers = sorted(
        [
            gene
            for gene, cancers in g2c.items()
            if len(cancers) == total_selected_cancers
            and total_selected_cancers > 1
        ]
    )


    genes_some_cancers = sorted(
        [
            gene
            for gene, cancers in g2c.items()
            if 1 < len(cancers) < total_selected_cancers
        ]
    )


    genes_one_cancer = sorted(
        [
            gene
            for gene, cancers in g2c.items()
            if len(cancers) == 1
        ]
    )


    analysis_tabs = st.tabs(
        [
            "Genes in All Selected Cancer Pathways",
            "Genes in Multiple Selected Cancer Pathways",
            "Genes in One Selected Cancer Pathway"
        ]
    )


    # -------------------------------------------------
    # 8A. ALL CANCER PATHWAYS
    # -------------------------------------------------

    with analysis_tabs[0]:

        st.subheader(
            f"Genes present in all {total_selected_cancers} "
            f"selected cancer pathways"
        )


        if total_selected_cancers <= 1:

            st.info(
                "Select at least two cancer pathways "
                "to perform a shared-gene comparison."
            )


        elif not genes_all_cancers:

            st.warning(
                "No genes were found to be common to all "
                "selected cancer pathways."
            )


        else:

            st.write(
                f"Found {len(genes_all_cancers)} genes "
                "common to all selected cancer pathways."
            )


            display_data = []


            for gene in genes_all_cancers:

                row = {

                    "Gene":
                        gene,

                    "Pathways":
                        ", ".join(
                            [
                                name
                                for name, genes
                                in p2g.items()
                                if gene in genes
                            ]
                        ),

                    **generate_icon_links(
                        gene
                    )
                }


                display_data.append(
                    row
                )


            df_display = pd.DataFrame(
                display_data
            ).sort_values(
                by="Gene"
            )


            if gene_filter:

                df_display = df_display[
                    df_display["Gene"]
                    .str.upper()
                    .str.contains(
                        gene_filter
                    )
                ]


            cols_display = [
                "Gene",
                "Pathways",
                "Gene Cards",
                "NCBI",
                "ENSEMBL",
                "GEO"
            ]


            st.write(
                df_display[
                    cols_display
                ].to_html(
                    escape=False,
                    index=False
                ),
                unsafe_allow_html=True
            )


    # -------------------------------------------------
    # 8B. MULTIPLE CANCER PATHWAYS
    # -------------------------------------------------

    with analysis_tabs[1]:

        st.subheader(
            "Genes present in multiple, but not all, "
            "selected cancer pathways"
        )


        if not genes_some_cancers:

            st.info(
                "No genes were found to be shared by "
                "multiple selected cancer pathways."
            )


        else:

            st.write(
                f"Found {len(genes_some_cancers)} genes "
                "shared by multiple selected cancer pathways."
            )


            display_data = []


            for gene in genes_some_cancers:

                row = {

                    "Gene":
                        gene,

                    "Cancer Pathways":
                        ", ".join(
                            sorted(
                                g2c[gene]
                            )
                        ),

                    **generate_icon_links(
                        gene
                    )
                }


                display_data.append(
                    row
                )


            df_display = pd.DataFrame(
                display_data
            ).sort_values(
                by="Gene"
            )


            if gene_filter:

                df_display = df_display[
                    df_display["Gene"]
                    .str.upper()
                    .str.contains(
                        gene_filter
                    )
                ]


            cols_display = [
                "Gene",
                "Cancer Pathways",
                "Gene Cards",
                "NCBI",
                "ENSEMBL",
                "GEO"
            ]


            st.write(
                df_display[
                    cols_display
                ].to_html(
                    escape=False,
                    index=False
                ),
                unsafe_allow_html=True
            )


    # -------------------------------------------------
    # 8C. ONE CANCER PATHWAY
    # -------------------------------------------------

    with analysis_tabs[2]:

        st.subheader(
            "Genes present in only one selected cancer pathway"
        )


        if not genes_one_cancer:

            st.info(
                "No genes were found to be unique to "
                "a single selected cancer pathway."
            )


        else:

            st.write(
                f"Found {len(genes_one_cancer)} genes "
                "present in only one selected cancer pathway."
            )


            display_data = []


            for gene in genes_one_cancer:

                row = {

                    "Gene":
                        gene,

                    "Cancer Pathway":
                        list(
                            g2c[gene]
                        )[0],

                    **generate_icon_links(
                        gene
                    )
                }


                display_data.append(
                    row
                )


            df_display = pd.DataFrame(
                display_data
            ).sort_values(
                by="Cancer Pathway"
            )


            if gene_filter:

                df_display = df_display[
                    df_display["Gene"]
                    .str.upper()
                    .str.contains(
                        gene_filter
                    )
                ]


            cols_display = [
                "Gene",
                "Cancer Pathway",
                "Gene Cards",
                "NCBI",
                "ENSEMBL",
                "GEO"
            ]


            st.write(
                df_display[
                    cols_display
                ].to_html(
                    escape=False,
                    index=False
                ),
                unsafe_allow_html=True
            )


    # -------------------------------------------------
    # 9. DOWNLOAD DATA
    # -------------------------------------------------

    st.markdown("---")

    st.subheader("Download Data")


    if genes_all_cancers:

        download_data = []


        for gene in genes_all_cancers:

            row = {

                "Gene Name":
                    gene,

                "Pathways":
                    ", ".join(
                        [
                            name
                            for name, genes
                            in p2g.items()
                            if gene in genes
                        ]
                    ),

                "Cancer Pathways":
                    ", ".join(
                        sorted(
                            g2c[gene]
                        )
                    ),

                "Frequency (All Selected Pathways)":
                    gene_counts[gene],

                **generate_url_links(
                    gene
                )
            }


            download_data.append(
                row
            )


        df_download = pd.DataFrame(
            download_data
        )


        column_order = [

            "Gene Name",

            "Pathways",

            "Cancer Pathways",

            "Frequency (All Selected Pathways)",

            "GeneCards_URL",

            "NCBI_URL",

            "ENSEMBL_URL",

            "GEO_URL"
        ]


        df_download = (
            df_download[
                column_order
            ]
            .sort_values(
                by="Gene Name"
            )
        )


        csv_string = (
            df_download
            .to_csv(
                index=False
            )
            .encode("utf-8")
        )


        st.download_button(

            label=(
                "Download Genes Present in All "
                "Selected Cancer Pathways"
            ),

            data=csv_string,

            file_name=(
                "shared_genes_all_selected_cancers.csv"
            ),

            mime="text/csv",

            key="dl_all_cancers"
        )


    else:

        st.write(
            "No genes present in all selected "
            "cancer pathways are available for download."
        )


# -----------------------------------------------------
# 10. NO PATHWAY SELECTED
# -----------------------------------------------------

else:

    st.info(
        "Please select at least one pathway from "
        "the sidebar to begin."
    )
