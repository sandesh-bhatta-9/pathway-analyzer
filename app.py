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

        if not isinstance(data, str) or not data.strip():
            raise ValueError("KEGG returned no pathway record")

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
        raise RuntimeError(f"Cannot retrieve or parse pathway {pathway_id}: {e}") from e


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


CHEMOTHERAPY_PATHWAYS = {
    "Platinum drug resistance": "hsa01524",
    "Antifolate resistance": "hsa01523",
    "EGFR tyrosine kinase inhibitor resistance": "hsa01521",
    "Pathways in cancer": "hsa05200",
    "Transcriptional misregulation in cancer": "hsa05202",
    "Chemical carcinogenesis - DNA adducts": "hsa05204",
    "Chemical carcinogenesis - reactive oxygen species": "hsa05208",
    "Chemical carcinogenesis - receptor activation": "hsa05207",
    "Base excision repair": "hsa03410",
    "Nucleotide excision repair": "hsa03420",
    "Mismatch repair": "hsa03430",
    "Homologous recombination": "hsa03440",
    "Non-homologous end-joining": "hsa03450",
    "p53 signaling pathway": "hsa04115",
    "Apoptosis": "hsa04210",
    "Cell cycle": "hsa04110",
    "Oocyte meiosis": "hsa04114",
    "DNA replication": "hsa03030",
    "Pyrimidine metabolism": "hsa00240",
    "Purine metabolism": "hsa00230",
    "One carbon pool by folate": "hsa00670",
    "Folate biosynthesis": "hsa00790",
    "Drug metabolism - cytochrome P450": "hsa00982",
    "Drug metabolism - other enzymes": "hsa00983",
    "Glutathione metabolism": "hsa00480",
    "Autophagy - animal": "hsa04140",
    "MAPK signaling pathway": "hsa04010",
    "NF-kappa B signaling pathway": "hsa04064",
    "HIF-1 signaling pathway": "hsa04066",
    "PI3K-Akt signaling pathway": "hsa04151",
    "TGF-beta signaling pathway": "hsa04350",
    "Rap1 signaling pathway": "hsa04015",
    "TNF signaling pathway": "hsa04668",
    "Ras signaling pathway": "hsa04014",
    "Central carbon metabolism in cancer": "hsa05230",
    "Choline metabolism in cancer": "hsa05231"
}


# Dataset supplied by Dr. Crasto. Associations are not individually verified.
# Embedded so app.py remains runnable without a separate CSV file.
NATURAL_PRODUCTS_CSV = 'Compound,Chemical_Class,Source_Organism,Kingdom_Domain,Pathway_Category,Primary_Molecular_Targets,KEGG_Pathway_ID,Mechanism_of_Action\nPhloretin,Dihydrochalcone (Flavonoid),Malus domestica (Apple),Plantae,Aerobic Glycolysis & Glucose Transport,"GLUT1, GLUT2 (via HNF6 downregulation)",hsa04931 / hsa00010,"Inhibits glucose transporters GLUT1/2, blocking glucose influx and triggering p53-mediated cell cycle arrest"\nGenistein,Isoflavone,Glycine max (Soybean),Plantae,Aerobic Glycolysis & Glucose Transport,"GLUT1, HK2, HIF-1alpha, Tyrosine kinases",hsa00010 / hsa04066,"Suppresses HIF-1alpha to downregulate GLUT1 and HK2, reducing glucose trapping and metastasis"\nCurcumin,Curcuminoid (Diarylheptanoid),Curcuma longa (Turmeric),Plantae,Aerobic Glycolysis & Lipogenesis,"HK2, ACLY, FASN, NF-kappaB, STAT3",hsa00010 / hsa00061,Detaches HK2 from outer mitochondrial membrane; inhibits ACLY and FASN to block de novo lipogenesis\nEpigallocatechin gallate (EGCG),Flavan-3-ol (Catechin),Camellia sinensis (Green tea),Plantae,Aerobic Glycolysis & Glucose Transport,"HK2, PFK1, DNMTs, VEGFR",hsa00010 / hsa05200,"Disrupts HK2-VDAC mitochondrial binding and inhibits PFK1, impairing glycolytic ATP generation"\nOleanolic Acid,Pentacyclic triterpenoid,Olea europaea (Olive leaves),Plantae,Aerobic Glycolysis & Glucose Transport,"HK2, PFK1, HIF-1alpha, AMPK",hsa00010 / hsa04066,Suppresses glycolytic flux and lactate production via HIF-1alpha downregulation and AMPK activation\nShikonin,Naphthoquinone,Lithospermum erythrorhizon,Plantae,Aerobic Glycolysis & Glucose Transport,PKM2 (Pyruvate kinase M2),hsa00010 / hsa05230,"Selectively inhibits PKM2 catalytic activity without affecting PKM1, starving tumor cells of ATP"\nLapachol,Naphthoquinone,Handroanthus impetiginosus (Lapacho tree),Plantae,Aerobic Glycolysis & Glucose Transport,"PKM2, Mitochondrial ETC",hsa00010 / hsa05230,Inhibits PKM2 and sensitizes cancer cells to mitochondrial uncoupling and ATP depletion\nGossypol,Polyphenolic binaphthyl dialdehyde,Gossypium spp. (Cottonseed),Plantae,Aerobic Glycolysis & Glucose Transport,"LDHA (Lactate dehydrogenase A), Bcl-2/Bcl-xL",hsa00010 / hsa00620,"Competitively inhibits NADH binding to LDHA, blocking lactate production and NAD+ regeneration"\nAstragaloside IV,Cycloartane-type triterpene glycoside,Astragalus membranaceus,Plantae,Aerobic Glycolysis & Glucose Transport,"c-Myc, GLUT1, HK2, LDHA",hsa05230,Downregulates c-Myc-driven transcription of core glycolytic transporters and enzymes\nWogonin,O-methylated flavone,Scutellaria baicalensis (Baikal skullcap),Plantae,Aerobic Glycolysis & Glucose Transport,"HK2, PDHK1, LDHA, p53",hsa00010 / hsa04115,Inhibits glycolysis via p53-dependent suppression of HK2 and PDHK1\nQuercetin,Flavonol,Allium cepa / Capparis spinosa,Plantae,Aerobic Glycolysis & Glucose Transport,"GLUT1, MCT1/4 (Monocarboxylate transporters), PKM2",hsa00010 / hsa05230,"Blocks glucose import and lactate export (MCT1/4), causing intracellular acidification"\nSilibinin (Silybin),Flavonolignan,Silybum marianum (Milk thistle),Plantae,Aerobic Glycolysis & Glucose Transport,"GLUT1, GLUT4, HIF-1alpha, STAT3",hsa04066 / hsa05230,Directly binds and blocks GLUT transporters and suppresses HIF-1alpha translation\nBerberine,Isoquinoline alkaloid,Coptis chinensis / Berberis vulgaris,Plantae,Lipid & Cholesterol Metabolism,"ACLY (ATP citrate lyase), AMPK, Complex I",hsa00020 / hsa04152,"Inhibits Mitochondrial Complex I and ACLY, depleting cytosolic acetyl-CoA and activating AMPK"\nLuteolin,Flavone,Apium graveolens / Petroselinum crispum,Plantae,Lipid & Cholesterol Metabolism,"FASN (Fatty acid synthase), ACC",hsa00061,"Inhibits beta-ketoacyl reductase domain of FASN, halting de novo palmitate synthesis"\n4-Methylumbelliferone,Coumarin,Apiaceae / Realgar-Coptis formulations,Plantae,Lipid & Cholesterol Metabolism,"FASN, HAS2, HAS3",hsa00061 / hsa00520,Dual inhibitor of de novo lipogenesis and UDP-glucuronic acid-dependent hyaluronan synthesis\nResveratrol,Stilbenoid,Vitis vinifera (Grapes),Plantae,Lipid & Cholesterol Metabolism,"SREBP-1c, FASN, SIRT1, AMPK",hsa04152 / hsa00061,Activates SIRT1/AMPK to suppress SREBP-1c proteolytic maturation and fatty acid synthesis\nNaringenin,Flavanone,Citrus spp.,Plantae,Lipid & Cholesterol Metabolism,"ACAT, HMG-CoA reductase, LXRalpha",hsa00100 / hsa04979,"Inhibits cholesterol esterification and promotes cholesterol efflux, disrupting oncogenic lipid rafts"\nCerulenin,Epoxy-carboxamide polyketide,Cephalosporium caerulens (Fungus),Fungi,Lipid & Cholesterol Metabolism,FASN (Beta-ketoacyl-ACP synthase domain),hsa00061,"Irreversibly alkylates the active-site cysteine of FASN, inducing malonyl-CoA accumulation and apoptosis"\nBetulin,Lupane-type triterpenoid,Betula spp. (Birch bark),Plantae,Lipid & Cholesterol Metabolism,SCAP / SREBP-1 & SREBP-2,hsa00061 / hsa00100,"Binds SCAP to retain SREBPs in the ER, blocking both fatty acid and cholesterol biosynthesis"\nLobetyolin,Polyacetylene glycoside,Codonopsis pilosula,Plantae,Amino Acid & Nucleotide Metabolism,ASCT2 (SLC1A5 glutamine transporter),hsa05230 / hsa00250,"Downregulates ASCT2, starving glutamine-addicted tumors of anaplerotic carbon and nitrogen"\nMorin,Flavonol,Maclura pomifera (Osage orange),Plantae,Amino Acid & Nucleotide Metabolism,"beta-catenin, c-Myc, GLS1",hsa05230 / hsa04310,Suppresses Wnt/beta-catenin and c-Myc signaling to co-inhibit glycolysis and glutaminolysis\nEsculetin,Dihydroxycoumarin,Fraxinus spp. (Ash tree),Plantae,Amino Acid & Nucleotide Metabolism,"beta-catenin, c-Myc, GLS1",hsa05230 / hsa04310,Disrupts c-Myc-mediated glutaminolysis and TCA cycle anaplerosis\n"3,3\'-Diindolylmethane (DIM)",Indole alkaloid derivative,Brassica spp. (Cruciferous vegetables),Plantae,Amino Acid & Nucleotide Metabolism,DHODH (Dihydroorotate dehydrogenase),hsa00240,"Inhibits mitochondrial DHODH, halting de novo pyrimidine synthesis and depleting UTP/CTP"\nRabdosianone I,ent-Kaurane diterpenoid,Isodon japonicus,Plantae,Amino Acid & Nucleotide Metabolism,"Thymidylate synthase (TS), ANT2, PHB2",hsa00240 / hsa00670,Promotes proteasomal degradation of thymidylate synthase and impairs mitochondrial ADP/ATP exchange\nHonokiol,Biphenolic lignan,Houpoea officinalis (Magnolia bark),Plantae,Amino Acid & Nucleotide Metabolism,"IDO1 / Tryptophan pathway, PPP, SIRT3",hsa00380 / hsa00030,"Activates mitochondrial SIRT3, suppresses pentose phosphate pathway flux, and blocks kynurenine immunosuppression"\n6-Diazo-5-oxo-L-norleucine (DON),Diazo amino acid analog,Streptomyces spp.,Bacteria,Amino Acid & Nucleotide Metabolism,"Glutaminase (GLS1/2), Amidotransferases",hsa00250 / hsa00230,Glutamine antagonist that irreversibly inhibits glutaminase and glutamine-utilizing purine/pyrimidine enzymes\nAzaserine,Diazoacetyl-L-serine,Streptomyces fragilis,Bacteria,Amino Acid & Nucleotide Metabolism,FGAR amidotransferase (Purine biosynthesis),hsa00230,Competitive glutamine antagonist blocking de novo purine nucleotide biosynthesis\nRhein,Anthraquinone,Rheum palmatum (Rhubarb),Plantae,"Mitochondrial OXPHOS, Redox & Ferroptosis","PDK1/2–PDH axis, Mitochondrial membrane potential",hsa00620 / hsa00190,"Inhibits PDK1/2 to reactivate PDH, forcing oxidative stress and cytochrome c-mediated apoptosis"\nBetulinic Acid,Lupane-type pentacyclic triterpenoid,Betula spp. (White birch bark),Plantae,"Mitochondrial OXPHOS, Redox & Ferroptosis","PDK–PDH axis, Mitochondrial PT pore",hsa00190 / hsa04210,Directly triggers mitochondrial membrane permeabilization and intrinsic apoptosis independent of p53\nGinsenoside Compound K,Dammarane-type triterpene saponin,Panax ginseng (intestinal metabolite),Plantae,"Mitochondrial OXPHOS, Redox & Ferroptosis","Mitochondrial ETC, SOD, Catalase",hsa00190 / hsa04210,Induces mitochondrial ROS burst exceeding endogenous antioxidant capacity\nArtemisinin / Dihydroartemisinin,Sesquiterpene lactone (Endoperoxide),Artemisia annua (Sweet wormwood),Plantae,"Mitochondrial OXPHOS, Redox & Ferroptosis","Intracellular Fe2+, GPX4, System Xc-, Akt/mTOR",hsa04216 / hsa04150,"Iron-dependent cleavage of endoperoxide bridge generates lipid peroxides, triggering ferroptosis"\nGliotoxin,Epipolythiodioxopiperazine (Mycotoxin),Trichoderma spp. / Marine fungi,Fungi,"Mitochondrial OXPHOS, Redox & Ferroptosis","HSP70–AR–MBOAT2 axis, GPX4",hsa04216,"Inhibits GPX4 and phospholipid remodeling, inducing lethal ferroptosis"\nPiperlongumine,Alkaloid amide,Piper longum (Long pepper),Plantae,"Mitochondrial OXPHOS, Redox & Ferroptosis","GSTP1 (Glutathione S-transferase Pi 1), TrxR1",hsa00480,"Depletes intracellular glutathione (GSH) and inhibits TrxR1, selectively elevating ROS in transformed cells"\nParthenolide,Germacranolide sesquiterpene lactone,Tanacetum parthenium (Feverfew),Plantae,"Mitochondrial OXPHOS, Redox & Ferroptosis","GSH depletion, IKKbeta / NF-kappaB, STAT3",hsa00480 / hsa04064,"Michael acceptor moiety covalently modifies thiol groups on GSH and IKKbeta, inducing oxidative apoptosis"\nWithaferin A,Withanolide (Steroidal lactone),Withania somnifera (Ashwagandha),Plantae,"Mitochondrial OXPHOS, Redox & Ferroptosis","KEAP1/Nrf2, GPX4, HSP90, Vimentin",hsa04216 / hsa04141,Inactivates GPX4 and targets KEAP1 to induce non-canonical ferroptosis and proteotoxic stress\nRotenone,Isoflavonoid (Rotenoid),Derris / Lonchocarpus spp.,Plantae,"Mitochondrial OXPHOS, Redox & Ferroptosis",Mitochondrial Complex I (NADH:ubiquinone oxidoreductase),hsa00190,"Blocks electron transfer from iron-sulfur centers in Complex I to ubiquinone, halting OXPHOS ATP"\nOligomycin A,Macrolide,Streptomyces diastatochromogenes,Bacteria,"Mitochondrial OXPHOS, Redox & Ferroptosis",ATP Synthase (Complex V Fo subunit),hsa00190,"Blocks the proton channel of ATP synthase, inhibiting mitochondrial oxidative phosphorylation"\nTanshinone I / IIA,Abietane diterpenoid quinone,Salvia miltiorrhiza (Danshen),Plantae,Master Metabolic Signaling Hubs,"PI3K, Akt, mTORC1, STAT3",hsa04151 / hsa04150,"Inhibits PI3K/Akt/mTORC1 phosphorylation, shutting down anabolic translation and survival signaling"\nHalofuginone,Quinazolinone alkaloid derivative,Dichroa febrifuga,Plantae,Master Metabolic Signaling Hubs,"EPRS (Prolyl-tRNA synthetase), GCN2, Akt/mTORC1",hsa00970 / hsa04150,Inhibits EPRS to trigger GCN2-mediated amino acid starvation response and mTORC1 suppression\nCapsaicin,Capsaicinoid (Vanilloid alkaloid),Capsicum annuum (Chili pepper),Plantae,Master Metabolic Signaling Hubs,"NF-kappaB, STAT3, TRPV1, Complex I/III",hsa04064 / hsa04630,Blocks NF-kappaB/STAT3 transcriptional activation and disrupts mitochondrial electron transport\nSulforaphane,Isothiocyanate,Brassica oleracea (Broccoli sprouts),Plantae,Master Metabolic Signaling Hubs,"HDACs, DNMTs, KEAP1/Nrf2, PKM2",hsa05202 / hsa04141,Inhibits HDACs/DNMTs to reactivate tumor-suppressor genes and modulates Nrf2 antioxidant response\nRapamycin (Sirolimus),Macrolide lactone,Streptomyces hygroscopicus,Bacteria,Master Metabolic Signaling Hubs,FKBP12–mTORC1 complex,hsa04150 / hsa04140,"Allosterically inhibits mTORC1, halting cap-dependent translation, lipogenesis, and inducing autophagy"\nWortmannin,Furanosteroid metabolite,Penicillium funiculosum,Fungi,Master Metabolic Signaling Hubs,PI3K (Phosphoinositide 3-kinase catalytic subunit p110),hsa04151,"Covalent, irreversible inhibitor of PI3K, blocking downstream Akt/mTOR metabolic signaling"\nTrichostatin A (TSA),Hydroxamic acid,Streptomyces hygroscopicus,Bacteria,Master Metabolic Signaling Hubs,Class I & II Histone Deacetylases (HDACs),hsa05202,"Chelates catalytic zinc ion in HDACs, reversing epigenetic silencing of cell-cycle and metabolic regulators"\nRomidepsin (FK228),Bicyclic depsipeptide,Chromobacterium violaceum,Bacteria,Master Metabolic Signaling Hubs,"Class I HDACs (HDAC1, HDAC2)",hsa05202,Prodrug reduced intracellularly to release a zinc-binding thiol that inhibits Class I HDACs\nPaclitaxel (Taxol),Diterpenoid taxane,Taxus brevifolia (Pacific yew),Plantae,Cytoskeletal & Mitotic Dynamics,Beta-tubulin,hsa04110 / hsa04540,"Binds beta-tubulin lumen to hyper-stabilize microtubules, preventing depolymerization and arresting mitosis"\nVinblastine & Vincristine,Bis-indole (Vinca) alkaloids,Catharanthus roseus (Madagascar periwinkle),Plantae,Cytoskeletal & Mitotic Dynamics,Alpha/Beta-tubulin heterodimers,hsa04110 / hsa04540,"Binds vinca domain on tubulin to inhibit microtubule polymerization, collapsing mitotic spindles"\nEribulin Mesylate,Macrocyclic ketone (Halichondrin B analog),Halichondria okadai (Marine sponge),Animalia (Porifera),Cytoskeletal & Mitotic Dynamics,Microtubule plus-ends (Beta-tubulin),hsa04110,"Inhibits microtubule growth phase without affecting shortening, sequestering tubulin into non-productive aggregates"\nColchicine,Tropolone alkaloid,Colchicum autumnale (Autumn crocus),Plantae,Cytoskeletal & Mitotic Dynamics,Beta-tubulin (Colchicine-binding site),hsa04110,"Binds soluble tubulin dimers to induce curved conformation, blocking microtubule assembly"\nEpothilone B (Ixabepilone),16-membered macrolide,Sorangium cellulosum (Myxobacterium),Bacteria,Cytoskeletal & Mitotic Dynamics,Beta-tubulin (Taxane pocket),hsa04110,Stabilizes microtubules similarly to paclitaxel while retaining activity in P-glycoprotein-overexpressing tumors\nCamptothecin,Pentacyclic quinoline alkaloid,Camptotheca acuminata (Happy tree),Plantae,DNA Topology & Replication Stress,DNA Topoisomerase I,hsa03030 / hsa03440,"Stabilizes covalent Topo I-DNA cleavage complexes, causing lethal double-strand breaks during S-phase"\nEtoposide & Teniposide,Podophyllotoxin lignan glycosides,Podophyllum peltatum (Mayapple),Plantae,DNA Topology & Replication Stress,DNA Topoisomerase II,hsa03030 / hsa03440,"Forms ternary complex with Topo II and DNA, preventing religation of cleaved DNA strands"\nDoxorubicin & Daunorubicin,Anthracycline polyketides,Streptomyces peucetius,Bacteria,DNA Topology & Replication Stress,"DNA intercalation, Topoisomerase II, Iron-mediated ROS",hsa03440 / hsa04210,"Intercalates between DNA base pairs, poisons Topo II, and generates cardiotoxic/cytotoxic semiquinone radicals"\nTrabectedin (Yondelis),Tetrahydroisoquinoline alkaloid,Ecteinascidia turbinata (Marine tunicate),Animalia (Chordata),DNA Topology & Replication Stress,"DNA minor groove (N2 of guanine), TC-NER machinery",hsa03420,"Alkylates guanine N2 in the DNA minor groove, trapping transcription-coupled nucleotide excision repair complexes"\nBleomycin,Glycopeptide antibiotic,Streptomyces verticillus,Bacteria,DNA Topology & Replication Stress,DNA (Fe2+/O2-dependent strand scission),hsa03450 / hsa04115,Chelates Fe(II) and reduces molecular oxygen to cleave DNA phosphodiester bonds\nMitomycin C,Mitomycin quinone,Streptomyces caespitosus,Bacteria,DNA Topology & Replication Stress,DNA (GpC inter-strand cross-linking),hsa03460,Bioreductively activated under hypoxia to form covalent inter-strand DNA cross-links\nActinomycin D (Dactinomycin),Chromopeptide lactone,Streptomyces parvulus,Bacteria,DNA Topology & Replication Stress,"DNA (GC-rich transcription initiation complexes), RNA Polymerase",hsa03020,"Intercalates into GC-rich DNA duplexes, blocking RNA polymerase elongation"\nSalinosporamide A (Marizomib),Beta-lactone-gamma-lactam,Salinispora tropica (Marine actinomycete),Bacteria,Proteostasis & Translation,"20S Proteasome (beta1, beta2, beta5 catalytic subunits)",hsa03050,Irreversibly acylates the N-terminal threonine of all three 20S proteasome proteolytic subunits\nHomoharringtonine (Omacetaxine),Cephalotaxine alkaloid ester,Cephalotaxus harringtonia (Plum yew),Plantae,Proteostasis & Translation,80S Ribosome (A-site cleft),hsa03010,"Inhibits initial elongation step of protein synthesis, rapidly depleting short-lived oncoproteins (Mcl-1, c-Myc)"\nPlitidepsin (Aplidine),Cyclic depsipeptide,Aplidium albicans (Marine tunicate),Animalia (Chordata),Proteostasis & Translation,eEF1A2 (Eukaryotic translation elongation factor 1A2),hsa03013 / hsa04010,"Targets eEF1A2 to disrupt protein folding/transport, activating Rac1/JNK oxidative apoptosis"\nGeldanamycin / Herbimycin A,Benzoquinone ansamycin,Streptomyces hygroscopicus,Bacteria,Proteostasis & Translation,HSP90 (N-terminal ATP-binding pocket),hsa04141,"Inhibits HSP90 chaperone ATPase activity, triggering ubiquitin-proteasomal degradation of client oncoproteins"\nCelastrol,Quinone methide triterpenoid,Tripterygium wilfordii (Thunder god vine),Plantae,Proteostasis & Translation,"HSP90-Cdc37 complex, 26S Proteasome, NF-kappaB",hsa04141 / hsa03050,Disrupts HSP90-Cdc37 protein-protein interaction and inhibits chymotrypsin-like proteasome activity\nBafilomycin A1,Plecomacrolide,Streptomyces griseus,Bacteria,Proteostasis & Translation,Vacuolar-type H+-ATPase (V-ATPase),hsa04140 / hsa04142,"Inhibits V-ATPase, preventing lysosomal acidification and autophagosome-lysosome fusion"\n'

@st.cache_data
def load_natural_products():
    from io import StringIO
    return pd.read_csv(StringIO(NATURAL_PRODUCTS_CSV)).fillna("")

def parse_human_gene_catalog(raw):
    """KEGG organism lists have ID, gene type, position, and description columns.
    Legacy two-column lists are accepted. Only parse the final description column.
    """
    import re
    lookup = {}
    primary = {}
    if not isinstance(raw, str):
        raise ValueError("Human gene catalog unavailable")
    for line in raw.splitlines():
        fields = line.split("\t")
        if len(fields) not in (2, 4) or not re.fullmatch(r"hsa:\d+", fields[0]):
            continue
        description = fields[-1]
        if ";" not in description:
            continue
        names = [v.strip() for v in description.split(";", 1)[0].split(",")]
        if not names or not re.fullmatch(r"[A-Za-z0-9_.-]+", names[0]):
            continue
        primary.setdefault(names[0].casefold(), set()).add(names[0])
        for name in names:
            if re.fullmatch(r"[A-Za-z0-9_.-]+", name):
                lookup.setdefault(name.casefold(), set()).add(names[0])
    lookup.update(primary)  # Primary symbols take precedence over historical aliases.
    if not lookup:
        raise ValueError("No valid human gene symbols parsed; comparison disabled")
    return lookup

@st.cache_data(ttl=86400)
def human_gene_lookup():
    return parse_human_gene_catalog(KEGG().list("hsa"))

# Explicit name normalization only; this does not verify biological effects.
TARGET_ALIASES = {
    "glut1": "SLC2A1", "glut2": "SLC2A2", "glut4": "SLC2A4",
    "hif-1alpha": "HIF1A", "c-myc": "MYC", "p53": "TP53",
    "beta-catenin": "CTNNB1", "asct2": "SLC1A5", "trxr1": "TXNRD1",
    "srebp-1c": "SREBF1", "lxralpha": "NR1H3", "gls1": "GLS",
    "thymidylate synthase": "TYMS", "dna topoisomerase i": "TOP1",
    "vimentin": "VIM", "eef1a2": "EEF1A2", "gcn2": "EIF2AK4",
}

def resolve_targets(text, lookup):
    import re
    rows = []
    # Commas within parentheses are annotations, not separate targets.
    tokens = re.split(r",(?![^()]*\))", text)
    for token in tokens:
        token = token.strip()
        name = re.sub(r"\s*\([^)]*\)", "", token).strip()
        candidate = TARGET_ALIASES.get(name.casefold(), name)
        # Family/complex names must not be collapsed to an arbitrary catalog alias.
        ambiguous = {"nf-kappab", "nf-kappab", "nf-kb", "nf-κb", "ampk", "akt", "pi3k",
                     "hsp90", "hsp70", "pfk1", "dnmts", "vegfr", "acc", "acat",
                     "sod", "hdacs", "beta-tubulin", "dna topoisomerase ii"}
        matches = set() if name.casefold() in ambiguous else lookup.get(candidate.casefold(), set())
        symbol = next(iter(matches)) if len(matches) == 1 else ""
        rows.append({"Reported target": token, "Gene": symbol,
                     "Mapping status": "Gene name matched (not evidence of compound activity)" if symbol
                     else "Unresolved: manual review required"})
    return rows

# Literature audit: primary-source abstracts/publisher record checked for these five claims.
# Support is experimental and context-specific, not clinical efficacy or direct binding.
CURCUMIN_EVIDENCE = {
    "HK2": {"Evidence status": "Supported: cancer cell study", "Cancer/model": "HCT116 and HT29 colorectal cancer cells",
        "Finding": "Reduced HKII expression/activity and mitochondrial dissociation reported",
        "Reference": "https://pubmed.ncbi.nlm.nih.gov/25229889/", "DOI": "10.1097/CAD.0000000000000132",
        "Review scope": "Primary research abstract checked; full-text/quality audit pending", "Eligible": True},
    "FASN": {"Evidence status": "Supported: cancer cell study", "Cancer/model": "MDA-MB-231 breast cancer cells",
        "Finding": "Reduced FAS activity, expression and mRNA; apoptosis reported",
        "Reference": "https://pubmed.ncbi.nlm.nih.gov/26985864/", "DOI": "10.3892/or.2016.4682",
        "Review scope": "Primary research abstract checked; full-text/quality audit pending", "Eligible": True},
    "STAT3": {"Evidence status": "Supported: cancer cell study", "Cancer/model": "Hodgkin lymphoma cell models",
        "Finding": "Suppression of STAT3 activation reported; direct binding not established here",
        "Reference": "https://pubmed.ncbi.nlm.nih.gov/18386790/", "DOI": "10.1002/ijc.23477",
        "Review scope": "Primary research abstract checked; full-text/quality audit pending", "Eligible": True},
    "NF-kappaB": {"Evidence status": "Supported signaling effect; gene identity unresolved", "Cancer/model": "Hodgkin lymphoma cell models",
        "Finding": "Suppression of NF-kappaB activation reported; family/complex not assigned to one gene",
        "Reference": "https://pubmed.ncbi.nlm.nih.gov/18386790/", "DOI": "10.1002/ijc.23477",
        "Review scope": "Primary research abstract checked; full-text/quality audit pending", "Eligible": False},
    "ACLY": {"Evidence status": "Context mismatch: cancer-specific claim not established by checked source",
        "Cancer/model": "NAFLD models; primary mouse hepatocytes, mice, HepG2/HEK293T and recombinant ACLY",
        "Finding": "Citrate metabolism/lipid effects studied for fatty liver disease; not validation of an anticancer ACLY mechanism",
        "Reference": "https://doi.org/10.1021/acs.jafc.1c03105", "DOI": "10.1021/acs.jafc.1c03105",
        "Review scope": "Primary publisher abstract checked; full-text/quality audit pending", "Eligible": False},
}

def annotate_evidence(compound, rows):
    for row in rows:
        evidence = CURCUMIN_EVIDENCE.get(row["Reported target"]) if compound == "Curcumin" else None
        row.update(evidence or {"Evidence status": "Pending: association not yet reviewed", "Cancer/model": "",
            "Finding": "Supplied dataset claim only", "Reference": "", "DOI": "",
            "Review scope": "Not reviewed", "Eligible": False})
    return rows

natural_df = load_natural_products()
natural_records = natural_df.set_index("Compound").to_dict("index")
# Values are listed pathway associations, never compound target gene sets.
NATURAL_PRODUCT_PATHWAYS = {
    name: row["KEGG_Pathway_ID"] for name, row in natural_records.items()
}



cancer_options = {
    name: pid
    for name, pid in all_paths.items()
    if any(
        keyword in name.lower()
        for keyword in CANCER_KEYWORDS
    )
}


chemo_options = CHEMOTHERAPY_PATHWAYS


natural_options = NATURAL_PRODUCT_PATHWAYS


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
            if k == "Curcumin"
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
    format_func=lambda name: (
        f"{name} ({cancer_options[name]})"
    ),
    key="cancer_key"
)


st.sidebar.subheader("Therapeutic Pathways")

sel_chemo = st.sidebar.multiselect(
    "Select pharmaceutical pathways:",
    list(chemo_options.keys()),
    format_func=lambda name: (
        f"{name} ({chemo_options[name]})"
    ),
    key="chemo_key"
)


sel_natural = st.sidebar.multiselect(
    "Select natural compounds:",
    list(natural_options.keys()),
    format_func=lambda name: (
        f"{name} ({natural_options[name]})"
    ),
    key="natural_key"
)


combined_options = {
    **cancer_options,
    **chemo_options,
    **natural_options
}


selected_pathways = (
    list(dict.fromkeys(sel_cancer + sel_chemo + sel_natural))
)


# -----------------------------------------------------
# 6. MAIN APP LOGIC
# -----------------------------------------------------

target_tables = {}
include_pending = st.sidebar.checkbox("Include mapped associations awaiting cancer-evidence review", value=False,
    help="Default comparisons use only currently reviewed cancer-cell associations. This is not clinical validation.")
if sel_natural:
    st.header("Natural compound details")
    st.caption("Reported associations from the supplied dataset, with a limited Curcumin evidence audit. "
               "Direct binding and clinical effectiveness are not established by these comparisons. Pathway membership is not a compound target claim.")
    try:
        lookup = human_gene_lookup()
    except Exception as exc:
        lookup = {}
        st.error(f"Cannot load the human gene catalog: {exc}. Target comparisons are unavailable; retry later.")
    for compound in sel_natural:
        record = natural_records[compound]
        target_tables[compound] = annotate_evidence(compound, resolve_targets(record["Primary_Molecular_Targets"], lookup))
        with st.expander(compound, expanded=True):
            st.write("Source:", record["Source_Organism"])
            st.write("Chemical class:", record["Chemical_Class"])
            st.write("Original supplied mechanism (not fully verified):", record["Mechanism_of_Action"])
            if compound == "Curcumin":
                st.info("Evidence audit: HK2, FASN and STAT3 effects are supported by cancer-cell studies. "
                        "NF-kappaB is retained as a signaling effect, not one gene. "
                        "The checked ACLY study concerns fatty liver disease; the anticancer ACLY claim remains unestablished. "
                        "These studies do not establish oral treatment effectiveness or melanoma-specific activity.")
            st.write("Reported pathway category:", record["Pathway_Category"])
            for pid in __import__("re").findall(r"hsa\d{5}", record["KEGG_Pathway_ID"]):
                st.markdown(f"Associated pathway: [{pid}](https://www.kegg.jp/pathway/{pid})")
            display_targets = pd.DataFrame(target_tables[compound]).drop(columns=["Eligible"])
            st.dataframe(display_targets, use_container_width=True)
            for evidence in target_tables[compound]:
                if evidence["Reference"]:
                    st.markdown(f"{evidence['Reported target']}: [supporting source]({evidence['Reference']}) — {evidence['Evidence status']}")
            st.caption("Evidence audit currently covers Curcumin only. Other entries remain pending. Gene-name matching and biological evidence are separate checks.")
    details = natural_df[natural_df["Compound"].isin(sel_natural)].copy()
    details["Verification_Status"] = details["Compound"].map(lambda c: "Partial primary-abstract audit; see target evidence export" if c == "Curcumin" else "Unreviewed supplied associations")
    st.download_button("Download selected compound details", details.to_csv(index=False).encode(),
                       "selected_compound_details.csv", "text/csv")
    with st.expander("Starting bibliography (not row-level verification)"):
        for label, doi in [
            ("Mangal et al., NPACT (2013 issue; online 2012)", "10.1093/nar/gks1047"),
            ("Chunarkar-Patil et al. (2024)", "10.3390/biomedicines12010201"),
            ("Koehn and Carter (2005)", "10.1038/nrd1657"),
            ("Newman and Cragg (2012)", "10.1021/np200906s"),
            ("Papadopoulou et al. (2026), corrected DOI", "10.3390/ph19060910"),
        ]:
            st.markdown(f"[{label}](https://doi.org/{doi})")

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

                for pathway in dict.fromkeys(sel_cancer + sel_chemo)
            }

            for future in concurrent.futures.as_completed(
                future_to_pathway
            ):

                pathway = future_to_pathway[future]

                try:
                    p2g[pathway] = future.result()
                except Exception as exc:
                    st.error(f"Pathway data unavailable: {exc}. Comparison stopped to avoid misleading zero overlaps.")
                    st.stop()


    for compound in sel_natural:
        p2g[compound] = {row["Gene"] for row in target_tables[compound]
                         if row["Gene"] and (include_pending or row["Eligible"])}
        if not p2g[compound]:
            st.warning(f"{compound}: no eligible resolved targets under the current evidence filter. "
                       "An empty comparison does not establish absence of biological overlap.")

    if sel_natural:
        st.header("Compound targets shared with cancer pathways")
        comparison = []
        for compound in sel_natural:
            for cancer in sel_cancer:
                shared = sorted(p2g[compound] & p2g.get(cancer, set()))
                comparison.append({"Compound": compound, "Cancer pathway": cancer,
                    "Matched reported targets": len(p2g[compound]), "Shared target count": len(shared) if p2g[compound] else None,
                    "Shared genes": ", ".join(shared),
                    "Status": ("No eligible targets; overlap unavailable" if not p2g[compound]
                               else "Descriptive overlap; pending associations included" if include_pending
                               else "Descriptive overlap of reviewed cell-study associations; cancer contexts may differ")})
        if comparison:
            comparison_df = pd.DataFrame(comparison)
            st.dataframe(comparison_df, use_container_width=True)
            st.download_button("Download compound–cancer overlaps", comparison_df.to_csv(index=False).encode(),
                               "compound_cancer_overlaps.csv", "text/csv")
        else:
            st.info("Select a cancer pathway to compare with the reported compound targets.")
        mapping_df = pd.DataFrame([{"Compound": c, **row} for c in sel_natural for row in target_tables[c]])
        st.download_button("Download target mapping and unresolved names", mapping_df.to_csv(index=False).encode(),
                           "compound_target_mapping.csv", "text/csv")

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

    st.header("Compound Target and Pathway-Gene Network")

    st.info(
        "The network combines KEGG pathway membership with unverified reported compound target associations. "
        "Green squares are compounds, not pathways. "
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
                "<br>Listed pathway association: "
                + combined_options[pathway]
                + "<br>Node type: "
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


        if cancers_for_gene and any(gene in p2g.get(c, set()) for c in sel_natural):
            color = "gold"
            category += "; also a matched reported compound target (association unverified)"
        gene_color.append(color)


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
            "Integrated Cancer, Natural Compound, "
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

        **Natural compounds**  
        Green square nodes represent compounds linked only to matched reported targets.
        Listed associated KEGG pathways are context, not expanded compound target sets.

        **Pharmaceutical pathways**  
        Blue square nodes represent pharmaceutical or chemotherapy-related pathways.


        ### Gene Nodes

        **Gold**  
        Matched reported compound target present in a selected cancer pathway.

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

        **Compound and pharmaceutical pathway connections**  
        Compound edges indicate unverified reported target associations. Pharmaceutical edges indicate pathway membership. These connections provide a
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
                        gene_filter, regex=False
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
                        gene_filter, regex=False
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
                        gene_filter, regex=False
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
