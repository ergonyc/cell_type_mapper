
# %%
import anndata as ad
import h5py
import json
import numpy as np
import pandas as pd
from pathlib import Path
import os
import ijson
import scanpy as sc
import anndata as ad
import scipy.sparse as sp


# %% [markdown]
# %%
from cell_type_mapper.cli.from_specified_markers import FromSpecifiedMarkersRunner

from cell_type_mapper.cli.map_to_on_the_fly_markers import OnTheFlyMapper

from cell_type_mapper.cli.precompute_stats_scrattch import PrecomputationScrattchRunner
from cell_type_mapper.cli.reference_markers import ReferenceMarkerRunner
from cell_type_mapper.cli.query_markers import QueryMarkerRunner



# %%


HOME = Path.home()

ROOT = HOME / "Projects/ASAP/cell_type_mapper"




CHUNK_SIZE = 40000
N_RUNNERS_UP = 5
RNG_SEED = 11235813
N_PROCESSORS = 6
MAX_GB = 32.0

#!/usr/bin/env python
# coding: utf-8

### set up train and test AnnData objects for LBL8R

# In[ ]:
### import local python functions in ../lbl8r


# In[ ]:
## e.g. for xylena data
XYLENA2_FULL = "xyl2_full.h5ad"

XYLENA2_TRAIN = "xyl2_train.h5ad"
XYLENA2_TEST = "xyl2_test.h5ad"
XYLENA2_QUERY_A = "xyl2_query.h5ad"
XYLENA2_QUERY_B = "xyl2_queryB.h5ad"


# linux
XYLENA2_PATH = "scdata/xylena"




## Make our ATLAS_DIR
ATLAS_NAME = "XYLENA2"
ATLAS_DIR = ROOT / ATLAS_NAME


# In[ ]:

# mac

# linux
data_path = ROOT / XYLENA2_PATH
# In[ ]: Load raw data
########################
# 0. LOAD TRAIN DATA
########################

train_filen = data_path / XYLENA2_TRAIN
# train_ad = ad.read_h5ad(train_filen, backed="r")
train_ad = ad.read_h5ad(train_filen)
# # In[ ]:  make test anndata

# test_filen = data_path / XYLENA2_TEST
# test_ad = ad.read_h5ad(test_filen, backed="r")
# queryA_filen = data_path / XYLENA2_QUERY_A
# query_ad = ad.read_h5ad(queryA_filen, backed="r")
# queryB_filen = data_path / XYLENA2_QUERY_B
# queryB_ad = ad.read_h5ad(queryB_filen, backed="r")
# %%  MAKE a local copy of the training data which is encoded with ensemble_ids


train_ad.var["gene_name"] = train_ad.var.index
train_ad.var_names = train_ad.var["gene_ids"]

# also coudl use lbl8r.model.utils._data.sparsify_adata

train_ad.X = sp.csr_matrix(train_ad.X)

# TODO: make recoded copies of these...
# test_filen = data_path / XYLENA2_TEST
# test_ad = ad.read_h5ad(test_filen, backed="r")
# queryA_filen = data_path / XYLENA2_QUERY_A
# query_ad = ad.read_h5ad(queryA_filen, backed="r")
# queryB_filen = data_path / XYLENA2_QUERY_B
# queryB_ad = ad.read_h5ad(queryB_filen, backed="r")



# In[ ]:  Make taxonomy

def class_assign(x):
    if x in ['ExN', 'InN']:
        return 'neuron'
    elif x in ['Oligo', 'VC', 'OPC', 'Astro', 'MG']:
        return 'non-neuronal'



# mapping = {
#     "oligo": "Oligo",
#     "opc": "OPC",
#     "glutamatergic": "ExN",
#     "gabergic": "InN",
#     "astrocyte": "Astro",
#     "immune": "MG",
#     "blood": "VC",
# }

subc_mapping = {
    "Oligo": "oligo",
    "OPC": "opc",
    "ExN": "glutamatergic",
    "InN": "gabaergic",
    "Astro": "astrocyte",
    "MG": "immune",
    "VC": "blood_vessel",
}


train_ad.obs['class'] = train_ad.obs['cell_type'].apply(class_assign)
train_ad.obs['subclass'] = train_ad.obs['cell_type'].map(subc_mapping)


heirarchy = ["class", "subclass"]

# %%

if not ATLAS_DIR.exists():
    ATLAS_DIR.mkdir()

atlas_filenm = ATLAS_DIR / f"{ATLAS_NAME}_atlas_adata.h5ad"
train_ad.write_h5ad(atlas_filenm)
# train_ad.write_h5ad(ATLAS_DIR / f"{ATLAS_NAME}_atlas_adata.h5ad")


# %% make scrattch reference


stats_path = ATLAS_DIR / f"{ATLAS_NAME}_precomputed_stats.h5"

precomputation_config = {
    'hierarchy': heirarchy,
    'h5ad_path': f"{atlas_filenm}",
    'output_path': f"{stats_path}",
    'n_processors': N_PROCESSORS,
    'normalization': 'raw',
}

# # %%
# !python -m cell_type_mapper.cli.precompute_stats_scrattch \
# --hierarchy f"{heirarchy}" \
# --h5ad_path f"{atlas_filenm}"\
# --n_processors f"{N_PROCESSORS}" \
# --normalization raw \
# --clobber True \
# --output_path f"{stats_path}" 
# %%

#


# the args=[] is important to prevent the PrecomputationScrattchRunner from grabbing
# any command line arguments when you invoke your python script
precomputation_runner = PrecomputationScrattchRunner(
    args=[], input_data=precomputation_config)
# %%
precomputation_runner.run()

print('\n====done with precomputed_stats====\n')
# %%

output_dir = ATLAS_DIR / "reference"

reference_config = {
    'precomputed_path_list': [stats_path],
    'n_valid': 20,
    'output_dir': output_dir,
    'clobber': True
}

reference_runner = ReferenceMarkerRunner(
    args=[], input_data=reference_config)
# %%
reference_runner.run()
print('\n====done with reference_markers====\n')

# %%

query_markers_path = ATLAS_DIR / f"references/query_markers.json"
reference_markers_path = ATLAS_DIR / f"references/reference_markers.h5"

query_config = {
    'output_path': f"{query_markers_path}",
    'reference_marker_path_list': [f"{reference_markers_path}"],
    'n_per_utility': 10,
    'n_processors': N_PROCESSORS,
}

query_runner = QueryMarkerRunner(
    args=[], input_data=query_config)

query_runner.run()


# # In[ ]:  Make taxonomy
# # %%
# # df = make_markers_table_v2(dfs, cell_types, all_genes)
# ### pack the Frontal Coretex defined taxa into a cell_assign table

# ######### from SEA-AD
# GABA_subtypes = ["Lamp5", "Sngc", "Sst", "Sst Chodl", "Pvalb", "Meis2"]
# glut_subtypes = [
#     "L2/3 IT",
#     "L4/5 IT",
#     "L5 IT",
#     "L6 IT",
#     "L6 IT Car3",
#     "L6 ET",
#     "L6 CT",
#     "L6b",
#     "L5/6 NP",
# ]
# non_neural_subtypes = [
#     "Endo",
#     "VLMC",
#     "SMC",
#     "Peri",
#     "Micro",
#     "PVM",
# ]  # Endothelial cells, vascular leptomeningeal cells, smooth muscle cells, pericytes, microglia, perivascular myeloid cells




# # %%    MAPPING CODE BELOW...


# TMP_DIR = ROOT / "tmp"
# if not TMP_DIR.exists():
#     TMP_DIR.mkdir()



# # %%  force CPU
# os.environ["AIBS_BKP_USE_TORCH"] = "false"

# os.environ["NUMEXPR_NUM_THREADS"] = "1"
# os.environ["MKL_NUM_THREADS"] = "1"
# os.environ["OMP_NUM_THREADS"] = "1"

# # %%
# FILE_ROOT = XYLENA2_TRAIN.rstrip(".h5ad")

# DATE = "20250207"


# ## SEA-AD taxonomy reference
# #%%  SEA-AD reference
# REFERENCE = "SEAAD"

# EXTENDED_RESULTS = f"{FILE_ROOT}.mmc.{REFERENCE}.{DATE}.json"
# LOG_FILE = f"{FILE_ROOT}.mmc.{REFERENCE}_log.{DATE}.txt"
# CSV_RESULTS = f"{FILE_ROOT}.mmc.{REFERENCE}_results.{DATE}.csv"


# RESULTS_DIR = ROOT / f"MMC.{REFERENCE}_RESULTS"
# if not RESULTS_DIR.exists():
#     RESULTS_DIR.mkdir()

# PRECOMPUTED_STATS = (
#     f"examples/data/abc_atlas_data/precomputed_stats.20231120.sea_ad.MTG.h5"
# )



# config = {
#     # "query_path": f"{ASAP_DATA}/{FILE_ROOT}.h5ad",
#     "query_path": f"{TMP_DIR}/_{FILE_ROOT}.h5ad",
#     "tmp_dir": f"{TMP_DIR}",
#     "extended_result_path": f"{RESULTS_DIR / EXTENDED_RESULTS}",
#     "csv_result_path": f"{RESULTS_DIR / CSV_RESULTS}",
#     "log_path": f"{RESULTS_DIR / LOG_FILE}",
#     "cloud_safe": False,
#     "verbose_csv": True,
#     "n_processors": N_PROCESSORS,
#     "max_gb": MAX_GB,
#     "map_to_ensembl": True,
#     "type_assignment": {
#         "normalization": "raw",
#         "bootstrap_iteration": 100,
#         "bootstrap_factor": 0.5,
#         "chunk_size": CHUNK_SIZE,
#         "n_runners_up": N_RUNNERS_UP,
#         "rng_seed": RNG_SEED,
#     },
#     "precomputed_stats": {"path": f"{ROOT / PRECOMPUTED_STATS}"},
#     "reference_markers": {"log2_fold_min_th": 0.5},
#     "query_markers": {"n_per_utility": 15, "genes_at_a_time": 1},
# }


# #%%
# runner = OnTheFlyMapper(args=[], input_data=config)
# runner.run()

# # %% WHB reference

# REFERENCE = "WMB"

# #REFERENCE = "SEAAD"

# EXTENDED_RESULTS = f"{FILE_ROOT}.mmc.{REFERENCE}.{DATE}.json"
# LOG_FILE = f"{FILE_ROOT}.mmc.{REFERENCE}_log.{DATE}.txt"
# CSV_RESULTS = f"{FILE_ROOT}.mmc.{REFERENCE}_results.{DATE}.csv"


# RESULTS_DIR = ROOT / f"MMC.{REFERENCE}_RESULTS"
# if not RESULTS_DIR.exists():
#     RESULTS_DIR.mkdir()

# PRECOMPUTED_STATS = f"examples/data/abc_atlas_data/precomputed_stats.siletti.training.h5"

# QUERY_MARKERS = f"examples/data/abc_atlas_data/query_markers.n10.20240221800.json"

# config = {
#     "query_path": f"{TMP_DIR}/_{FILE_ROOT}.h5ad",
#     "tmp_dir": f"{TMP_DIR}",
#     "extended_result_path": f"{RESULTS_DIR / EXTENDED_RESULTS}",
#     "csv_result_path": f"{RESULTS_DIR / CSV_RESULTS}",
#     "log_path": f"{RESULTS_DIR / LOG_FILE}",
#     "map_to_ensembl": True,
#     "cloud_safe": False,
#     "verbose_csv": True,
#     # "n_processors": N_PROCESSORS,
#     "max_gb": MAX_GB,
#     "type_assignment": {
#         "normalization": "raw",
#         "bootstrap_iteration": 100,
#         "bootstrap_factor": 0.5,
#         "chunk_size": CHUNK_SIZE,
#         "n_processors": N_PROCESSORS,
#         # "n_runners_up": N_RUNNERS_UP,
#         # "rng_seed": RNG_SEED,
#     },
#     # "reference_markers": {"log2_fold_min_th": 0.5},
#     "precomputed_stats": {"path": f"{ROOT / PRECOMPUTED_STATS}"},
#     "query_markers": {"serialized_lookup": f"{ROOT / QUERY_MARKERS}"},
# }

# runner = FromSpecifiedMarkersRunner(
#     args=[],
#     input_data=config
# )

# runner.run()


# #%%%

# HOME = Path.home()

# CWD = Path.cwd()

# tmp_dir = CWD.parent / "tmp"
# if not tmp_dir.exists():
#     tmp_dir.mkdir()


# print(f"temporary directory is {tmp_dir}")

# # %% [markdown]
# # Refernce taxonomy is best defined in Extended Table 2 / "glossary" 
# # from https://doi.org/10.1038/s41586-021-03950-0 https://www.nature.com/articles/s41586-021-03950-0
# # 
# # 
# # In the cell below, we define a taxonomy with three classes, which are broken up into a total of 
# # five subclasses, which are broken up into a total of eight clusters. The inheritance structure of the 
# # taxonomy is encoded in the names of the taxon, so that `cluster_ABC` is a descendent of `subclass_AB` 
# # which is a descendant of `class_A` (so, `cluster_101`is a descendant of subclass `subclass_10` 
# # which is a descendant of `class_1`). We will create a characteristic gene profile across a panel of 
# # 90 genes for each of our classes, subclasses, and clusters. These characteristic profiles are not 
# # meant to be realistic. They are engineered so that there genes that are good at discriminating 
# # between all pairs of taxons at all levels in our taxonomy.
# # 
# # **You should just run this cell without bothering to really understand what it is doing. It only 
# # exists to create the reference data for our example.**

# # %%
# ASAP_DATA = (
#     HOME
#     / "Projects/ASAP/data/asap-curated-cohort-pmdbs-sc-rnaseq/pmdbs_sc_rnaseq/cohort_analysis"
# )
# ASAP_DATA = HOME / "Projects/ASAP/harmonized-wf-dev/data/cohort_analysis"

# fname = ASAP_DATA / "asap-cohort.merged_adata_object.h5ad"
# # fname = "/home/ergonyc/Projects/ASAP/data/asap-curated-cohort-pmdbs-sc-rnaseq/pmdbs_sc_rnaseq/cohort_analysis/asap-cohort.merged_adata_object.h5ad"
# adata = ad.read_h5ad(fname)
# adata.var["gene_name"] = adata.var.index
# adata.var_names = adata.var["gene_id"]

# # %%
# adata.write_h5ad(tmp_dir / "_asap-cohort.merged_adata_object.h5ad")

# # %% [markdown]
# # ## TODO:
# # Load the asap-cohort.final_adata_object.h5ad 
# # 
# # Merge in cell_type, umaps, etc into "full" adata....
# # 
# # Also make a "full_gene" adata subset to redact the filtered cells.  
# # - compare full vs. filtered typing for SEA-AD, WHOLE BRAIN types (in case "marker genes" change per type)
# #     - compare marker genes
# # - for filtered cells. compare full genes vs. filtered for HVG only  
# # 
# # 
# # 
# # ### CARD taxonomy
# # create CARD taxonomy and reference (mean expression per type + "marker genes")

# # %%
# adata.var_names

# # %%
# adata.var["gene_id"]

# # %%
# tax_path = CWD / "data/abc_atlas_data"

# cluster_annotation_path = tax_path / "cluster_annotation_term.csv"
# cluster_membership_path = tax_path / "cluster_to_cluster_annotation_membership.csv"

# cluster_path = tax_path / "cluster.csv"
# cluster_term_set_path = tax_path / "cluster_annotation_term_set.csv"

# cluster = pd.read_csv(cluster_path, index_col=0)
# cluster_term_set = pd.read_csv(cluster_term_set_path, index_col=0)

# cluster_annotation = pd.read_csv(cluster_annotation_path, index_col=0)
# cluster_membership = pd.read_csv(cluster_membership_path, index_col=0)

# # %%
# cluster_membership

# # %%
# cluster_term_set

# # %%
# cluster_annotation

# # %%
# cluster

# # %%
# cluster_annotation = data_dir / "WMB-taxonomy/20230630/cluster_annotation_term.csv"
# assert cluster_annotation.is_file()

# cluster_membership = cluster_annotation.parent / "cluster_to_cluster_annotation_membership.csv"
# assert cluster_membership.is_file()

# cell_metadata = data_dir / "WMB-10X/20230630/cell_metadata.csv"
# assert cell_metadata.is_file()

# # %%
# from cell_type_mapper.taxonomy.taxonomy_tree import TaxonomyTree

# taxonomy_tree = TaxonomyTree.from_data_release(
#     cell_metadata_path=self.args["cell_metadata_path"],
#     cluster_annotation_path=self.args["cluster_annotation_path"],
#     cluster_membership_path=self.args["cluster_membership_path"],
#     hierarchy=self.args["hierarchy"],
#     do_pruning=self.args["do_pruning"],
# )

# # %%
# ROOT = HOME / "Projects/ASAP/cell_type_mapper"
# RESULTS_DIR = ROOT / "MMC.SEAAD_RESULTS"
# FILE_ROOT = "asap-cohort.merged_adata_object"
# DATE = "20250129"

# EXTENDED_RESULTS = f"{FILE_ROOT}.mmc.seaad_results.{DATE}.json"

# # json_path = "/Users/ergonyc/Projects/ASAP/cell_type_mapper/examples/data/asap-cohort.merged_adata_object.mmc.seaad_results.20250123b.json"
# json_path = RESULTS_DIR / EXTENDED_RESULTS


# # with open(json_path, 'r') as f:
# #     parser = ijson.parse(f)
# #     for prefix, event, value in parser:
# #         if prefix == 'item.name':
# #             print(value)


# def read_key_from_json(filename, key):
#     with open(filename, "r") as f:
#         object = ijson.items(f, key)
#         return next(object, {})


# tax_tree = read_key_from_json(json_path, "taxonomy_tree")

# # %%
# tax_tree.keys()

# # %%
# # df = make_markers_table_v2(dfs, cell_types, all_genes)
# ### pack the Frontal Coretex defined taxa into a cell_assign table
# taxa = dict()


# heirarchy = ["type", "subtype"]

# # define top level types.
# top_class = ["neuronal", "non-neuneuronalral"]
# neuronal_subclass = ["glutamatergic", "gabaergic"]

# # FULL TAXONOMY (SEA-AD)
# non_neuronal_subclass = [
#     "oligo",
#     "astrocyte",
#     "OPC",
# ]  # oligodendrocytes, astrocytes, oligodendrocyte precursor cells


# GABA_subtypes = ["Lamp5", "Sngc", "Sst", "Sst Chodl", "Pvalb", "Meis2"]
# glut_subtypes = [
#     "L2/3 IT",
#     "L4/5 IT",
#     "L5 IT",
#     "L6 IT",
#     "L6 IT Car3",
#     "L6 ET",
#     "L6 CT",
#     "L6b",
#     "L5/6 NP",
# ]
# non_neural_subtypes = [
#     "Endo",
#     "VLMC",
#     "SMC",
#     "Peri",
#     "Micro",
#     "PVM",
# ]  # Endothelial cells, vascular leptomeningeal cells, smooth muscle cells, pericytes, microglia, perivascular myeloid cells

# # CARD Taxa conversion...
# blood_vessel = ["Endo", "Peri"]

# taxon_list = [
#     "type_0",
#     "type_1",
#     "type_2",
#     "type_3",
#     "type_4",
#     "type_5",
#     "cluster_00",
#     "cluster_01",
#     "cluster_10",
#     "cluster_20",
#     "cluster_30",
#     "cluster_40",
#     "cluster_50",
# ]

# neuron = dict(alias="neurons", neurons=neuron_subtypes)

# taxonomy = dict(
#     type_0=neuron,
#     type_1="astrocyte",
#     type_2="oligo",
#     type_3="opc",
#     type_4="immune",
#     type_5="blood_vessel",
# )

# astrocyte = ASTROCYTE
# protoplasmic_astrocyte = ASTROCYTE + astro_sub["protoplasmic"]
# fibrous_astrocyte = ASTROCYTE + astro_sub["fibrous"]

# immune = IMMUNE
# microglia = IMMUNE + immune_sub["microglia"]
# t_cell = IMMUNE + immune_sub["t_cell"]
# b_cell = IMMUNE + immune_sub["b_cell"]

# blood = BLOOD_VESSEL
# pericyte = BLOOD_VESSEL + blood_sub["pericytes"]
# endothelial = BLOOD_VESSEL + blood_sub["endothelial"]

# oligo = OLIGO
# opc = OPC
# unknown = []
# ###

# full_cell_types = [
#     "oligo",
#     "opc",
#     # neuron
#     "glutamatergic",
#     "gabergic",
#     # astrocyte
#     "protoplasmic_astrocyte",
#     "fibrous_astrocyte",
#     # immune
#     "microglia",
#     "t_cell",
#     "b_cell",
#     # blood_vessel
#     "pericyte",
#     "endothelial",
#     # unknown
#     "unknown",
# ]

# # In[ ]:
# # remap columns as above
# mapping = {
#     "oligo": "Oligo",
#     "opc": "OPC",
#     "glutamatergic": "ExN",
#     "gabergic": "InN",
#     "astrocyte": "Astro",
#     "immune": "MG",
#     "blood": "VC",
# }

# # %% [markdown]
# # 
# # from https://alleninstitute.github.io/abc_atlas_access/notebooks/WHB_cluster_annotation_tutorial.html
# # 
# # subcluster:
# # 
# # > The finest level of cell type definition in the human whole brain taxonomy. It is defined by applying Louvain clustering independently for each cluster. Cells within a subcluster share similar characteristics and belong to the same cluster. In some cases subclusters are distinct, and in other cases represent cell states or gradients within the same cluster.
# # 
# # cluster:
# # 
# # > An intermediate level of cell type definitions in the human whole brain taxonomy. It is defined by applying Louvain clustering independently for each supercluster. All cells within a subcluster belong to the same cluster.   Clusters group together similar subclusters and are highly distinct from one another.
# # 
# # supercluster:
# # 
# # > The top level of cell type definition in the human whole brain taxonomy. It is defined using Louvain clustering (Cytograph) with most superclusters determined by broad brain region and neurotransmitter type.  All cells within a cluster belong to the same supercluster. Supercluster provides a broader characterization of cell types.
# # 
# # neurotransmitter:
# # 
# # > Neurotransmitter terms are assigned to clusters based on the expression (or lack of expression) of one or more neurotransmitter-associated transporter or enzyme genes

# # %%
# rng = np.random.default_rng(22112)

# n_genes = 90

# characteristic_profiles = dict()

# taxon_list = [
#     "class_0",
#     "class_1",
#     "class_2",
#     "subclass_00",
#     "subclass_01",
#     "subclass_10",
#     "subclass_11",
#     "subclass_20",
#     "cluster_000",
#     "cluster_001",
#     "cluster_010",
#     "cluster_011",
#     "cluster_100",
#     "cluster_101",
#     "cluster_110",
#     "cluster_200",
# ]

# cluster_names = [t for t in taxon_list if t.startswith("cluster")]

# for taxon in taxon_list:
#     characteristic_profiles[taxon] = np.zeros(n_genes, dtype=float)

# # create characteristic profiles for classes
# scratch = np.zeros(n_genes, dtype=float)
# scratch[:20] = np.linspace(0.0, 200.0, 20)
# characteristic_profiles["class_0"] += scratch
# characteristic_profiles["class_1"] += scratch

# scratch = np.zeros(n_genes, dtype=float)
# scratch[20:30] = np.linspace(0.0, 175.0, 10)[-1::-1]
# characteristic_profiles["class_0"] += scratch

# scratch = np.zeros(n_genes, dtype=float)
# scratch[10:30] = 2.0 * (np.arange(10, 30, dtype=float) - 20) ** 2
# characteristic_profiles["class_2"] += scratch

# # create characteristic profiles for subclasses

# vec0 = np.zeros(n_genes, dtype=float)
# vec0[30:60] = np.linspace(0.0, 200.0, 30)

# vec1 = np.zeros(n_genes, dtype=float)
# vec1[30:60] = np.linspace(0, 200.0, 30)[-1::-1]

# for taxon in taxon_list:
#     if "subclass" not in taxon:
#         continue
#     tag = taxon.split("_")[-1]
#     class_tag = tag[0]
#     characteristic_profiles[taxon] += characteristic_profiles[f"class_{class_tag}"]
#     subclass_tag = tag[-1]
#     if subclass_tag == "0":
#         characteristic_profiles[taxon] += vec0
#     else:
#         assert subclass_tag == "1"
#         characteristic_profiles[taxon] += vec1

#     params = rng.choice(np.arange(1, 30, dtype=int), 2)
#     wave = 0.5 * (1.0 + np.sin(2.0 * np.pi * (np.arange(30) + params[1]) / params[0]))
#     characteristic_profiles[taxon][30:60] *= wave

# # create characteristic profiles for clusters

# for taxon in taxon_list:
#     if "cluster" not in taxon:
#         continue

#     params = 1.0 + rng.random(2) * 6.0
#     wave = 100.0 * (1.0 + np.sin(2.0 * np.pi * np.arange(15) / params[0] + params[1]))

#     tag = taxon.split("_")[-1]
#     subclass_tag = tag[:2]
#     characteristic_profiles[taxon] += characteristic_profiles[
#         f"subclass_{subclass_tag}"
#     ]
#     cluster_tag = tag[-1]
#     if cluster_tag == "0":
#         characteristic_profiles[taxon][60:75] += wave
#     else:
#         assert cluster_tag == "1"
#         characteristic_profiles[taxon][75:90] += wave

#     params = rng.choice(np.arange(1, 10, dtype=int), 2)

# # %%
# cluster_name_list[0], cluster_name_list[1]

# # %% [markdown]
# # The cell below defines a helper function for generating noisy gene expression profiles drawn from our taxonomy.
# # 
# # **Again, you do not need to understand what this does in detail if you do not want to. It is not designed to be a realistic simulation.**

# # %%
# def create_noisy_cell(taxon_to_wgt, rng, n_dropouts=10):
#     """
#     A helper function to create a cell which is a member of one of our cartoon taxons, but with a noisy gene expression profile.

#     Parameters
#     ----------
#     taxon_to_wgt:
#         A dict mapping cell type taxon to the numeric weight multiplied by the taxon's characteristic gene expression profile
#         when assembling our simulated cell's gene expression profile (i.e. the simulated cell will be a linear combination
#         of characteristic gene expression profiles).
#     rng:
#         a numpy random number generator
#     n_dropouts:
#         select this many genes to randomly drop out of the cell

#     Returns
#     -------
#     A numpy array of raw gene expression (ints) for the simulated cell.
#     """

#     # assemble the linear combination of characteristic gene expression profiles
#     norm = 0.0
#     this_cell = np.zeros(n_genes, dtype=float)
#     for taxon in taxon_to_wgt:
#         this_cell += taxon_to_wgt[taxon] * characteristic_profiles[taxon]
#         norm += taxon_to_wgt[taxon]
#     this_cell = this_cell / norm

#     # for each gene, draw from a Poisson distribution centered on the value resulting
#     # from the linear combination above
#     nonzero = np.where(this_cell > 0.0)
#     this_cell[nonzero] = rng.poisson(this_cell[nonzero])

#     # set dropout genes to 0.0
#     if n_dropouts > 0:
#         sorted_dex = np.argsort(this_cell)
#         valid_for_dropout = np.ones(len(this_cell), dtype=bool)
#         valid_for_dropout[this_cell < 0.5] = False

#         if valid_for_dropout.sum() > n_dropouts:
#             dropout_genes = rng.choice(
#                 np.where(valid_for_dropout)[0], n_dropouts, replace=False
#             )
#             this_cell[dropout_genes] = 0.0

#     # clip negative values and convert to integers
#     this_cell = np.clip(this_cell, a_min=0.0, a_max=None)
#     this_cell = np.round(this_cell).astype(int)

#     return this_cell

# # %% [markdown]
# # The cell below creates an h5ad file containing the reference data we will use for our mapping.
# # 
# # **Just run this cell. We will look at its output below.**

# # %%
# rng = np.random.default_rng(6611222)

# n_labeled_cells = 2000

# var = pd.DataFrame([{"gene_id": f"g{ii}"} for ii in range(n_genes)]).set_index(
#     "gene_id"
# )

# obs_data = []
# cell_by_gene = np.zeros((n_labeled_cells, n_genes))
# for i_cell in range(n_labeled_cells):
#     # chosen_cluster = rng.choice(cluster_names)

#     cluster_name_list = rng.choice(cluster_names, 2, replace=True)
#     chosen_cluster = cluster_name_list[0]
#     other_cluster = cluster_name_list[1]

#     if other_cluster == chosen_cluster:
#         taxon_to_wgt = {chosen_cluster: 1.0}
#     else:
#         if chosen_cluster[:-1] == other_cluster[:-1]:
#             a0 = rng.normal(0.8, 0.05)
#         elif chosen_cluster[:-2] == other_cluster[:-2]:
#             a0 = rng.normal(0.9, 0.05)
#         else:
#             a0 = rng.normal(0.95, 0.02)

#         if a0 > 1.0:
#             a0 = 1.0
#         elif a0 < 0.0:
#             a0 = 0.0
#         a1 = 1.0 - a0

#         if a1 > a0:
#             _t = other_cluster
#             other_cluster = chosen_cluster
#             chosen_cluster = _t
#             _t = a1
#             a1 = a0
#             a1 = _t

#         taxon_to_wgt = {chosen_cluster: a0, other_cluster: a1}

#     this_cell = create_noisy_cell(taxon_to_wgt=taxon_to_wgt, rng=rng)

#     assert this_cell.sum() > 0
#     cell_by_gene[i_cell, :] = this_cell
#     tag = chosen_cluster.split("_")[-1]
#     obs_entry = {
#         "cell_id": f"cell_{i_cell}",
#         "cluster": chosen_cluster,
#         "subclass": f"subclass_{tag[:2]}",
#         "class": f"class_{tag[0]}",
#     }
#     obs_data.append(obs_entry)

# obs = pd.DataFrame(obs_data).set_index("cell_id")
# a_data = anndata.AnnData(X=cell_by_gene, obs=obs, var=var)


# training_data_path = tmp_dir / "training_data.h5ad"
# if training_data_path.exists():
#     training_data_path.unlink()
# a_data.write_h5ad(training_data_path)
# print(f"training_data is at {training_data_path}")

# # %%
# obs

# # %% [markdown]
# # We just created an h5ad file at `pipeline_example_data/training_data.h5ad`. This h5ad file contains everything we need to select marker genes with the `cell_type_mapper` code base. Let's examine the contents of that h5ad file.

# # %%
# src = anndata.read_h5ad("pipeline_example_data/training_data.h5ad", backed="r")
# print(src.shape)

# # %% [markdown]
# # The h5ad file is `n_cells` by `n_genes` in size. The `var` dataframe is just an index with the identifiers of the measured genes. Note that, when we finally do the cell type mapping, the code will directly compare these gene identifiers with the index of the `var` dataframe in the unlabeled data. If, for instance, the reference data identifies genes with ENSEMBL IDs and the unlabeled data identifies genes with NCBI IDs, the mapping code will determine that the two datasets have non-overlapping sets of genes, and the mapping will fail.

# # %%
# src.var

# # %%
# src.obs

# # %% [markdown]
# # The `obs` dataframe contains a unique label for each cell as well as its assigned cell type at each level in our taxonomy.

# # %% [markdown]
# # # 2 -- Select marker genes from our cartoon reference data
# # 
# # The entire purpose of this section of the notebook is to describe how to produce a JSON lookup table of marker genes from the cartoon reference data we just created. Marker gene selection is an ongoing research question, if you have a preferred way of selecting marker genes, all you have to do is encode your marker genes in a JSON file that comports with the schema documented [on this page](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/input_data_files/marker_gene_lookup.md) and you will be able to run your markers through the `cell_type_mapper`. This section of the notebook is only useful if you want to use the tools in the `cell_type_mapper` to do the marker gene selection for you.
# # 
# # ## Selecting marker genes using command line tools
# # 
# # In the cells below we will invoke the command line tools to create a set of marker genes for mapping to our reference dataset. The `!` at the start of each cell tells Jupyter that the cell is a command to be executed in the shell, i.e. these shells contain the commands you should invoke from the command line to select marker genes.
# # 
# # **Note:** to see the full call signatures of any of the tools below, replace the command line arguments with `--help`, i.e.
# # ```
# # python -m cell_type_mapper.cli.reference_markers --help
# # ```
# # will show documentation for the full call signature for the `reference_markers` tool. Most of the command line arguments for the tools in question have sensible defaults. We will explain the arguments we set to non-default values below.
# # 
# # ### Precompute statistics
# # 
# # The first step is to create the precomputed statisitcs file for the reference data. This is an HDF5 file that contains, essentially, the mean gene expression profile of every cell type in the taxonomy. It also contains an encoding of the cell type taxonomy. The detailed contents of this file are documented [on this page](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/input_data_files/precomputed_stats_file.md).
# # 
# # The call signature for the tool used below can be accessed by running `python -m cell_type_mapper.cli.precompute_stats_scrattch --help`. Most of the arguments are set to sensible defaults. The arguments we are explicitly setting below are:
# # 
# # - `heirarchy`: a list specifying the inheritance order from most gross to most fine of the levels in our cell type taxonomy (note the nested set of quotations around the cell types and the square brackets). These must be the names of columns in the `obs` dataframe of the reference data.
# # - `h5ad_path`: the path to the h5ad file containing the reference data for the taxonomy
# # - `n_processors`: the number of independent worker processes to spin up
# # - `normalization`: either `raw` or `log2CPM` indicating whether the data in the input h5ad file is in raw counts or `log2(CPM+1)` normalization.
# # - `clobber`: if False and the file specified by `output_path` already exists, the code will fail before doing anything. If `True`, and the output file already exists, the output file willb e overwritten
# # - `output_path`: the path to the HDF5 file that will be written out
# # 

# # %%
# !python -m cell_type_mapper.cli.precompute_stats_scrattch \
# --hierarchy '["class", "subclass", "cluster"]' \
# --h5ad_path pipeline_example_data/training_data.h5ad \
# --n_processors 3 \
# --normalization raw \
# --clobber True \
# --output_path pipeline_example_data/precomputed_stats.h5 

# # %% [markdown]
# # ### Select reference markers
# # 
# # The next step is to select "reference markers." This is essentially identifying all of the marker genes between all of the pairs of cell types that exist in our taxonomy. This step can be expensive (taking several hours for a taxonomy with several thousand taxons in it) and creates a very large file (10s of GB for a taxonomy with several thousand taxons in it). However, in theory, it only needs to be run once for a specific reference dataset.
# # 
# # A later step, query marker selection, downselects this set of marker genes to find an optimal and minial combination of marker genes given the unlabeled dataset at hand. That step takes the output of this step as an input.
# # 
# # One quirk of the `reference_marker` tool is that it does not let you specify an output file, only an output directory. This design is meant to support use cases where multiple experimental apparatuses (e.g. 10Xv3, 10Xv2, 10XMulti) are used to take data for the same reference dataset. The idea is that the data from each apparatus would be treated indpendently through the `precomputed_stats` and `reference_markers` steps, and then joined during query marker selection (the next step) into a single set of query markers. This level of detail is likely overkill for most users. The name of the file written by this tool is printed out to `stdout` during running. Look for a file named like `reference_markers.h5` in your specified output directory. For those worried about overwriting previous results, if you do not run with `--clobber True`, the tool will crash rather than overwrite an existing file. 
# # 
# # There are many unexercised arguments to this function (again, because they default to reasonable values) in the example below. These control the statistical definition of what counts as a marker gene. They are documented [on this page](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/algorithms/marker_gene_selection.md). The arguments we are specifying below are
# # - `precomputed_path_list`: a list of paths to precomputed statistics files (i.e. the output of `precompute_stats_scrattch` above). For each of these files, a different reference marker file will be written out for reasons described above (multiple apparatuses; one taxonomy). Note the nested quotation marks surrounding the list.
# # - `n_valid`: the code has strict criteria which a gene must meet to be considered a marker for discriminating between two cell types. Some data is so noisy, however, that, for some pairs of cell types, very few genes will satisfy these strict criteria. [The algorithm](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/algorithms/marker_gene_selection.md) will gradually relax the criteria so that any given pair of cell types has at least `n_valid` markers discriminating between them.
# # - `output_dir`: the directory to which the reference marker files will be written (see discussion above for why the output destination is specified this way)
# # - `clobber`: if running the code would overwrite an existing file, the code will fail before doing anything. Set `--clobber True` to go ahead and overwrite te file anyway.
# # 
# # There is one argument we did not set that we would like to draw attention to
# # - `query_path`: the path to the h5ad file specifying the unlabeled data that will ultimately be mapped. If you specify this argument, then the marker genes will only be selected from those genes that are present in the `var` dataframe of the corresponding file. If you do not specify this argument, then any genes in the precomputed stats file are considered valid choices as marker genes. Specifying the `query_path` can be helpful in preventing the code from placing too much stock in a very distinct set of marker genes that end up not being useful because they are not present in the unlabeled data.

# # %%
# !python -m cell_type_mapper.cli.reference_markers \
# --precomputed_path_list "['pipeline_example_data/precomputed_stats.h5']" \
# --n_valid 20 \
# --n_processors 3 \
# --output_dir pipeline_example_data \
# --clobber True

# # %% [markdown]
# # ### Select query markers
# # 
# # The final step in marker selection is to downsample the `reference_markers` file to find a minimal set of markers to be used during the cell type mapping. This is a much quicker step than `reference_marker` finding, and produces a JSON file that maps between cell types in our taxonomy tree and marker genes that are used to discriminate between those cell types' children as discussed [on this page](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/input_data_files/marker_gene_lookup.md).
# # 
# # The arguments we are using below are
# # - `reference_marker_path_list`: the list of the paths of the reference marker files (i.e. the HDF5 files produced by the `reference_marker` tool above) from which these query markers are to be sampled. For the same reason that the reference marker finder specifies only an output directory, this function takes a list of reference marker files. This is the step where multiple reference markers for multiple experimental apparatuses could theoretically be joined into a single marker gene lookup table.
# # - `n_per_utility`: for each pair of cell type taxons `A` and `B`, there are two "utilities" a marker gene can serve: "up regulated in A" and "up regulated in B." The query marker finder will do its best to select `n_per_utility` genes that are markers *and* are up regulated in A and `n_per_utility` genes that are markers *and* are up regulated in B (assuming that many markers exist in the reference marker file).
# # - `n_processors`: the number of independent worker processes to spin up while selecting markers.
# # 
# # There is two areguments we are not using below which it is worth taking note of:
# # - `search_for_stats_file`: technically, the query marker selecter needs the reference marker file and the precomputed stats file associated with it. The path to the precomputed stats file used in generating the reference marker file is encoded in the metadata stored in the reference marker file. If, for some reason, you files have moved from their original locations (for instance, if you are working in a cloud computing environment), you can set `--search_for_stats_file True` and the query marker finder will look for the precomputed stats file in the same directory that contains the reference marker file.
# # - `query_path`: the path to the h5ad file specifying the unlabeled data that will ultimately be mapped. If you specify this argument, then the marker genes will only be selected from those genes that are present in the `var` dataframe of the corresponding file. If you do not specify this argument, then any genes in the reference marker file are considered valid choices as marker genes.

# # %%
# !python -m cell_type_mapper.cli.query_markers \
# --output_path pipeline_example_data/query_markers.json \
# --reference_marker_path_list '["pipeline_example_data/reference_markers.h5"]' \
# --n_per_utility 10 \
# --n_processors 3

# # %% [markdown]
# # ## Selecting marker genes programmatically
# # 
# # In the cell below, we will demonstrate how to run the same commands from within a python script, for those who would rather work directly in python than invoke command-line tools.

# # %%
# from cell_type_mapper.cli.precompute_stats_scrattch import PrecomputationScrattchRunner

# from cell_type_mapper.cli.reference_markers import ReferenceMarkerRunner

# from cell_type_mapper.cli.query_markers import QueryMarkerRunner

# precomputation_config = {
#     "hierarchy": ["class", "subclass", "cluster"],
#     "h5ad_path": "pipeline_example_data/training_data.h5ad",
#     "output_path": "pipeline_example_data/precomputed_stats.h5ad",
#     "n_processors": 3,
#     "normalization": "raw",
#     "clobber": True,
# }

# # the args=[] is important to prevent the PrecomputationScrattchRunner from grabbing
# # any command line arguments when you invoke your python script
# precomputation_runner = PrecomputationScrattchRunner(
#     args=[], input_data=precomputation_config
# )

# precomputation_runner.run()

# print("\n====done with precomputed_stats====\n")

# reference_config = {
#     "precomputed_path_list": ["pipeline_example_data/precomputed_stats.h5"],
#     "n_valid": 20,
#     "output_dir": "pipeline_example_data",
#     "clobber": True,
# }

# reference_runner = ReferenceMarkerRunner(args=[], input_data=reference_config)

# reference_runner.run()

# print("\n====done with reference_markers====\n")


# query_config = {
#     "output_path": "pipeline_example_data/query_markers.json",
#     "reference_marker_path_list": ["pipeline_example_data/reference_markers.h5"],
#     "n_per_utility": 10,
#     "n_processors": 3,
# }

# query_runner = QueryMarkerRunner(args=[], input_data=query_config)

# query_runner.run()

# # %% [markdown]
# # # 3 -- Create unlabeled data
# # 
# # Let's create a small h5ad file of unlabeled data. For each cluster, we will create one cell which is a linear combination of the chosen cluster with a small (0.05) contribution from another cluster. This will introduce some noise into the mapping an allow us to examine the ways in which the `cell_type_mapper` encodes uncertainty in its mappings. Since we are only creating one cartoon cell per cluster, we will identify each cell according to the cluster we know it "should" belong to.

# # %%
# rng = np.random.default_rng(77121133)

# n_cells = 12
# obs_data = []
# cell_by_gene = []
# for cluster in cluster_names:
#     this_obs = {"cell_id": cluster}
#     obs_data.append(this_obs)
#     if cluster != "cluster_200":
#         other_cluster = "cluster_200"
#     else:
#         other_cluster = "cluster_000"

#     this_cell = create_noisy_cell(
#         taxon_to_wgt={cluster: 0.95, other_cluster: 0.05}, rng=rng
#     )
#     cell_by_gene.append(this_cell)

# cell_by_gene = np.array(cell_by_gene)
# obs = pd.DataFrame(obs_data).set_index("cell_id")

# a_data = anndata.AnnData(obs=obs, var=var, X=cell_by_gene)

# unlabeled_path = pathlib.Path("pipeline_example_data/unlabeled.h5ad")
# if unlabeled_path.exists():
#     unlabeled_path.unlink()
# a_data.write_h5ad(unlabeled_path)

# # %% [markdown]
# # # 4 -- Map the unlabeled data using the command line tool
# # 
# # We wrote our unlabeled data to `pipeline_example/unlabeled.h5ad`.
# # 
# # Below is the command to perform the cell type mapping using the command line tool. The arguments for that tool are documented [on this page](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/mapping_cells.md#mapping-unlabeled-data-onto-the-reference-taxonomy). The arguments we are using below are
# # - `precomputed_stats.path`: the path to the HDF5 file containing the precomputed statistics for our taxonomy (i.e. the output of the `precompute_stats_scrattch` tool above)
# # - `query_markers.serialized_lookup`: the path to the JSON file containing the marker gene lookup table (i.e. the output of the `query_markers` tool above)
# # - `type_assignment.bootstrap_factor`: the factor by which to randomly downsample the full marker gene list at each [maapping iteration](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/algorithms/hierarchical_mapping.md)
# # - `type_assignment.bootstrap_iteration`: the number of mapping iterations to perform at each level of the taxonomy
# # - `type_assignment.rng_seed`: the seed for the random number generator. This is specified as an argument so that two runs with the same configuration produce identical results (rather than depending on the vagaries of a random number generator seeded from the system clock)
# # - `type_assignment.n_processors`: the number of independent worker processes to spin up while mapping.
# # - `query_path`: the path to the h5ad file containing the unlabeled data (i.e. the cell-by-gene data for the cells you would like to map onto this taxonomy).
# # - `extended_result_path`: the path to the JSON file where mapping results will be written.
# # - `csv_result_path`: the path to the CSV file where mapping results will be written.

# # %%
# !python -m cell_type_mapper.cli.from_specified_markers \
# --precomputed_stats.path pipeline_example_data/precomputed_stats.h5 \
# --query_markers.serialized_lookup pipeline_example_data/query_markers.json \
# --type_assignment.bootstrap_factor 0.5 \
# --type_assignment.bootstrap_iteration 100 \
# --type_assignment.rng_seed 661123 \
# --type_assignment.n_processors 1 \
# --type_assignment.normalization raw \
# --query_path pipeline_example_data/unlabeled.h5ad \
# --extended_result_path pipeline_example_data/mapping_output.json \
# --csv_result_path pipeline_example_data/mapping_output.csv


# # %% [markdown]
# # # 4.1 -- Map the unlabeled dataset programmatically
# # 
# # Here we show how to perform the same mapping from within a python script

# # %%
# from cell_type_mapper.cli.from_specified_markers import FromSpecifiedMarkersRunner

# config = {
#     "precomputed_stats": {"path": "pipeline_example_data/precomputed_stats.h5"},
#     "query_markers": {"serialized_lookup": "pipeline_example_data/query_markers.json"},
#     "type_assignment": {
#         "bootstrap_factor": 0.5,
#         "bootstrap_iteration": 100,
#         "rng_seed": 661123,
#         "n_processors": 1,
#         "normalization": "raw",
#     },
#     "query_path": "pipeline_example_data/unlabeled.h5ad",
#     "extended_result_path": "pipeline_example_data/mapping_output.json",
#     "csv_result_path": "pipeline_example_data/mapping_output.csv",
# }

# mapping_runner = FromSpecifiedMarkersRunner(args=[], input_data=config)

# mapping_runner.run()

# # %% [markdown]
# # # 5 -- Examining the outputs of the mapping
# # 
# # We wrote the outputs of the mapping in two forms, a CSV file at `pipeline_example_data/mapping_output.csv` and a JSON file at `pipeline_example_data/mapping_output.json`. The contents of these files are documented [on this page](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/output.md).
# # 
# # ## The CSV output file
# # 
# # The CSV file is the simplest output file, so we will examine that first.
# # 
# # **Note:** most of the work below is redundant with what is discussed in [this notebook](https://github.com/AllenInstitute/cell_type_mapper/blob/main/examples/explore_mapping_results.ipynb), which creates an h5ad file of real data and shows users how to run it through the on-line MapMyCells interface, download the results, and inspect them. Here, we will go into a little more discussion about the probabilistic nature of the cell type mapping, since we had direct control over both the cell type taxonomy and the unlabeled data.
# # 
# # The CSV file can be read in with pands, provided you instruct pandas to ignore all lines that start with `#`.

# # %%
# results_df = pd.read_csv("pipeline_example_data/mapping_output.csv", comment="#")

# # %%
# results_df

# # %% [markdown]
# # Inspecting this CSV we see that
# # 
# # - Every cell was ultimately assigned to the expected cluster (i.e. `cell_id` is identical to `cluster_label`.
# # - In some cases, the `cell_type_mapper` was not absolutely confident in its assignments (`bootstrapping_probability < 1.0`), this being a reflection of the noise we introduced by making each cell a linear combination of two cell type characteristic profiles.

# # %% [markdown]
# # ## The JSON output file

# # %%
# with open("pipeline_example_data/mapping_output.json", "rb") as src:
#     json_results = json.load(src)

# # %% [markdown]
# # As discussed [in this notebook](https://github.com/AllenInstitute/cell_type_mapper/blob/main/examples/explore_mapping_results.ipynb), the JSON output file contains a great deal of metadata about the specific `cell_type_mapper` run. You should consult that notebook for a guided tour of those metadata products. Here we will focus on the actual cell type mapping results contained under the `'results'` key. Let's look at one cell in particular.

# # %%
# print(json.dumps(json_results["results"][2], indent=2))

# # %% [markdown]
# # The cell is represented by a dict whose keys are the levels in our cartoon cell type taxonomy (except for `cell_id`, which maps to the unique identifier of this cell). For each level of the taxonomy, we see have
# # - The cell type to which the cell was actually assigned
# # - The fraction of mapping iterations that made that assignment (the `bootstrapping_probability`)
# # - The average correlation coefficient between the cell and the assigned cell type (the `avg_correlation`)
# # - A list of the cell types to which the cell might have been assigned but wasn't, along with their `bootstrapping_probability` and `avg_correlation` scores.
# # 
# # For instance, looking at the `class` level for our example cell

# # %%
# print(json.dumps(json_results["results"][2]["class"], indent=2))

# # %% [markdown]
# # We see that the cell was assigned to `class_0`, that `0.85` of the random mapping iterations made that determination, and that the average correlation coefficient between the cell and the average gene expression profile for `class_0` in the space of marker genes was `0.672`. However, `0.11` of the random iterations chose `class_2` (average correlation coefficient `0.582`) and `0.04` of the iterations chose `class_1` (average correlation coefficient `0.614`). These latter figures can be taken as an indication of how trustworthy the mapping was.
# # 
# # The `aggregate_probability` field is a running product of the `bootstrapping_probability` of cell at the given level in the taxonomy, i.e. because our taxonomy is structured `class -> subclass -> cluster`, the `aggregate_probability` at the subclass level is the product of the `bootstrapping_probability` at the subclass level with the `bootstrapping_probability` at the class level. The `aggregate_probability` at the cluster level is the product of `bootstrapping_probability` at all three levels in the taxonomy. A way to think about it is that the `bootstrapping_probability` at a given level is a conditional probability (conditioned on the assumption that the prior cell type assignments were correct). The `aggregate_probability` is the unconditional probability that the cell is a member of the given cell type.
# # 
# # Finally `directly_assigned` is a boolean indicating whether or not the code actually mapped to that level in the taxonomy. The cell type mapper gives you the option to skip over levels of the taxonomy, in which case those levels are inferred from the assignments made at lower levels in the taxonomy. Let us explore what this means by using the `--flatten` command line argument, which tells the cell type mapper to skip all intermediate levels in the taxonomy and assign the cell directly to the leaf level (`cluster` in our taxonomy).

# # %%
# !python -m cell_type_mapper.cli.from_specified_markers \
# --precomputed_stats.path pipeline_example_data/precomputed_stats.h5 \
# --query_markers.serialized_lookup pipeline_example_data/query_markers.json \
# --type_assignment.bootstrap_factor 0.5 \
# --type_assignment.bootstrap_iteration 100 \
# --type_assignment.rng_seed 661123 \
# --type_assignment.n_processors 1 \
# --type_assignment.normalization raw \
# --query_path pipeline_example_data/unlabeled.h5ad \
# --extended_result_path pipeline_example_data/flat_mapping_output.json \
# --csv_result_path pipeline_example_data/flat_mapping_output.csv \
# --flatten True

# # %%
# with open("pipeline_example_data/flat_mapping_output.json", "rb") as src:
#     flat_json_results = json.load(src)

# print(json.dumps(flat_json_results["results"][2], indent=2))

# # %% [markdown]
# # In the case of mapping with `--flatten True`, the mapper directly assigns the cell to the cluster level and infers the cell's subclass and class membership based on which subclasses and classes are ancestors of the assigned cluster. You can tell that this is what happened because `directly_assigned` is `false` in the class and subclass levels of the assignment. For book keeping purposes, `bootstrapping_probabiltiy`, `avg_correlation` and `aggregate_probability` for the class and subclass levels are copied directly from the cluster level, even though no bootstrapping iteration (or correlation with marker genes) was done at these levels.
# # 
# # Let us look at one more cell in our original, hierarchical mapping example.

# # %%
# print(json.dumps(json_results["results"][7], indent=2))

# # %% [markdown]
# # This cell was assigned to `class_2`. You will recall that only one subclass (`subclass_20`) and only one cluster (`cluster_200`) descend from `class_2`. In this case, though subclass and cluster were both directly assigned, the `bootstrapping_probability` at those levels is `1.0`. At the point where the cell was assigned to `class_2`, it **must** also belong to `subclass_20` and `cluster_200` according to our cell type taxonomy. No actual correlation was done. `bootstrapping_probability` was set to 1.0 and `avg_correlation` for the subclass and cluster were copied from the `avg_correlation` at the class level.

# # %% [markdown]
# # # 6 -- Special case: small taxonomies
# # 
# # We broke up the mapping and marker gene selection steps above because marker gene selection (especially reference marker selection) can be time consuming for large taxonomies. The code must compare all combinations of leaf nodes. For taxonomies with a few thousand nodes, this process can take hours. However, it only needs to be run once (in theory) for any given taxonomy, and the results can be resuded for different unlabeled datasets. For smaller taxonomies (taxonomies with a few hundred nodes) the marker gene selection process can be accomplished in minutes and there is no harm in re-running it for each unlabeled dataset (depending on your tolerance for "wasting" a few minutes). The advantage to this approach is that, if you end up mapping several datasets, each with a different set of measured genes in it, you can be certain that the marker genes will be selected with prior knowledge of the genes that are available in your unlabeled dataset and the code will not depend on the presence of a very good marker that is sadly unavailable in the unlabeled dataset.
# # 
# # We have provided a command line tool for selecting the marker genes "on the fly" and mapping with those internally selected marker genes. Note that, for this tool, all of the parameters passed to the reference and query marker finder are passed to the single tool with the `query_markers` and `reference_markers` prefix as appropriate.

# # %%
# !python -m cell_type_mapper.cli.map_to_on_the_fly_markers \
# --precomputed_stats.path pipeline_example_data/precomputed_stats.h5 \
# --type_assignment.bootstrap_factor 0.5 \
# --type_assignment.bootstrap_iteration 100 \
# --type_assignment.rng_seed 661123 \
# --type_assignment.normalization raw \
# --query_path pipeline_example_data/unlabeled.h5ad \
# --extended_result_path pipeline_example_data/otf_mapping_output.json \
# --csv_result_path pipeline_example_data/otf_mapping_output.csv \
# --query_markers.n_per_utility 10 \
# --reference_markers.n_valid 20 \
# --n_processors 1

# # %% [markdown]
# # Let's load the results of our "on the fly marker mapping" and show that they are identical to the results we got when we split up marker selection and mapping into different steps.

# # %%
# with open("pipeline_example_data/otf_mapping_output.json", "rb") as src:
#     otf_json_results = json.load(src)


# assert otf_json_results["results"] == json_results["results"]

# # %% [markdown]
# # Here is how to call the same on-the-fly mapping programmatically

# # %%
# from cell_type_mapper.cli.map_to_on_the_fly_markers import OnTheFlyMapper


# config = {
#     "precomputed_stats": {"path": "pipeline_example_data/precomputed_stats.h5"},
#     "type_assignment": {
#         "bootstrap_factor": 0.5,
#         "bootstrap_iteration": 100,
#         "rng_seed": 661123,
#         "normalization": "raw",
#     },
#     "reference_markers": {"n_valid": 20},
#     "query_markers": {"n_per_utility": 10},
#     "query_path": "pipeline_example_data/unlabeled.h5ad",
#     "extended_result_path": "pipeline_example_data/otf_mapping_output.json",
#     "csv_result_path": "pipeline_example_data/otf_mapping_output.csv",
#     "n_processors": 1,
# }


# runner = OnTheFlyMapper(args=[], input_data=config)
# runner.run()

# # %% [markdown]
# # One potential disadvantage to this approach is that there is less of a "paper trail." The reference marker and query marker files are not written out to disk. Only the final mappings are written to disk. However, you can find the marker gene lookup table that was actually used to map your data in the JSON results file here:

# # %%
# otf_json_results["marker_genes"]

# # %%



