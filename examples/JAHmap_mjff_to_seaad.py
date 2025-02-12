from cell_type_mapper.cli.map_to_on_the_fly_markers import OnTheFlyMapper
from pathlib import Path
import os


os.environ["AIBS_BKP_USE_TORCH"] = "false"

# os.environ["NUMEXPR_NUM_THREADS"] = "1"
# os.environ["MKL_NUM_THREADS"] = "1"
# os.environ["OMP_NUM_THREADS"] = "1"

HOME = Path.home()

ROOT = HOME / "Projects/ASAP/cell_type_mapper"


TMP_DIR = ROOT / "tmp"
if not TMP_DIR.exists():
    TMP_DIR.mkdir()

ASAP_DATA = (
    HOME
    / "Projects/ASAP/data/asap-curated-cohort-pmdbs-sc-rnaseq/pmdbs_sc_rnaseq/cohort_analysis"
)
ASAP_DATA = HOME / "Projects/ASAP/harmonized-wf-dev/data/cohort_analysis"

CHUNK_SIZE = 40000
N_RUNNERS_UP = 5
RNG_SEED = 11235813
N_PROCESSORS = 8
MAX_GB = 48.0


#
FILE_ROOT = "asap-cohort.merged_adata_object"
DATE = "20250129"


REFERENCE = "SEAAD"

EXTENDED_RESULTS = f"{FILE_ROOT}.mmc.{REFERENCE}.{DATE}.json"
LOG_FILE = f"{FILE_ROOT}.mmc.{REFERENCE}_log.{DATE}.txt"
CSV_RESULTS = f"{FILE_ROOT}.mmc.{REFERENCE}_results.{DATE}.csv"


RESULTS_DIR = ROOT / f"MMC.{REFERENCE}_RESULTS"
if not RESULTS_DIR.exists():
    RESULTS_DIR.mkdir()

PRECOMPUTED_STATS = (
    f"examples/data/abc_atlas_data/precomputed_stats.20231120.sea_ad.MTG.h5"
)


def main():

    config = {
        # "query_path": f"{ASAP_DATA}/{FILE_ROOT}.h5ad",
        "query_path": f"{TMP_DIR}/_{FILE_ROOT}.h5ad",
        "tmp_dir": f"{TMP_DIR}",
        "extended_result_path": f"{RESULTS_DIR / EXTENDED_RESULTS}",
        "csv_result_path": f"{RESULTS_DIR / CSV_RESULTS}",
        "log_path": f"{RESULTS_DIR / LOG_FILE}",
        "cloud_safe": False,
        "verbose_csv": True,
        "n_processors": N_PROCESSORS,
        "max_gb": MAX_GB,
        "map_to_ensembl": True,
        "type_assignment": {
            "normalization": "raw",
            "bootstrap_iteration": 100,
            "bootstrap_factor": 0.5,
            "chunk_size": CHUNK_SIZE,
            "n_runners_up": N_RUNNERS_UP,
            "rng_seed": RNG_SEED,
        },
        "precomputed_stats": {"path": f"{ROOT / PRECOMPUTED_STATS}"},
        "reference_markers": {"log2_fold_min_th": 0.5},
        "query_markers": {"n_per_utility": 15, "genes_at_a_time": 1},
    }

    # config = {
    #     "query_path": "/Users/ergonyc/Projects/ASAP/data/asap-curated-cohort-pmdbs-sc-rnaseq/pmdbs_sc_rnaseq/cohort_analysis/asap-cohort.merged_adata_object.h5ad",
    #     "tmp_dir": "/Users/ergonyc/Projects/ASAP/cell_type_mapper/tmp",
    #     "extended_result_path": "/Users/ergonyc/Projects/ASAP/cell_type_mapper/examples/data/asap-cohort.merged_adata_object.mmc.seaad_results.20250123b.json",
    #     "csv_result_path": "/Users/ergonyc/Projects/ASAP/cell_type_mapper/examples/data/asap-cohort.merged_adata_object.mmc.seaad_results.20250123b.csv",
    #     "log_path": "/Users/ergonyc/Projects/ASAP/cell_type_mapper/examples/data/asap-cohort.merged_adata_object.mmc.seaad_log.20250123b.txt",
    #     "cloud_safe": False,
    #     "verbose_csv": True,
    #     "n_processors": 4,
    #     "max_gb": 32.0,
    #     "map_to_ensembl": True,
    #     "type_assignment": {
    #         "normalization": "raw",
    #         "bootstrap_iteration": 100,
    #         "bootstrap_factor": 0.5,
    #         "chunk_size": 40000,
    #         "n_runners_up": 5,
    #         "rng_seed": 11235813,
    #     },
    #     "precomputed_stats": {
    #         "path": "/Users/ergonyc/Projects/ASAP/cell_type_mapper/precomputed_stats.20231120.sea_ad.MTG.h5"
    #     },
    #     "reference_markers": {"log2_fold_min_th": 0.5},
    #     "query_markers": {"n_per_utility": 15, "genes_at_a_time": 1},
    # }

    runner = OnTheFlyMapper(args=[], input_data=config)

    runner.run()


if __name__ == "__main__":
    main()


# # SEA-AD "references"
# MTG
# DLPFC
# MTG_DLPFC

# look at discrepancies.

# ## xylena reference...
# xylena_train for Taxonomy

# test on xylena_test + xylena_query
