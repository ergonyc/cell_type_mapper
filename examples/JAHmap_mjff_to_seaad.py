from cell_type_mapper.cli.map_to_on_the_fly_markers import OnTheFlyMapper


def main():

    config = {
        "query_path": "/Users/ergonyc/Projects/ASAP/data/asap-curated-cohort-pmdbs-sc-rnaseq/pmdbs_sc_rnaseq/cohort_analysis/asap-cohort.merged_adata_object.h5ad",
        "tmp_dir": "/Users/ergonyc/Projects/ASAP/cell_type_mapper/tmp",
        "extended_result_path": "/Users/ergonyc/Projects/ASAP/cell_type_mapper/examples/data/asap-cohort.merged_adata_object.mmc.seaad_results.20250123b.json",
        "csv_result_path": "/Users/ergonyc/Projects/ASAP/cell_type_mapper/examples/data/asap-cohort.merged_adata_object.mmc.seaad_results.20250123b.csv",
        "log_path": "/Users/ergonyc/Projects/ASAP/cell_type_mapper/examples/data/asap-cohort.merged_adata_object.mmc.seaad_log.20250123b.txt",
        "cloud_safe": False,
        "verbose_csv": True,
        "n_processors": 4,
        "max_gb": 32.0,
        "map_to_ensembl": True,
        "type_assignment": {
            "normalization": "raw",
            "bootstrap_iteration": 100,
            "bootstrap_factor": 0.5,
            "chunk_size": 40000,
            "n_runners_up": 5,
            "rng_seed": 11235813,
        },
        "precomputed_stats": {
            "path": "/Users/ergonyc/Projects/ASAP/cell_type_mapper/precomputed_stats.20231120.sea_ad.MTG.h5"
        },
        "reference_markers": {"log2_fold_min_th": 0.5},
        "query_markers": {"n_per_utility": 15, "genes_at_a_time": 1},
    }

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



