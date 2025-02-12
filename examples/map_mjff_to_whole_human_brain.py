from cell_type_mapper.cli.from_specified_markers import (
    FromSpecifiedMarkersRunner
)


def main():
    config = {
      "query_path": "/allen/scratch/aibstemp/chris.morrison/MJFF_pmdbs-sc-rnaseq/asap-cohort.merged_adata_object.h5ad",
      "tmp_dir": "/local1/scott_daniel/scratch",
      "extended_result_path": "/allen/scratch/aibstemp/chris.morrison/MJFF_MapMyCells_results/asap-cohort.merged_adata_object.mmc_results.20250123.json",
      "csv_result_path": "/allen/scratch/aibstemp/chris.morrison/MJFF_MapMyCells_results/asap-cohort.merged_adata_object.mmc_results.20250123.csv",
      "log_path": "/allen/scratch/aibstemp/chris.morrison/MJFF_MapMyCells_results/asap-cohort.merged_adata_object.mmc_results.20250123.validation_log.txt",
      "map_to_ensembl": True,
      "cloud_safe": False,
      "verbose_csv": True,
      "type_assignment": {
        "normalization": "raw",
        "n_processors": 6,
        "bootstrap_iteration": 100,
        "bootstrap_factor": 0.5,
        "chunk_size": 10000,
        "n_runners_up": 5
      },
      "precomputed_stats": {
        "path": "/allen/aibs/technology/danielsf/knowledge_base/mmc_official_runs/data_files/whb/precomputed_stats.siletti.training.h5"
      },
      "query_markers": {
        "serialized_lookup": "/allen/aibs/technology/danielsf/knowledge_base/mmc_official_runs/data_files/whb/query_markers.n10.20240221800.json"
      }
    }

    runner = FromSpecifiedMarkersRunner(
        args=[],
        input_data=config
    )

    runner.run()


if __name__ == "__main__":
    main()
