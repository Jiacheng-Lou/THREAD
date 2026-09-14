DPM AD example dataset

Generated: 2026-07-18 10:38:45 UTC
DPM package version: 0.1.0
Chromosome: 1
GWAS matching mode: snp_name
LD ancestry: EUR
LD symmetrization: reciprocal pair recovery
PSD handling: trace-preserving eigenvalue clipping

Genes: 200
Cellular labels: 31
Unique SNPs: 20605
Directed LD records: 364382

Files:
  example_x.rds
    Gene-by-cellular-label expression matrix.
  example_response.rds
    Frozen LD-corrected gene-level response.
  example_gwas.tsv.gz
    GWAS summary-statistics subset.
  example_gene_snp_matches.rds
    SNP-to-gene mapping subset.
  example_ld_chr1.rds
    Chromosome-specific LD subset containing only required SNP pairs.
  example_data_manifest.csv
    File sizes and MD5 checksums.

These bundled files are small demonstration and test fixtures prepared by the
DPM authors. The example follows the project's simulation-data generation
workflow and does not redistribute the manuscript's real disease datasets.
The gene-SNP mapping was generated for this example, and
example_response.rds is a frozen derived output from the example workflow.
example_ld_chr1.rds is a reduced chromosome-1 LD fixture derived from the
1000 Genomes Project Phase 3 European (EUR) reference panel. THREAD does not
distribute the full genome-wide LD reference panel; users performing real
analyses must provide a compatible ancestry-matched external LD reference.
