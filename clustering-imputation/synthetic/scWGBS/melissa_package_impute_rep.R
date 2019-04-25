# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()) {cur.dir <- dirname(parent.frame(2)$ofile); setwd(cur.dir)}
suppressPackageStartupMessages(library(Melissa))

io <- list()
io$met_dir <- "~/datasets/encode/scWGBS/parsed/binarised/"
io$anno_file <- "~/datasets/encode/annotations/promoter_hg38_chr1_6.bed"
melissa_data <- create_melissa_data_obj(met_dir = io$met_dir, anno_file = io$anno_file,
                                        cov = 3, no_cores = 3)

saveRDS(file = "~/datasets/encode/scWGBS/processed/unfiltered/prom10k_package.rds", melissa_data)

melissa_data <- filter_by_cpg_coverage(melissa_data, min_cpgcov = 15)
melissa_data <- filter_by_variability(melissa_data, min_var = 0.2)
melissa_data <- filter_by_coverage_across_cells(melissa_data, min_cell_cov_prcg = 0.5)

# Partition to training and test set
part_dt <- partition_dataset(dt_obj = melissa_data, data_train_prcg = 0.4,
                             cpg_train_prcg = 0.5, is_synth = FALSE)
# Run Melissa with K = 3 clusters
K <- 3
melissa_obj <- melissa(X = part_dt$met, K = K, basis = create_rbf_object(M = 9),
                       delta_0 = rep(3, K) + rbeta(K, 1e-1, 1e1),
                       vb_max_iter = 500, is_parallel = TRUE, no_cores = 10)
