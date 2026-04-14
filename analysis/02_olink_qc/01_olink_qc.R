# Olink QC and IQR Outlier Detection
# Pipeline for MS and ALS Proteomics Replications

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(stringr)
    library(here)
})

REPO_ROOT <- here::here()
CADASIL_DIR <- file.path(dirname(REPO_ROOT), "CADASIL_Proteome_ML_Keller_2024_Rebuttal")
DATA_DIR <- file.path(CADASIL_DIR, "data", "ukb")
OUT_DIR <- file.path(REPO_ROOT, "results", "qc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Data paths
OLINK_PATH <- file.path(DATA_DIR, "olink", "i0", "olink_instance_0_extracted_data.csv")
MS_COHORT_PATH <- file.path(REPO_ROOT, "data", "ukb", "cohort", "ms_cohort.csv")
ALS_COHORT_PATH <- file.path(REPO_ROOT, "data", "ukb", "cohort", "als_cohort.csv")

message("Loading Olink NPX Data (Instance 0)...")
if (!file.exists(OLINK_PATH)) {
    stop("Olink data not found at: ", OLINK_PATH)
}
olink_dt <- fread(OLINK_PATH)
# Rename columns nicely (remove olink_instance_0. prefix if present)
names(olink_dt) <- str_replace(names(olink_dt), "^olink_instance_0\\.", "")

message("Loading cohorts...")
ms_cohort <- fread(MS_COHORT_PATH)
als_cohort <- fread(ALS_COHORT_PATH)

run_qc <- function(disease_name, cohort_dt) {
    message(paste("Running QC for", disease_name, "cohort..."))
    
    # Merge cohort info with Olink
    dt <- merge(cohort_dt, olink_dt, by = "eid", all.x = FALSE, all.y = FALSE)
    n_total <- nrow(dt)
    
    if (n_total == 0) {
        message("  -> No participants with Olink data in ", disease_name, " cohort.")
        return(NULL)
    }
    
    disease_status_col <- paste0(tolower(disease_name), "_status")
    
    # Calculate IQR-based outliers on global NPX sum to replicate standard QC
    protein_cols <- setdiff(names(dt), names(cohort_dt))
    
    # Simple median per participant across all proteins
    dt[, global_median := rowMedians(as.matrix(.SD), na.rm=TRUE), .SDcols = protein_cols]
    
    Q1 <- quantile(dt$global_median, 0.25, na.rm=TRUE)
    Q3 <- quantile(dt$global_median, 0.75, na.rm=TRUE)
    IQR_val <- Q3 - Q1
    lower_bound <- Q1 - 3 * IQR_val
    upper_bound <- Q3 + 3 * IQR_val
    
    dt[, qc_outlier := global_median < lower_bound | global_median > upper_bound]
    n_outliers <- sum(dt$qc_outlier, na.rm=TRUE)
    
    message(sprintf("  -> Detected %d IQR outliers (3 IQR from median NPX).", n_outliers))
    
    # Plot global NPX distribution colored by case/control to ensure no batch effects
    p <- ggplot(dt, aes_string(x = disease_status_col, y = "global_median", fill = disease_status_col)) +
        geom_violin(alpha = 0.6) +
        geom_boxplot(width = 0.2, alpha = 0.8) +
        theme_minimal() +
        labs(title = paste(disease_name, "Global NPX Median"),
             y = "Median NPX Across All Proteins",
             x = "Status") +
        theme(legend.position = "none")
    
    plot_path <- file.path(OUT_DIR, paste0(tolower(disease_name), "_npx_distribution.pdf"))
    ggsave(plot_path, p, width = 6, height = 5)
    
    # Filter out outliers and save cleaned Olink subset for this cohort
    dt_clean <- dt[qc_outlier == FALSE]
    dt_clean[, c("global_median", "qc_outlier") := NULL]
    
    clean_path <- file.path(REPO_ROOT, "data", "ukb", "olink", "processed", 
                            paste0(tolower(disease_name), "_olink_qc.csv"))
    dir.create(dirname(clean_path), showWarnings = FALSE, recursive = TRUE)
    fwrite(dt_clean, clean_path)
    
    message(paste("  -> Saved QC'd data to:", clean_path))
    return(dt_clean)
}

if ("eid" %in% names(ms_cohort)) {
    run_qc("MS", ms_cohort)
}

if ("eid" %in% names(als_cohort)) {
    run_qc("ALS", als_cohort)
}
