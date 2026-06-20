# VizMut-lolliplot — Installation script
# Run this script once before using VizMut-lolliplot
#
# First, install system dependencies (Ubuntu/Debian):
# sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev \
#   libgit2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
#   libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
#   libuv1-dev libcairo2-dev libxt-dev
#
# Then run: Rscript install_dependencies.R

message("=== VizMut-lolliplot — Dependency installer ===\n")

# limpiar locks huerfanos
lock_dir <- file.path(.libPaths()[1])
locks <- list.files(lock_dir, pattern="^00LOCK", full.names=TRUE)
if (length(locks) > 0) {
  message("Removing stale lock files...")
  unlink(locks, recursive=TRUE)
}

#---------------------------
# CRAN packages
#---------------------------
cran_packages <- c(
  "optparse", "dplyr", "stringr", "ggplot2", "ggrepel",
  "RColorBrewer", "ggnewscale", "httr", "jsonlite",
  "xml2", "R.utils", "scales", "patchwork", "testthat"
)

message("Checking CRAN packages...")
missing_cran <- cran_packages[!cran_packages %in%
                               installed.packages()[, "Package"]]

if (length(missing_cran) > 0) {
  message("Installing: ", paste(missing_cran, collapse=", "))
  install.packages(missing_cran,
                   repos   = "https://cloud.r-project.org",
                   Ncpus   = 2)
} else {
  message("All CRAN packages already installed.")
}

#---------------------------
# Bioconductor packages
#---------------------------
message("\nChecking Bioconductor packages...")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")

bioc_packages <- c("GenomicAlignments", "GenomicRanges", "rtracklayer")

missing_bioc <- bioc_packages[!bioc_packages %in%
                               installed.packages()[, "Package"]]

if (length(missing_bioc) > 0) {
  message("Installing: ", paste(missing_bioc, collapse=", "))
  for (pkg in missing_bioc) {
    message("  Installing ", pkg, "...")
    tryCatch(
      BiocManager::install(pkg, ask=FALSE, update=FALSE),
      error = function(e) message("  WARNING: failed to install ", pkg,
                                  " — ", e$message)
    )
  }
} else {
  message("All Bioconductor packages already installed.")
}

#---------------------------
# Verificacion final
#---------------------------
message("\n=== Verification ===")
all_packages <- c(cran_packages, "GenomicRanges", "rtracklayer")
for (pkg in all_packages) {
  status <- ifelse(pkg %in% installed.packages()[, "Package"], "OK", "MISSING")
  message(sprintf("  %-25s %s", pkg, status))
}

missing_final <- all_packages[!all_packages %in%
                               installed.packages()[, "Package"]]
if (length(missing_final) == 0) {
  message("\nAll dependencies installed successfully!")
  message("You can now run: Rscript main.R --help")
} else {
  message("\nWARNING: some packages failed to install:")
  for (p in missing_final) message("  - ", p)
  message("\nPlease install system dependencies first:")
  message("  sudo apt-get install -y libcurl4-openssl-dev libssl-dev \\")
  message("    libxml2-dev libgit2-dev libfontconfig1-dev libharfbuzz-dev \\")
  message("    libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev \\")
  message("    libjpeg-dev libuv1-dev libcairo2-dev libxt-dev")
  message("Then run: Rscript install_dependencies.R again")
}
