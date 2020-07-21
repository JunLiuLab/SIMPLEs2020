import::from(here, here)
import::from(Rtsne, Rtsne)
library(ggplot2)
library(tidyverse)
library(ggpubr)


scale_level <- 10000
sccnt <- readRDS(here::here("10xGenomics",
                            "scRNAseq", "sc_10x_5cl_rgmm_cnt.rds"))
sc_genes <- rownames(sccnt)
cell_cluster <- gsub(".*:", "", colnames(sccnt))


get_tpm <- function(obs_counts_gbc, scale_level) {
  scale(obs_counts_gbc, center = F, scale = colSums(obs_counts_gbc)) *
    scale_level
}

log_tpm_rseq <- function(tpm_gbc, margin = 1) {
  tmp <- log(tpm_gbc + margin)
  ## add an small value to first posion if that row are equals
  ## to the same value.
  eps <- 0.1
  tmp_rowvar <- matrixStats::rowVars(tmp)
  modify_rows <- which(dplyr::near(tmp_rowvar, 0))
  tmp[modify_rows, 1] <- tmp[modify_rows, 1] + eps
  return(tmp)
}

tsneplot <- function(tsne_res, types, fnm) {
  dat <- data.frame(cbind(tsne_res), types)
  colnames(dat) <- c("TSNE1", "TSNE2", "Type")
  p <- ggplot(dat, aes(x = TSNE1, y = TSNE2, color = Type)) +
    geom_point() +
    theme_bw()
  pdf(file=fnm)
  plot(p)
  dev.off()
  return(p)
}

log_tpm <- log_tpm_rseq(get_tpm(sccnt, scale_level))

set.seed(0)
tsne_10x <- Rtsne(t(log_tpm),
  pca_scale = T, pca_center = T, initial_dims = 50,
  pca = T, check_duplicates = FALSE
)

p_10x <- tsneplot(tsne_10x$Y, cell_cluster, "tsne_10xGenomics.pdf")
