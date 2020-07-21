import::from(here, here)
import::from(Rtsne, Rtsne)
library(ggplot2)
library(tidyverse)
library(ggpubr)

## * for SymSim
p <- list(
  npop = 5, ncell = 300, ngene = 1000, scale_level_umi = 10000,
  scale_level_nonumi = 1e+06
)

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
  pdf(file = fnm)
  plot(p)
  dev.off()
  return(p)
}

tsne_symsim_groundtruth <- function(rep) {
  load(file = here::here(
    "simutool", "jobs", "symsim_data",
    stringr::str_c("sim", p$ncell, p$ngene, rep,
      ".RData",
      sep = "_"
    )
  ))
  logtpm_true <- log_tpm_rseq(get_tpm(
    true_counts$counts, p$scale_level_nonumi))
  set.seed(0)
  tsne_true <- Rtsne(t(logtpm_true),
    pca_scale = T, pca_center = T, initial_dims = 50,
    pca = T, check_duplicates = FALSE
  )
  p_true <- tsneplot(
    tsne_true$Y, as.character(cc),
    here::here(
      "deep_dive_simple", "tsne_symsim",
      stringr::str_glue("tsne_symsim_{rep}.pdf")
    )
  )
}

p_list <- lapply(1:20, tsne_symsim_groundtruth)
ggarrange(plotlist = p_list, nrow = 4, ncol = 5)
