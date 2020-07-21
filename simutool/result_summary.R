library(here)
library(matrixStats)
library(tidyverse)

p <- list(
  nrep = 20, npop = 5, ncell = 300, ngene = 1000,
  pcs = seq(2, 10, 2)
)

load_model_result <- function(modelnm, k = 2, m = 1,
                              dirnm = here::here("simutool", "jobs", "result"),
                              outmodelnm = toupper(modelnm), ...) {
  if (modelnm != "simple") {
    modelnm <- modelnm
    fnm <- paste0(dirnm, "/", str_c(modelnm, 2, 1, ".rds", sep = "_"))
  } else {
    modelnm <- str_c(modelnm, k, m, sep = "_")
    fnm <- paste0(dirnm, "/", modelnm, "_.rds")
  }
  npc <- length(p$pcs)
  r <- readRDS(fnm)
  ## clari short for cluster ari
  ## explm short for experiment platform: umi or nonumi
  ## aspect for deauc or clari
  get_auc_sums <- function(aucs, explm) {
    data.frame(
      mean = mean(aucs), std = sd(rowMeans(aucs)) / sqrt(p$nrep),
      explm = explm, aspect = "deauc"
    )
  }
  ## deprecated
  get_ari_sums_matrix <- function(aris, explm) {
    data.frame(
      mean = apply(aris, 2, mean), std = apply(aris, 2, sd) / sqrt(p$nrep),
      explm = rep(explm, npc), aspect = rep("clari", npc)
    )
  }

  get_ari_sums_array <- function(aris, explm, col = npc) {
    data.frame(
      mean = mean(aris[, col]), std = sd(aris[, col]) / sqrt(p$nrep),
      explm = explm, aspect = "clari"
    )
  }

  de_auc_umi <- get_auc_sums(r$de_auc_umi, "umi")
  de_auc_nonumi <- get_auc_sums(r$de_auc_nonumi, "nonumi")
  clust_ari_umi <- get_ari_sums_array(r$clust_ari_umi, explm = "umi", col = npc)
  clust_ari_nonumi <- get_ari_sums_array(r$clust_ari_nonumi,
    explm = "nonumi", col = npc
  )

  mysum <- rbind.data.frame(de_auc_umi, de_auc_nonumi,
    clust_ari_umi, clust_ari_nonumi,
    make.row.names = FALSE, stringsAsFactors = FALSE
  )
  mysum[["model"]] <- c(outmodelnm)
  mysum[["pcs"]] <- c(0, 0, p$pcs[npc], p$pcs[npc])
  return(mysum)
}


draw_bar_plot <- function(e = "nonumi", a = "deauc",
                          data = all_res, lpos = "none",
                          mytitle = "AUC", myangle = 0) {
  p <- data %>%
    filter(explm == e, aspect == a) %>%
    ggplot(data = ., aes(x = model, y = mean, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = mean - std, ymax = mean + std),
      width = .2, position = position_dodge(.9)
    ) +
    theme(
      axis.text.x = element_text(
        angle = myangle, color = "black",
        size = 25, family = "Helvetica"
      ),
      axis.text.y = element_text(color = "black",
                                 family = "Helvetica", size = 13),
      axis.title.x = element_blank(),
      axis.title.y = element_text(color = "black",
                                  family = "Helvetica", size = 25),
      title = element_text(
        family = "Helvetica", color = "black", size = 24,
        face = "bold"
      ),
      plot.title = element_text(hjust = 0.5),
      legend.position = lpos
    ) +
    ylab(mytitle)
  return(p)
}

## * load_all result
control_res <- load_model_result("control", outmodelnm = "Control")
rmagic_res <- load_model_result("rmagic", outmodelnm = "MAGIC")
saver_res <- load_model_result("saver", outmodelnm = "SAVER")
scimpute_res <- load_model_result("scimpute", outmodelnm = "scImpute")
scrabble_res <- load_model_result("scrabble", outmodelnm = "SCRABBLE")
viper_res <- load_model_result("viper", outmodelnm = "VIPER")
scvi_res <- load_model_result("scvi", outmodelnm = "scVI")

simple21_res <- load_model_result("simple", 2, 1)
simple41_res <- load_model_result("simple", 4, 1)
simple81_res <- load_model_result("simple", 8, 1)
simple101_res <- load_model_result("simple", 10, 1)
simple22_res <- load_model_result("simple", 2, 2)
simple42_res <- load_model_result("simple", 4, 2)
simple82_res <- load_model_result("simple", 8, 2)
simple102_res <- load_model_result("simple", 10, 2)
simple24_res <- load_model_result("simple", 2, 4)
simple44_res <- load_model_result("simple", 4, 4)
simple84_res <- load_model_result("simple", 8, 4)
simple104_res <- load_model_result("simple", 10, 4)
simple25_res <- load_model_result("simple", 2, 5)
simple45_res <- load_model_result("simple", 4, 5)
simple85_res <- load_model_result("simple", 8, 5)
simple105_res <- load_model_result("simple", 10, 5)


simple_res <- rbind(
  simple21_res,
  simple41_res,
  simple81_res,
  simple101_res,
  simple22_res,
  simple42_res,
  simple82_res,
  simple102_res,
  simple24_res,
  simple44_res,
  simple84_res,
  simple104_res,
  simple25_res,
  simple45_res,
  simple85_res,
  simple105_res
)

all_res <- rbind(
  control_res,
  rmagic_res,
  saver_res,
  scimpute_res,
  scrabble_res,
  viper_res,
  scvi_res,
  simple101_res,
  simple41_res,
  simple25_res
)

## * Figure
deauc_nonumi_p <- draw_bar_plot("nonumi", "deauc", myangle = 90)
ggsave("deauc_nonumi_all.pdf",
  device = "pdf",
  plot = deauc_nonumi_p, width = 11, height = 11
  )

deauc_umi_p <- draw_bar_plot("umi", "deauc", myangle = 90)
ggsave("deauc_umi_all.pdf",
  device = "pdf",
  plot = deauc_umi_p, width = 11, height = 11
  )

ari_nonumi_p <- draw_bar_plot("nonumi", "clari",
  mytitle = "aRI", myangle = 90
)
ggsave("ari_nonumi_all.pdf",
  device = "pdf",
  plot = ari_nonumi_p, width = 11, height = 11
  )

ari_umi_p <- draw_bar_plot("umi", "clari",
  mytitle = "aRI", myangle = 90
)
ggsave("ari_umi_all.pdf",
  device = "pdf",
  plot = ari_umi_p, width = 11, height = 11
)

deauc_nonumi_p <- draw_bar_plot("nonumi", "deauc",
  data = simple_res, "none", myangle = 90
)
ggsave("deauc_nonumi_simples.pdf",
  device = "pdf", plot = deauc_nonumi_p, width = 11, height = 11
)

deauc_umi_p <- draw_bar_plot("umi", "deauc",
  data = simple_res, "none", myangle = 90
)
ggsave("deauc_umi_simples.pdf",
  device = "pdf", plot = deauc_umi_p, width = 11, height = 11
)

ari_nonumi_p <- draw_bar_plot("nonumi", "clari",
  data = simple_res, lpos = "none", mytitle = "aRI", myangle = 90
)
ggsave("ari_nonumi_simples.pdf",
  device = "pdf", plot = ari_nonumi_p, width = 11, height = 11
)


ari_umi_p <- draw_bar_plot("umi", "clari",
  data = simple_res, lpos = "none", mytitle = "aRI", myangle = 90
)
ggsave("ari_umi_simples.pdf",
  device = "pdf", plot = ari_umi_p, width = 11, height = 11
)
