#!/usr/bin/env Rscript
library(tidyverse)
source("helper.R")
source("deg_variance.R")


args = commandArgs(trailingOnly=TRUE)
csv_dir = args[1]
out_prefix = args[2]

if (length(args) != 2) {
    exit("Usage: plotting.R csv_dir out_prefix")
}

cat("Reading csvs from", csv_dir, "out_prefix", out_prefix)

summary_tbl <- read.csv(paste(csv_dir, 'summary.csv', sep='/'))

deg_tbl <- read.csv(paste(csv_dir, "degrees.csv", sep='/'))
degvar_norm_tbl <- deg_var_norm(deg_tbl)
degvar_norm_tbl <- degvar_norm_tbl[, c("graph", "deg_var_norm")]

summary_tbl <- merge(x=summary_tbl, y=degvar_norm_tbl, by="graph")

print(colnames(summary_tbl))

p1 <-  ggplot(summary_tbl, aes(x=n, y=m)) +
    geom_point() +
    xlab("n") +
    ylab("m")

p2 <- ggplot(summary_tbl, aes(summary_tbl$avg_deg)) +
    geom_histogram() +
    xlab("average degree d") +
    ylab("number of graphs with average degree d")

p4 <- ggplot(summary_tbl, aes(summary_tbl$deg_var_norm)) +
    geom_histogram() +
    xlab("normalized degree variance") +
    ylab("number of graphs")
p5 <- ggplot(summary_tbl, aes(summary_tbl$clustering)) +
    geom_histogram() +
    xlab("clustering coefficient c") +
    ylab("number of graphs with clustering coefficient c")
p6 <- ggplot(summary_tbl, aes(x=avg_deg, y=deg_var_norm, size=n)) +
    geom_point() +
    xlab("Average degree") +
    ylab("Normalized degree variance")
p7 <- ggplot(summary_tbl, aes(x=avg_deg, y=clustering, size=n)) +
    geom_point() +
    xlab("Average degree") +
    ylab("Clustering coefficient")
p8 <- ggplot(summary_tbl, aes(x=deg_var_norm, y=clustering, size=n)) +
    geom_point() +
    xlab("Normalized degree variance") +
    ylab("Clustering coefficient")
p9 <- ggplot(summary_tbl, aes(x=deg_var_norm, y=clustering, size=avg_deg)) +
    geom_point() +
    xlab("Normalized degree variance") +
    ylab("Clustering coefficient")

p <- ggarrange(p1,p2,p4,p5,p6,p7,p8,p9, nrow=9, ncol=1)

ggsave(paste("output.pdf/", out_prefix, "_summary.pdf", sep=""), plot = p,
       width = 10, height = 50, limitsize = FALSE)
