#!/usr/bin/env Rscript
library(tidyverse)
source("helper.R")

args = commandArgs(trailingOnly=TRUE)
csv_dir = args[1]
out_prefix = args[2]

if (length(args) != 2) {
    exit("Usage: deg_dist.R csv_dir out_prefix")
}

######################################################################
## loading the table and creating some additional columns

tbl <- read.csv(paste(csv_dir, "degrees.csv", sep='/'))

## frequencies relative to number of vertices
tbl$frequency_rel <- tbl$frequency / tbl$n

## complementary cumulative distribution function of the frequencies
tbl$frequency_ccdf <- ave(tbl$frequency, tbl$graph,
                          FUN = function(x) rev(cumsum(rev(x))))

## ccdf of the relative frequencies
tbl$frequency_ccdf_rel <- tbl$frequency_ccdf / tbl$n

## maximum degree
tbl$max_deg <- ave(tbl$degree, tbl$graph, FUN = max)

## degrees relative to the maximum degree
tbl$degree_rel <- tbl$degree / tbl$max_deg

## restrict to subset of instances
#graphs <- levels(as.factor(tbl$graph))
#graphs
#tbl <- tbl[tbl$graph %in% 
#           c(graphs[37+72],
#             graphs[17+72],
#             graphs[7],
#             graphs[30+72],
#             graphs[32+72],
#             graphs[33+72]
#           ),]
## tbl <- tbl[tbl$graph %in% 
##            c(graphs[13],
##              graphs[37],
##              graphs[61]
##            ),]

######################################################################
## creating the plots

## relative frequencies depending on the degree as log-log plot
p1 <- ggplot(tbl, aes(x = degree, y = frequency_rel, color = graph)) +
    geom_point() +
    log_log_plot() +
    facet_wrap(~ graph) +
    xlab("degree d") +
    ylab("fraction of vertices with degree d") +
    theme(legend.position = "none")

## same as p1 but cumulative
p2 <- ggplot(tbl, aes(x = degree, y = frequency_ccdf_rel, color = graph)) +
    geom_step(direction = "vh") +
    log_log_plot() +
    facet_wrap(~ graph) +
    xlab("degree d") +
    ylab("fraction of vertices with degree > d") +
    theme(legend.position = "none")

## same as p1 but with normal axis (not log) and relative degrees
p3 <- ggplot(tbl, aes(x = degree_rel, y = frequency_rel, color = graph)) +
    geom_point() +
    facet_wrap(~ graph) +
    xlab("degree d / maximum degree") +
    ylab("fraction of vertices with degree d") +
    theme(legend.position = "none")

## same as p3 but cumulative
p4 <- ggplot(tbl, aes(x = degree_rel, y = frequency_ccdf_rel, color = graph)) +
    geom_step(direction = "vh") +
    facet_wrap(~ graph) +
    xlab("degree d / maximum degree") +
    ylab("fraction of vertices with degree > d") +
    theme(legend.position = "none")

######################################################################
## arranging the plots in a grid:
##   p1 p3
##   p2 p4
p <- ggarrange(p1 + rremove("x.title") + rremove("x.text"),
               p3 + rremove("xy.title") + rremove("x.text"),
               p2,
               p4 + rremove("y.title"),
               nrow = 2, ncol = 2)

######################################################################
## save to pdf
ggsave(paste("output.pdf/", out_prefix,"_degree_distribution.pdf", sep=''), plot = p, 
       width = 30, height = 30)
