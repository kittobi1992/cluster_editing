## Variance of the degree distribution.  The input must be a histogram
## with columns 'graph', 'degree', and 'frequency'.  The output has
## columns 'graph', 'n', 'm', 'avg_deg', and 'deg_var'.
deg_var <- function (tbl) {
    tbl$n <- ave(tbl$frequency, tbl$graph, FUN = sum)
    tbl$m <- ave(tbl$frequency * tbl$degree / 2, tbl$graph, FUN = sum)
    tbl$avg_deg <- tbl$m * 2 / tbl$n
    res <- aggregate(
        list(deg_var = (tbl$degree - tbl$avg_deg)^2 / tbl$n),
        list(graph = tbl$graph, n = tbl$n, m = tbl$m, avg_deg = tbl$avg_deg),
        FUN = sum)
    return(res)
}

## Normalized degree variance (https://arxiv.org/abs/1803.03057). The
## input table must be a histogram with columns 'graph', 'degree', and
## 'frequency'.  The output has columns 'graph', 'n', 'm', and
## 'deg_var_norm'.
deg_var_norm <- function (tbl) {
    res <- aggregate(
        list(n = tbl$frequency,
             m = tbl$frequency * tbl$degree / 2,
             sum_of_sq_deg = tbl$degree^2 * tbl$frequency),
        list(graph = tbl$graph), FUN = sum)

    ## first summand of v(G) (see paper)
    res$v_first_summand <- res$sum_of_sq_deg / (res$n - 1)

    ## second summand of v(G) (see paper)
    res$v_secnd_summand <- res$m / res$n * res$m / (res$n - 1) * 4

    ## density d (see paper)
    res$d <- 2 * res$m / (res$n * (res$n - 1))

    ## normalization factor
    res$factor <- (res$n - 1) / ((1 - res$d) * res$n * res$m)

    ## normalized variance
    res$deg_var_norm <- res$factor * (res$v_first_summand - res$v_secnd_summand)
    
    return(res[, c("graph", "n", "m", "deg_var_norm")])
}

