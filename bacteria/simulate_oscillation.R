#!/usr/bin/env Rscript
library(data.table)
library(parallel)
library(argparser)
library(ggplot2, quietly=TRUE)
source("abundance.R")
source("simulate.R")

GENOMES <- list(
  EColiMG1655=list(name="E. Coli MG1655", L=4641652, O=3925975)
)

MODELS <- list(
  vt=vt.harmonic.model,
  vd=vd.harmonic.model,
  vx=vx.harmonic.model
)

# Decode arguments
p <- arg_parser("simulate_oscillations")
p <- add_argument(p, "genome", help="genome to use")
p <- add_argument(p, "model", help="model to use")
p <- add_argument(p, "output", help="output file")
p <- add_argument(p, "--coverage", help="coverage", type="numeric", default=1000)
p <- add_argument(p, "--window", help="window size", type="numeric", default=1e4)
p <- add_argument(p, "--ngenomes", help="window size", type="numeric", default=1e5)
p <- add_argument(p, "--plot", help="output abundance plot")
p <- add_argument(p, "--cores", help="number of cores to use", type="numeric", default=1)
args <- parse_args(p)
if (!(args$genome %in% names(GENOMES)))
  stop("Unknown genome ", args$genome)
genome <- GENOMES[[args$genome]]
if (!(args$model %in% names(MODELS)))
  stop("Unknown model ", args$model)
model <- MODELS[[args$model]]
if (file.exists(args$output))
  stop("Output file ", args$output, " already exist")
cat(paste0("Genome: ", genome$name, "\n"))
cat(paste0("Model: ", model$name, "\n"))
cat(paste0("Output: ", args$output, "\n"))
cat(paste0("Plot: ", if (!is.null(args$plot)) args$plot else "-", "\n"))
cat(paste0("Coverage: ", args$coverage, "\n"))
cat(paste0("Window: ", args$window, "\n"))
cat(paste0("Num. genomes to simulate: ", args$ngenomes, "\n"))
cat(paste0("Cores: ", args$cores, "\n"))

message("*** Generating random parameter vector")
p.true <- list(
  L=genome$L,
  lambda=1/3600,
  model=model
)
p.true$v0 <- 10^runif(1, min=2, max=log10(2000))
p.true$D <- 10^runif(1, min=3, max=8)
p.true$delta <- runif(1, min=0, max=0.75)
p.true$omega <- if (args$model == "vt") {
  2^runif(1, min=-2, max=2) *2*pi * p.true$v0 / genome$L
} else {
  2^runif(1, min=-2, max=2) *2*pi / genome$L
}
p.true$phi <- runif(1, min=0, max=2*pi)
pn <- c("L", "lambda", "v0", "D", "delta", "omega", "phi")
cat(paste0("  ", pn, ": ", lapply(p.true[pn], signif, 3), collapse="\n"))
cat("\n")

message("*** Simulating ", args$ngenomes, " genomes")
cluster <- if (exists("CLUSTER")) {
  CLUSTER
} else if (args$cores > 1) {
  makeCluster(args$cores)
} else NULL
if (!is.null(cluster))
  clusterExport(cluster, ls())
sim <- simulate.abundance(p.true, dx=args$window, N=args$ngenomes, cluster=cluster)
if (!is.null(cluster))
  stopCluster(cluster)

message("*** Computing sample and reference counts")
ref <- rpois(length(sim$x), lambda=args$coverage)
sample.true <- args$coverage * length(sim$a) * sim$a / sum(sim$a)
sample <- rpois(length(sample.true), lambda=sample.true)
if (!is.na(args$plot)) {
  message("*** Plotting simulated abundances to ", args$plot)
  plot <- ggplot(data.table(x=sim$x, y=sample, ytrue=sample.true)) +
    geom_line(aes(x=x, y=y)) +
    geom_line(aes(x=x, y=ytrue), linetype="dashed")
  ggsave(args$plot, plot)
}

message("*** Re-estimating parameters")
d <- normalize.counts(sim$x, sample, ref)
ps <- fit.abundance.stepwise(d$x, d$y, d$w, L=genome$L, lambda=p.true$lambda, model=model, mc.cores=args$cores)
tab <- rbindlist(lapply(c("p.true", names(ps)), function(n) {
  if (n == "p.true")  {
    p <- p.true
    m <- evaluate.model(p.true, d)
    res <- m$residual
    logl <- m$logl
  } else {
    p <- ps[[n]]
    res <- attr(p, "residual")
    logl <- d$logl.const - res 
  }
  as.data.table(c(list(type=n), p[c("lambda", "D", "v0", "delta", "omega", "phi")], list(residual=res, logl=logl)))
}))

cat("Results:\n")
print(tab)

# Write output
message("Saving results to ", args$output)
write.table(file=args$output, tab, sep="\t", col.names=TRUE)
