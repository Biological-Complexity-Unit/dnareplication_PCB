#!/usr/bin/env Rscript
library(data.table)
library(parallel)
source("abundance.R")

GENOMES <- list(
  EColiMG1655=list(name="E. Coli MG1655", L=4641652, O=3925975)
)

MODELS <- list(
  vt=vt.harmonic.model,
  vd=vd.harmonic.model,
  vx=vx.harmonic.model
)

# Decode arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("fit_oscillation.R <genome> <model> <sample> <lambda> <reference> <output> <cores>\n", file=stderr())
  exit(1)
}
argi <- 1
if (!(args[[argi]] %in% names(GENOMES)))
  stop("Unknown genome ", args[[argi]])
genome <- GENOMES[[args[[argi]]]]
argi <- argi + 1
if (!(args[[argi]] %in% names(MODELS)))
  stop("Unknown model ", args[[argi]])
model <- MODELS[[args[[argi]]]]
argi <- argi + 1
sample.file <- args[[argi]]
if (!file.exists(sample.file))
  stop("Sample file ", sample.file, " does not exist")
argi <- argi + 1
sample.lambda <- as.numeric(args[[argi]])
argi <- argi + 1
reference.file <- args[[argi]]
if (!file.exists(reference.file))
  stop("Reference file ", reference.file, " does not exist")
argi <- argi + 1
output.file <- args[[argi]]
if (file.exists(output.file))
  stop("Output file ", output.file, " already exist")
argi <- argi + 1
if (length(args) >= argi) {
  cores <- as.integer(args[[argi]])
} else {
  cores <- detectCores()
}
cat(paste0("Genome: ", genome$name, "\n"))
cat(paste0("Model: ", model$name, "\n"))
cat(paste0("Sample: ", sample.file, "\n"))
cat(paste0("Sample growth rate (1/h): ", sample.lambda, "\n"))
cat(paste0("Reference: ", reference.file, "\n"))
cat(paste0("Cores: ", cores, "\n"))

# Read inputs
sample <- read.counts.deepak(sample.file, genome$L, genome$O)[order(x)]
reference <- read.counts.deepak(reference.file, genome$L, genome$O)[order(x)]

# Remove bins missing in sample or reference
data <- sample[reference, list(x, sample=x.count, reference=i.count), on=.(x), nomatch=0]
cat(paste0("Overlapping bins: ", nrow(data), ", (",
           round(100*nrow(data)/nrow(sample)), "% of sample, ",
           round(100*nrow(data)/nrow(reference)), "% of reference)\n"))

# Compute normalized abundances
d <- normalize.counts(data$x, data$sample, data$reference)

# Fit
ps <- fit.abundance.stepwise(d$x, d$y, d$w, L=genome$L, model=model,
                             lambda=sample.lambda/3600, mc.cores=cores)

tab <- rbindlist(lapply(names(ps), function(n) {
  p <- ps[[n]]
  res <- attr(p, "residual")
  logl <- d$logl.const - res 
  as.data.table(c(list(type=n), p[c("lambda", "D", "v0", "delta", "omega", "phi")], list(residual=res, logl=logl)))
}))

cat("Results:\n")
print(tab)

# Write output
write.table(file=output.file, tab, sep="\t", col.names=TRUE)
