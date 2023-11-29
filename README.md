# Description

Implementation of the bacterial and eukaryotic model of Pflug et al, Genome replication in
asynchronously growing microbial populations, *PLOS Computational Biology*.

# Contents

| File | Description |
| ---- | ----------- |
| ./bacteria/parmlapply.R | Parallelization helper |
| ./bacteria/abundance.R | Computes abundance profiles under different bacterial models |
| ./bacteria/simulate.R | Simulates abundance profiles under different bacterial models |
| ./bacteria/simulate_oscillation.R | Simulates an abundance profile using random parameters and recovers them by fitting |
| ./bacteria/fit_oscillation.R | Fit the bacterial model |
| ./bacteria/fit_oscillation.vt.sh | Fit the bacterial v(t) model to the data of Bhat et al., 2022 |
| ./bacteria/fit_oscillation.vd.sh | Fit the bacterial v(d) model to the data of Bhat et al., 2022 |
| ./bacteria/data/bhat2022/est.vt.aic.dat.bz2 | Table of fitted parameters for v(t) |
| ./bacteria/data/bhat2022/est.vt | Per-sample fitting results for v(t) |
| ./bacteria/data/bhat2022/est.vd.aic.dat.bz2 | Table of fitted parameters for v(d) |
| ./bacteria/data/bhat2022/est.vd | Per-sample fitting results for v(d) |
| ./bacteria/data/bhat2022/raw | Raw abundance data of Bhat et al., 2023 |
| ./yeast/abundance.R | Computes abundance profiles under the Eukaryotic model |
| ./yeast/SimulateData.Rmd | Simulates an abundance profile under the Eukaryotic model |
| ./yeast/data/s_cerevisiae.W303.fna.gz | The reference genome for S. Cerevisiae W303 |
| ./yeast/data/s_cerevisiae.W303.gff.gz | Genome annotations for S. Cerevisiae W303 lifted from S228C | 
| ./yeast/data/simulated/CV-0.04/ | Inference using simulated data |
| ./yeast/data/simulated/CV-0.04/Truth/chr4.dat.bz2 | Ground truth for Chr. IV |
| ./yeast/data/simulated/CV-0.04/Inferred/cluster.sh | Inference script |
| ./yeast/data/simulated/CV-0.04/Inferred/origins.dat.bz2 | Inferred origins for Chr. IV |
| ./yeast/data/simulated/CV-0.04/Inferred/data.dat.bz2 | Inferrence results |
| ./yeast/data/simulated/CV-0.04/Inferred/code.py | Inference algorithm |
| ./yeast/data/simulated/CV-0.04/Inferred/abundance.dat.bz2 | Fitted abundance profile for Chr. IV|
| ./yeast/data/simulated/CV-0.04/Inferred/parameters.dat.bz2 | Inferred parameters |
| ./yeast/data/simulated/CV-0.04/Inferred/anneal.log.bz2 | Inference log file |
| ./yeast/data/mueller2014/Experiments | Raw data of Mueller et al., 2014 |
| ./yeast/data/mueller2014/Inferred/all.varest | Inference results for the data of Mueller et al., 2014 |
| ./yeast/data/mueller2014/Inferred/all.varest/cluster.sh | Inference script |
| ./yeast/data/mueller2014/Inferred/all.varest/origins.dat.bz2 | Inferred origins for Chr. IV |
| ./yeast/data/mueller2014/Inferred/all.varest/data.dat.bz2 | Inferrence results |
| ./yeast/data/mueller2014/Inferred/all.varest/code.py | Inference algorithm |
| ./yeast/data/mueller2014/Inferred/all.varest/abundance.dat.bz2 | Fitted abundance profile |
| ./yeast/data/mueller2014/Inferred/all.varest/cluster.reps.sh | Inference scripts for replicates |
| ./yeast/data/mueller2014/Inferred/all.varest/parameters.dat.bz2 | Inferred parameters |
| ./yeast/data/mueller2014/Inferred/all.varest/anneal.log.bz2 | Inference log file |
| ./yeast/data/mueller2014/Inferred/all.varest/reps | Replicated inference runs |
