#!/bin/bash

CORES=32

run() {
	name=$1; shift
	echo "$name: " $*
	exec $*
}

run 17-1-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/17-1-10kb.dat.bz2 0.183853 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/17-1-17.dat $CORES
run 17-1-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/17-1-10kb.dat.bz2 0.183853 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/17-1-27.dat $CORES
run 17-1-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/17-1-10kb.dat.bz2 0.183853 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/17-1-37.dat $CORES
run 17-2-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/17-2-10kb.dat.bz2 0.207921 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/17-2-17.dat $CORES
run 17-2-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/17-2-10kb.dat.bz2 0.207921 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/17-2-27.dat $CORES
run 17-2-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/17-2-10kb.dat.bz2 0.207921 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/17-2-37.dat $CORES
run 17-3-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/17-3-10kb.dat.bz2 0.214417 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/17-3-17.dat $CORES
run 17-3-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/17-3-10kb.dat.bz2 0.214417 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/17-3-27.dat $CORES
run 17-3-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/17-3-10kb.dat.bz2 0.214417 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/17-3-37.dat $CORES
run 22-1-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/22-1-10kb.dat.bz2 0.446480 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/22-1-17.dat $CORES
run 22-1-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/22-1-10kb.dat.bz2 0.446480 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/22-1-27.dat $CORES
run 22-1-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/22-1-10kb.dat.bz2 0.446480 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/22-1-37.dat $CORES
run 22-2-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/22-2-10kb.dat.bz2 0.501008 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/22-2-17.dat $CORES
run 22-2-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/22-2-10kb.dat.bz2 0.501008 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/22-2-27.dat $CORES
run 22-2-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/22-2-10kb.dat.bz2 0.501008 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/22-2-37.dat $CORES
run 22-3-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/22-3-10kb.dat.bz2 0.475542 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/22-3-17.dat $CORES
run 22-3-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/22-3-10kb.dat.bz2 0.475542 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/22-3-27.dat $CORES
run 22-3-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/22-3-10kb.dat.bz2 0.475542 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/22-3-37.dat $CORES
run 27-1-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/27-1-10kb.dat.bz2 0.894485 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/27-1-17.dat $CORES
run 27-1-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/27-1-10kb.dat.bz2 0.894485 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/27-1-27.dat $CORES
run 27-1-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/27-1-10kb.dat.bz2 0.894485 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/27-1-37.dat $CORES
run 27-2-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/27-2-10kb.dat.bz2 0.898088 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/27-2-17.dat $CORES
run 27-2-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/27-2-10kb.dat.bz2 0.898088 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/27-2-27.dat $CORES
run 27-2-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/27-2-10kb.dat.bz2 0.898088 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/27-2-37.dat $CORES
run 27-3-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/27-3-10kb.dat.bz2 0.965962 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/27-3-17.dat $CORES
run 27-3-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/27-3-10kb.dat.bz2 0.965962 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/27-3-27.dat $CORES
run 27-3-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/27-3-10kb.dat.bz2 0.965962 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/27-3-37.dat $CORES
run 32-1-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/32-1-10kb.dat.bz2 1.698130 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/32-1-17.dat $CORES
run 32-1-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/32-1-10kb.dat.bz2 1.698130 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/32-1-27.dat $CORES
run 32-1-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/32-1-10kb.dat.bz2 1.698130 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/32-1-37.dat $CORES
run 32-2-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/32-2-10kb.dat.bz2 1.645930 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/32-2-17.dat $CORES
run 32-2-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/32-2-10kb.dat.bz2 1.645930 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/32-2-27.dat $CORES
run 32-2-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/32-2-10kb.dat.bz2 1.645930 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/32-2-37.dat $CORES
run 32-3-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/32-3-10kb.dat.bz2 1.459090 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/32-3-17.dat $CORES
run 32-3-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/32-3-10kb.dat.bz2 1.459090 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/32-3-27.dat $CORES
run 32-3-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/32-3-10kb.dat.bz2 1.459090 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/32-3-37.dat $CORES
run 37-1-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/37-1-10kb.dat.bz2 2.021380 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/37-1-17.dat $CORES
run 37-1-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/37-1-10kb.dat.bz2 2.021380 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/37-1-27.dat $CORES
run 37-1-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/37-1-10kb.dat.bz2 2.021380 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/37-1-37.dat $CORES
run 37-2-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/37-2-10kb.dat.bz2 1.976010 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/37-2-17.dat $CORES
run 37-2-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/37-2-10kb.dat.bz2 1.976010 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/37-2-27.dat $CORES
run 37-2-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/37-2-10kb.dat.bz2 1.976010 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/37-2-37.dat $CORES
run 37-3-17 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/37-3-10kb.dat.bz2 1.937030 data/bhat/raw/17-stationary-10kb.dat.bz2 data/bhat/est.vd/37-3-17.dat $CORES
run 37-3-27 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/37-3-10kb.dat.bz2 1.937030 data/bhat/raw/27-stationary-10kb.dat.bz2 data/bhat/est.vd/37-3-27.dat $CORES
run 37-3-37 ./fit_oscillation.R EColiMG1655 vd data/bhat/raw/37-3-10kb.dat.bz2 1.937030 data/bhat/raw/37-stationary-10kb.dat.bz2 data/bhat/est.vd/37-3-37.dat $CORES
