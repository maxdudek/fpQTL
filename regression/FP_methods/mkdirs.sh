#!/bin/bash

for FP_METHOD in PRINT_no_gaussian; do
    mkdir -p $FP_METHOD/data
    mkdir -p $FP_METHOD/figures/missing_data
    mkdir -p $FP_METHOD/figures/regression_results_no_covariates
    mkdir -p $FP_METHOD/figures/regression_results_with_covariates
    mkdir -p $FP_METHOD/regression_results
done
