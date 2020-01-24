context("testing ASCAT")
library(testthat)

test_that("ASCAT is loading files as expected",{
            expect_output(loadASCAT("../testdata/results"),"ASCAT files:\n") 
            expect_output(loadASCAT("../testdata/results"),"../testdata/results/VariantCalling/ICGC_T_vs_ICGC_N/ASCAT/ICGC_T.LogR\n")
            expect_output(loadASCAT("../testdata/results"),"../testdata/results/VariantCalling/ICGC_T_vs_ICGC_N/ASCAT/ICGC_N.LogR\n")
            expect_output(loadASCAT("../testdata/results"),"../testdata/results/VariantCalling/ICGC_T_vs_ICGC_N/ASCAT/ICGC_T.BAF\n")
            expect_output(loadASCAT("../testdata/results"),"../testdata/results/VariantCalling/ICGC_T_vs_ICGC_N/ASCAT/ICGC_N.BAF\n")
            expect_output(loadASCAT("../testdata/results"),"../testdata/results/VariantCalling/ICGC_T_vs_ICGC_N/ASCAT/ICGC_T.cnvs.txt")
})

test_that("calculateASCAT",{
            expect_output(calculateASCAT(),"../testdata/results/VariantCalling/ICGC_T_vs_ICGC_N/ASCAT/ICGC_T.LogR\n")
})
