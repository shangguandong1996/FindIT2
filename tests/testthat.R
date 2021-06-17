library(testthat)
library(FindIT2)
if (!requireNamespace("TxDb.Athaliana.BioMart.plantsmart28")) {
    stop("unable to load TxDb.Athaliana.BioMart.plantsmart28")
}

test_check("FindIT2")
