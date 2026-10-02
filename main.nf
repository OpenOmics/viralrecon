#!/usr/bin/env nextflow
// Stock viralrecon remains available. The TCR branch accepts POD5 inputs only.
if (params.tcr_pod5) {
    include { TCR_POD5 } from './workflows/tcr_pod5'
} else {
    include { VIRALRECON_STOCK_ENTRY } from './main_viralrecon_stock'
}
workflow {
    if (params.tcr_pod5) {
        TCR_POD5()
    } else {
        VIRALRECON_STOCK_ENTRY()
    }
}
