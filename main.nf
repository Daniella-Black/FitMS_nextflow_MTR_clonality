#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey, all, clonal_any, clonal_NA, clonal_early, clonal_late, subclonal -> [tumour_sample_platekey, file(all), file(clonal_any), file(clonal_NA), file(clonal_early), file(clonal_late), file(subclonal)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    container = 'dockeraccountdani/fitms2:latest' 
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}/$tumour_sample_platekey", mode: 'copy'
    
    input:
    set val(tumour_sample_platekey), file(all), file(clonal_any), file(clonal_NA), file(clonal_early), file(clonal_late), file(subclonal) from ch_input

    output:
    //file "*_SNV_catalogues.pdf"
    //file "*_clonality_state_catalogue.csv"
    //file "exposures.tsv"
    //path "results/*"
    file "clonality.txt"
    file "clonality_print.txt"
    file "clonality_in.txt"
    
    script:
    """
    fitms_nf.R '$tumour_sample_platekey' '$all' '$clonal_any' '$clonal_NA' '$clonal_early' '$clonal_late' '$subclonal'
    """ 
}
