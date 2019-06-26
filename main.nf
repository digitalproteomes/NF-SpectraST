// Create pseudoMS2 with DIAUmpire
// Searches with Comet and Tandem
// Pooled PeptideProphet for Comet and Tandem separately
// InterProphet on pooled PeptideProphets
// ProteinProphet on InterProphet
// Mayu on InterProphet

// This branch filters proteins at 1%FDR and includes all
// corresponding peptides at 1%FDR

if(params.help) {

    log.info "PQP generation workflow"
    log.info "----------------------"
    log.info ""
    log.info "	  DIA files"
    log.info "	      +"
    log.info "	      |"
    log.info "	      v"
    log.info "	to_umpire.txt"
    log.info "	      +"
    log.info "	      |"
    log.info "	+-----v------+"
    log.info "	| DIA-Umpire | +---> pseudo DDA files"
    log.info "	+------------+       +"
    log.info "	                     |"
    log.info "	                     v"
    log.info "	  DDA files  +-----> |"
    log.info "	                     |"
    log.info "	                     |  +-----------------+  +----------------+"
    log.info "	                     +--> search engine 1 +--> PeptideProphet +--+"
    log.info "	                     |  +-----------------+  +----------------+  |"
    log.info "	                     |                                           |"
    log.info "	                     |          ....                          +--+"
    log.info "	                     |                                           |"
    log.info "	                     |  +-----------------+  +----------------+  |"
    log.info "	                     +--> search engine 1 +--> PeptideProphet +--+"
    log.info "	                        +-----------------+  +----------------+  |"
    log.info "	+----------------+                                               |"
    log.info "	| ProteinProphet <------+                                        |"
    log.info "	+----------------+      |                                        |"
    log.info "	                        |                                        |"
    log.info "	+-------------+         |                 +----------+           |"
    log.info "	| calctppstat <---------------------------+ iProphet <-----------+"
    log.info "	+-------------+         |                 +----+-----+"
    log.info "	                        |                      |"
    log.info "	+------+                |                      |"
    log.info "	| Mayu <----------------+                      |"
    log.info "	+--+---+                                       |"
    log.info "	   |                                           |"
    log.info "	   |  1% FDR Protein  +-----------+            |"
    log.info "	   +------------------> SpectraST <------------+"
    log.info "	      Level           +-----+-----+"
    log.info "	                            |"
    log.info "	                            v"
    log.info "	                Consensus Spectral Library"
    log.info "	                            +"
    log.info "	                            |"
    log.info "	               +------------v------------+"
    log.info "	               | OpenSwathAssayGenerator |"
    log.info "	               +------------+------------+"
    log.info "	                            |"
    log.info "	               +------------v------------+"
    log.info "	               | OpenSwathDecoyGenerator |"
    log.info "	               +------------+------------+"
    log.info "	                            |"
    log.info "	                            v"
    log.info "	                PQP assay library with decoys"
    log.info ""
    log.info ""
    log.info "Options:"
    log.info "  --help:              show this message and exit"
    log.info "  --dia_folder:        file to be DIAUmpired (default: $params.dia_folder)"
    log.info "  --diau_se_params:    parameter file for DIAUmpire (default: $params.diau_se_params)"
    log.info "  --dda_folder:        files to be searched (default: $params.dda_folder)"
    log.info "  --comet_params:      comet parameter file (default: $params.comet_params)"
    log.info "  --protein_db:        comet parameter file (default: $params.protein_db)"
    log.info "  --tpp:               TPP options (default: $params.tpp)"
    log.info "  --decoy:             decoy prefix (default: $params.decoy)"
    log.info "  --rt_file:           RT normalization file (default: $params.rt_file"
    log.info "  --st_fragmentation:  fragmentation type to use in SpectraST (default: $params.st_fragmentation)"
    log.info "  --st_fix_mods:       sed script file for conversion of TPP modifications to UniMod modifications (default: $params.st_fix_mods)"
    log.info "  --swath_window_file: file with the definition of the SWATH windows (default: $params.swath_window_file)"
    log.info ""
    log.info "The final library (SpecLib_opt_dec.pqp) will be in Results/SpectraST"
    log.info ""
    exit 1
}


process diaUmpire {
    scratch 'ram-disk'
    stageInMode "copy"
    tag "$dia_file"

    cpus params.diau_threads
        
    // Generate pseudoDDA files for each DIA file
    input:
    file dia_file from file("${params.dia_folder}/*.mzXML")
    file diau_se_params from file(params.diau_se_params)

    output:
    file '*.mgf' into diaUmpireOut

    script:
    """
    sed -i 's,Thread = 0,Thread = $params.diau_threads,' $diau_se_params
    java -jar -Xmx32G /home/biodocker/bin/DIA-Umpire/v2.1.2/DIA_Umpire_SE.jar $dia_file $diau_se_params
    """

}
    // script:
    // """
    // sed -i 's,Thread = 0,Thread = $params.diau_threads,' $diau_se_params
    // dia_umpire_se -XX:MaxRAMPercentage=80 -XX:+UseContainerSupport $dia_file $diau_se_params
    // """


process mgf2mzxml {
    // Convert mgfs generated by DIAUmpire to mzXML
    tag "$mgf_file"

    publishDir 'Results/DIAUmpire', mode: 'link'

    input:
    // For each DIA file DIA-Umpire will generate 3 pseudo mgf files
    file mgf_file from diaUmpireOut.flatten()

    output:
    file '*.mzXML' into mgf2mzxmlOut
    
    """
    msconvert $mgf_file --mzXML
    """
}


mgf2mzxmlOut.into{mgf2mzxmlOut1; mgf2mzxmlOut2}


process msfraggerIndex {
    tag "$mzXML_fragger"
    
    cpus params.fragger_threads
    
    input:
    file protein_db from file(params.protein_db)
    file mzXML_fragger from file(params.empty_mzxml)
    file fragger_params from file(params.fragger_params)

    output:
    file '*.pepindex' into msfraggerIndexOut
    
    script:
    """
    sed -i 's,num_threads = 0,num_threads = $params.fragger_threads,' $fragger_params
    sed -i 's,db_path,$protein_db,' $fragger_params

    java -XX:MaxRAMPercentage=80 -XX:+UseContainerSupport \
    -jar /usr/local/bin/MSFragger.jar $fragger_params $mzXML_fragger

    """
}


process msfraggerSearch {
    tag "$mzXML_fragger"
    
    cpus params.fragger_threads

    publishDir 'Results/MSFragger', mode: 'link'
    
    input:
    file protein_db from file(params.protein_db)
    file fragger_index from msfraggerIndexOut
    file mzXML_fragger from Channel.fromPath("${params.dda_folder}/*.mzXML").concat(mgf2mzxmlOut1)
    file fragger_params from file(params.fragger_params)

    output:
    set file("${mzXML_fragger}"), file('*.pepXML') into msfraggerSearchOut
    
    script:
    """
    sed -i 's,num_threads = 0,num_threads = $params.fragger_threads,' $fragger_params
    sed -i 's,db_path,$protein_db,' $fragger_params

    java -XX:MaxRAMPercentage=80 -XX:+UseContainerSupport \
    -jar /usr/local/bin/MSFragger.jar $fragger_params $mzXML_fragger
    """
}


//subsample_ratio = (1 / (file("${params.dda_folder}/*.mzXML").size() + mgf2mzxmlOut2.count().val) )
//subsample_ratio = (2 / file("${params.dda_folder}/*.mzXML").size())
subsample_ratio = 0.5
process msfraggerConvert {
    tag "$pepXML"
    
    input:
    set mzXML , pepXML from msfraggerSearchOut
    file unimod from file(params.msfraggerconvert_unimod)
    
    output:
    file '*_psms.tsv'
    file '*_subpsms.tsv' into msfraggerConvertPsmsOut
    file '*.peakpkl' into msfraggerConvertPklOut
    
    script:
    """
    easypqp convert --unimod $unimod \
    --pepxml $pepXML \
    --mzxml $mzXML \
    --subsample_fraction $subsample_ratio
    """
}


// Remove searches with less than 10 decoys or 10 target groups (required by pyprophetScore process)
// and duplicate resulting channel
msfraggerConvertPsmsFilteredOut = msfraggerConvertPsmsOut
    .filter{ it.readLines().grep( ~/.*True.*/ ).size() > 9 }
    .filter{ it.readLines().grep( ~/.*False.*/ ).size() > 9 }


msfraggerConvertPsmsFilteredOut.into{msfraggerConvertPsmsOut1; msfraggerConvertPsmsOut2}


process pyprophetMerge {
    tag "$subpsms"
    
    input:
    file subpsms from msfraggerConvertPsmsOut1.collect()

    output:
    file 'pyprophet_learn_Q0.tsv' into pyprophetMergeQ0Out
    file 'pyprophet_learn_Q1.tsv' into pyprophetMergeQ1Out
    file 'pyprophet_learn_Q2.tsv' into pyprophetMergeQ2Out
    file 'pyprophet_learn_Q3.tsv' into pyprophetMergeQ3Out
    
    script:
    """
    awk 'BEGIN {{ FS=\"\t\"; OFS=\"\t\" }} (FNR>1 && \$22==0) || NR==1 {{print \$0}}' $subpsms > pyprophet_learn_Q0.tsv; \
    awk 'BEGIN {{ FS=\"\t\"; OFS=\"\t\" }} (FNR>1 && \$22==1) || NR==1 {{print \$0}}' $subpsms > pyprophet_learn_Q1.tsv; \
    awk 'BEGIN {{ FS=\"\t\"; OFS=\"\t\" }} (FNR>1 && \$22==2) || NR==1 {{print \$0}}' $subpsms > pyprophet_learn_Q2.tsv; \
    awk 'BEGIN {{ FS=\"\t\"; OFS=\"\t\" }} (FNR>1 && \$22==3) || NR==1 {{print \$0}}' $subpsms > pyprophet_learn_Q3.tsv
    """
}


process pyprophetLearn {
    scratch 'ram-disk'
    stageInMode "copy"
    tag "$pyprophetMergeQ0Out - $pyprophetMergeQ1Out - $pyprophetMergeQ2Out - $pyprophetMergeQ3Out"

    input:
    file q0 from pyprophetMergeQ0Out
    file q1 from pyprophetMergeQ1Out
    file q2 from pyprophetMergeQ2Out
    file q3 from pyprophetMergeQ3Out    

    output:
    file "pyprophet_learn_Q0_weights.csv" into pyprophetLeanQ0Out
    file "pyprophet_learn_Q0_summary_stat.csv"
    file "pyprophet_learn_Q0_full_stat.csv"
    file "pyprophet_learn_Q0_scored.tsv"
    file "pyprophet_learn_Q0_report.pdf"
    file "pyprophet_learn_Q1_weights.csv" into pyprophetLeanQ1Out
    file "pyprophet_learn_Q1_summary_stat.csv"
    file "pyprophet_learn_Q1_full_stat.csv"
    file "pyprophet_learn_Q1_scored.tsv"
    file "pyprophet_learn_Q1_report.pdf"
    file "pyprophet_learn_Q2_weights.csv" into pyprophetLeanQ2Out
    file "pyprophet_learn_Q2_summary_stat.csv"
    file "pyprophet_learn_Q2_full_stat.csv"
    file "pyprophet_learn_Q2_scored.tsv"
    file "pyprophet_learn_Q2_report.pdf"
    file "pyprophet_learn_Q3_weights.csv" into pyprophetLeanQ3Out
    file "pyprophet_learn_Q3_summary_stat.csv"
    file "pyprophet_learn_Q3_full_stat.csv"
    file "pyprophet_learn_Q3_scored.tsv"
    file "pyprophet_learn_Q3_report.pdf"


    shell:
    '''
    LINES=$(wc -l !{q0} | cut -f1 -d' ')
    if [ "$LINES" -gt "1" ]
    then
        pyprophet score --in !{q0} --threads=!{params.pyprophet_learn_threads};
    fi

    LINES=$(wc -l !{q1} | cut -f1 -d' ')
    if [ "$LINES" -gt "1" ]
    then
        pyprophet score --in !{q1} --threads=!{params.pyprophet_learn_threads};
    fi

    LINES=$(wc -l !{q2} | cut -f1 -d' ')
    if [ "$LINES" -gt "1" ]
    then
        pyprophet score --in !{q2} --threads=!{params.pyprophet_learn_threads};
    fi

    LINES=$(wc -l !{q3} | cut -f1 -d' ')
    if [ "$LINES" -gt "1" ]
    then
        pyprophet score --in !{q3} --threads=!{params.pyprophet_learn_threads};
    fi
    '''
}


process pyprophetScore {
    scratch 'ram-disk'
    stageInMode "copy"
    tag "$subpsms"
    
    input:
    file subpsms from msfraggerConvertPsmsOut2
    file q0 from pyprophetLeanQ0Out
    file q1 from pyprophetLeanQ1Out
    file q2 from pyprophetLeanQ2Out
    file q3 from pyprophetLeanQ3Out
        
    output:
    file "*_subpsms_scored.tsv" into pyprophetScoreOut
    file "*_subpsms_summary_stat.csv"
    file "*_subpsms_full_stat.csv"
    file "*_subpsms_report.pdf"
    
    script:
    if( subpsms.contains('_Q1_') )
    """
    pyprophet score --in $subpsms --apply_weights=$q1 --pi0_lambda=$params.pyprophetscore_pi0_lambda
    """
    else if( subpsms.contains('_Q2_') )
    """
    pyprophet score --in $subpsms --apply_weights=$q2 --pi0_lambda=$params.pyprophetscore_pi0_lambda
    """
    else if( subpsms.contains('_Q3_') )
    """
    pyprophet score --in $subpsms --apply_weights=$q3 --pi0_lambda=$params.pyprophetscore_pi0_lambda

    """
    else
    """
    pyprophet score --in $subpsms --apply_weights=$q0 --pi0_lambda=$params.pyprophetscore_pi0_lambda
    """
}


// Create Spectral Library
process easypqp {
//    scratch 'true'
//    stageInMode "copy"
    tag "$psms"
    
    publishDir 'Results/easypqpLib', mode: 'link'
    
    input:
    file psms from pyprophetScoreOut.collect()
    file peakpkl from msfraggerConvertPklOut.collect()

    output:
    file "*_global_peaks.tsv" into easypqpOut
    file "pyprophet_peptide_report.pdf"
    file "pyprophet_protein_report.pdf"
    
    script:
    """
    easypqp library --psm_fdr_threshold=$params.easypqp_psm_fdr_threshold \
    --peptide_fdr_threshold=$params.easypqp_peptide_fdr_threshold \
    --protein_fdr_threshold=$params.easypqp_protein_fdr_threshold \
    --pi0_lambda=$params.easypqp_pi0_lambda \
    --peptide_plot=pyprophet_peptide_report.pdf \
    --protein_plot=pyprophet_protein_report.pdf \
    $psms \
    $peakpkl
    """
}


process globalTargetPqp {
    tag "$peaks"

    input:
    file peaks from easypqpOut.flatten()
    
    output:
    file "*_global_peaks_pqp.tsv" into assayGeneratorOut
    
    script:
    """
    OpenSwathAssayGenerator -in $peaks \
    -out `basename ${peaks} .tsv`_pqp.tsv  \
    -precursor_upper_mz_limit $params.globaltargetpqp_precursor_upper_mz_limit \
    -product_lower_mz_limit $params.globaltargetpqp_product_lower_mz_limit \
    -min_transitions $params.globaltargetpqp_min_transitions \
    -max_transitions $params.globaltargetpqp_max_transitions
    """
}


process globalCombinedPqp {
    tag "$global_targets"
    
    input:
    file global_targets from assayGeneratorOut.collect()

    output:
    file 'combined_global_pqp.tsv' into globalCombinedPqpOut
    
    script:
    """
    awk 'BEGIN {{ FS=\"\t\"; OFS=\"\t\" }} FNR>1 || NR==1 {{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$15,\$18,\$19,\$20,\$25,\$26,\$27,\$28,\$29}}' $global_targets > combined_global_pqp.tsv
    """
}


process globalCombinedDecoyPqp {
    scratch 'ram-disk'
    stageInMode "copy"
    tag "global_pqp"
    
    publishDir 'Results/easypqpLib', mode: 'link'

    input:
    file global_pqp from globalCombinedPqpOut

    output:
    file 'library.pqp'
    
    script:
    """
    OpenSwathDecoyGenerator -in $global_pqp -out library.pqp
    """
}

