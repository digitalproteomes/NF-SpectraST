// Create pseudoMS2 with DIAUmpire
// Searches with Comet and Tandem
// Pooled PeptideProphet for Comet and Tandem separately
// InterProphet on pooled PeptideProphets
// ProteinProphet on InterProphet
// Mayu on InterProphet

// TODO update help

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

// TODO: Implement user defined publishDir

// NOTE if you already have run DIA-Umpire separately, you can add the
// pseudo-DDA files to the DDA folder and this step will be skipped
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


mgf2mzxmlOut.into{mgf2mzxmlOut1; mgf2mzxmlOut2; mgf2mzxmlOut3}


process msfraggerSearch {
    tag "$mzXML_fragger"
    
    cpus params.fragger_threads

    publishDir 'Results/MSFragger', mode: 'link'
    
    input:
    file protein_db from file(params.protein_db)
//    file fragger_index from msfraggerIndexOut
    file mzXML_fragger from Channel.fromPath("${params.dda_folder}/*.mzXML").concat(mgf2mzxmlOut1).collect()
    file fragger_params from file(params.fragger_params)

    output:
    //set file("${mzXML_fragger}"), file('*.pepXML') into msfraggerSearchOut
    file '*.pepXML' into msfraggerSearchOutPep
    
    
    script:
    """
    sed -i 's,num_threads = 0,num_threads = $params.fragger_threads,' $fragger_params
    sed -i 's,db_path,$protein_db,' $fragger_params

    java -XX:MaxRAMPercentage=80 -XX:+UseContainerSupport \
    -jar /usr/local/bin/MSFragger.jar $fragger_params $mzXML_fragger
    """
}


// TODO: Consider if this step should be split for different Qs

// xinteract all search results as one
process pooledTpp {
    publishDir 'Results/fragger', mode: 'link'
    
    input:
    file pepxmls from msfraggerSearchOutPep.collect()
    file mzXML_fragger from Channel.fromPath("${params.dda_folder}/*.mzXML").concat(mgf2mzxmlOut2).collect()
    file protein_db from file(params.protein_db)

    output:
    file 'fragger_merged.pep.xml' into tppPepOut
    file 'fragger_merged.pep-MODELS.html'
    file 'fragger_merged.pep.xml.index'
    file 'fragger_merged.pep.xml.pIstats'
    file 'fragger_merged.prot-MODELS.html'
    file 'fragger_merged.prot.xml' into tppProtOut
    file(protein_db) // Required for ProteinProphet visualization

    // xinteract and refactor links in prot.xml 
    """
    xinteract $params.tpp -d$params.decoy -Nfragger_merged.pep.xml $pepxmls
    sed -ri 's|/work/.{2}/.{30}|/Results/Fragger|'g fragger_merged.prot.xml
    """
}


// This needs to run once for each mzXML file we have (even if we have
// pooled all search results into a single pepXML file)
process easypqpConvert {
    input:
    file pepxml from tppPepOut
    file mzxml from Channel.fromPath("${params.dda_folder}/*.mzXML").concat(mgf2mzxmlOut3)
    file unimod from file(params.msfraggerconvert_unimod)

    output:
    file '*_psms.tsv' into pepxmlConvertPsmsOut
    file '*.peakpkl' into pepxmlConvertPklOut    
    
    script:
    """
    easypqp convert --pepxml $pepxml --mzxml $mzxml --unimod $unimod
    """
}


// Create Spectral Library
process easypqp {
//    scratch 'true'
//    stageInMode "copy"
    tag "$psms"
    
    publishDir 'Results/easypqpLib', mode: 'link'
    
    input:
    file psms from pepxmlConvertPsmsOut.collect()
    file peakpkl from pepxmlConvertPklOut.collect()

    output:
    file "library.tsv" into easypqpOut
    file "pyprophet_peptide_report.pdf"
    file "pyprophet_protein_report.pdf"
    
    script:
    """
    easypqp library --out library.tsv \
    --psm_fdr_threshold=$params.easypqp_psm_fdr_threshold \
    --peptide_fdr_threshold=$params.easypqp_peptide_fdr_threshold \
    --protein_fdr_threshold=$params.easypqp_protein_fdr_threshold \
    --pi0_lambda=$params.easypqp_pi0_lambda \
    --peptide_plot=pyprophet_peptide_report.pdf \
    --protein_plot=pyprophet_protein_report.pdf \
    $psms \
    $peakpkl
    """
}


process oswAssayGenerator {
    tag "$peaks"

    publishDir 'Results/easypqpLib', mode: 'link'
    
    input:
    file library from easypqpOut
    
    output:
    file "pqp.tsv" into assayGeneratorOut
    
    script:
    if( $params.oswAssayGenerator_mode == 'OSW' )
        """
        OpenSwathAssayGenerator -in $library \
        -out pqp.tsv  \
        -precursor_upper_mz_limit $params.oswAssayGenerator_precursor_upper_mz_limit \
        -product_lower_mz_limit $params.oswAssayGenerator_product_lower_mz_limit \
        -min_transitions $params.oswAssayGenerator_min_transitions \
        -max_transitions $params.oswAssayGenerator_max_transitions
        """
    else if( $params.oswAssayGenerator_mode == 'IPF' )
	"""
        OpenSwathAssayGenerator -in $library
        -out pqp.tsv
        -enable_ipf 
        -unimod_file $params.oswAssayGenerator_unimod
        -disable_identification_ms2_precursors 
        -disable_identification_specific_losses
        """
    else
	error "Invalid assay generation mode: ${params.mode}"
}


process oswDecoyGenerator {
    scratch 'ram-disk'
    stageInMode "copy"
    tag "$pqp"
    
    publishDir 'Results/easypqpLib', mode: 'link'

    input:
    file pqp from assayGeneratorOut

    output:
    file 'library.pqp'
    
    script:
    """
    OpenSwathDecoyGenerator -in $pqp -out library.pqp
    """
}

