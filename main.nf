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
    log.info "	                     v"
    log.info "	                +-----------+"
    log.info "	                | MSFragger |"
    log.info "	                +----+------+"
    log.info "	                     |"
    log.info "	                     v"
    log.info "	                +-----------+"
    log.info "	                | Xinteract |"
    log.info "	                +----+------+"
    log.info "	                     |"
    log.info "	                     v"
    log.info "	                +-----------+"
    log.info "	                |  easypqp  |"
    log.info "	                +----+------+"
    log.info "	                     |"
    log.info "	                     v"
    log.info "	        Consensus Spectral Library"
    log.info "	                     |"
    log.info "	                     v"
    log.info "	       +------------------------+"
    log.info "	       | OpenSwathAssayGenerator|"
    log.info "	       +------------------------+"
    log.info "	                     |"
    log.info "	                     v"
    log.info "	       +------------------------+"
    log.info "	       | OpenSwathDecoyGenerator|"
    log.info "	       +------------------------+"
    log.info "	                     |"
    log.info "	                     v"
    log.info "	        PQP assay library with decoys"
    log.info ""
    log.info ""
    log.info "Options:"
    log.info " ** General parameters **"
    log.info "  --help:		Show this message and exit"
    log.info "  --dda_folder:	Path to DDA files ($params.dda_folder)"
    log.info "  --dia_folder:	Path to DIA files ($params.dia_folder)"
    log.info "  --protein_db:	Protein database to search against ($params.protein_db)"
    log.info ""
    log.info " ** DIA-Umpire parameters **"
    log.info "  --diau_threads:		Number of threads to use for DIA-Umpire SE ($params.diau_threads)"
    log.info "  --diau_se_params:	DIA-Umpire parameter file for SE ($params.diau_se_params)"
    log.info ""
    log.info " ** MSFragger parameters **"
    log.info "  --fragger_threads:		Number of threads to use for MS-Fragger search ($params.fragger_threads)"
    log.info "  --fragger_params:		MS-Fragger search parameter file ($params.fragger_params)"
    log.info "  --msfraggerconvert_unimod:	Unimod modifications file ($params.msfraggerconvert_unimod)"
    log.info ""
    log.info " ** TPP parameters **"
    log.info "  --tpp:		xinteract parameters for tpp analysis ($params.tpp)"
    log.info "  --decoy:	Decoy prefix in protein database ($params.decoy)"
    log.info ""
    log.info " ** easypqp library generation parameters **"
    log.info "  --easypqp_pi0_lambda:			pi0 lamdba interval for library generation ($params.easypqp_pi0_lambda)"
    log.info "  --easypqp_psm_fdr_threshold:		PSM FDR for library generation ($params.easypqp_psm_fdr_threshold)"
    log.info "  --easypqp_peptide_fdr_threshold:	peptide FDR for library generation ($params.easypqp_peptide_fdr_threshold)"
    log.info "  --easypqp_protein_fdr_threshold:	protein FDR for library generation ($params.easypqp_protein_fdr_threshold)"
    log.info ""
    log.info " ** oswAssayGeneartor parameters **"
    log.info "  --oswAssayGenerator_mode:			OSW or IPF ($params.oswAssayGenerator_mode)"
    log.info "  --oswAssayGenerator_product_lower_mz_limit	low m/z boundary for library inclusion ($params.oswAssayGenerator_product_lower_mz_limit)"
    log.info "  --oswAssayGenerator_precursor_upper_mz_limit	high m/z boundary for library inclusion ($params.oswAssayGenerator_precursor_upper_mz_limit)"
    log.info "  --oswAssayGenerator_min_transitions		min number of transitions for library inclusion ($params.oswAssayGenerator_min_transitions)"
    log.info "  --oswAssayGenerator_max_transitions		max number of transitions for library inclusion ($params.oswAssayGenerator_max_transitions)"
    log.info "  --oswAssayGenerator_unimod			Unimod modifications file for IPF ($params.oswAssayGenerator_unimod)"
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


mgf2mzxmlOut.into{mgf2mzxmlOut1; mgf2mzxmlOut2}


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


process peptideProphet {
    tag "$pepxml"
    
    input:
    file pepxml from msfraggerSearchOutPep

    output:
    file '*.pep.xml' into peptideProphetOut

    script:
    """
    InteractParser ${pepxml.baseName}.pep.xml ${pepxml} DECOY=$params.decoy ACCMASS PPM NONPARAM DECOYPROBS 
    PeptideProphetParser ${pepxml.baseName}.pep.xml 
    """
}


process iprophet {
    cpus params.iprophet_threads
    
    publishDir 'Resutls/iProphet', mode: 'link'

    input:
    file pepxmls from peptideProphetOut.collect()

    output:
    file 'iprophet.pep.xml' into iProphetOut

    script:
    """
    InterProphetParser THREADS=$params.iprophet_threads DECOY=$params.decoy ${pepxmls} iprophet.pep.xml
    """
}



// This needs to run once for each mzXML file we have (even if we have
// pooled all search results into a single pepXML file)
process easypqpConvert {
    memory = 50.GB
    
    input:
    file pepxml from iProphetOut
    file mzxml from Channel.fromPath("${params.dda_folder}/*.mzXML").concat(mgf2mzxmlOut2)
    file unimod from file(params.unimod)

    output:
    file '*_psms.tsv' into pepxmlConvertPsmsOut
    file '*.peakpkl' into pepxmlConvertPklOut    
    
    script:
    """
    easypqp convert --unimod $unimod --pepxml $pepxml --psms ${mzxml.baseName}_psms.tsv --spectra $mzxml --peaks ${mzxml.baseName}.peakpkl
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
    --rt_lowess_fraction=$params_lowess_fraction \
    --pi0_lambda=$params.easypqp_pi0_lambda \
    --peptide_plot=pyprophet_peptide_report.pdf \
    --protein_plot=pyprophet_protein_report.pdf \
    $psms \
    $peakpkl
    """
}


process oswAssayGenerator {
    tag "$library"

    publishDir 'Results/easypqpLib', mode: 'link'
    
    input:
    file library from easypqpOut
    
    output:
    file "pqp.tsv" into assayGeneratorOut
    
    script:
    if( params.oswAssayGenerator_mode == 'OSW' )
        """
        OpenSwathAssayGenerator -in $library \
        -out pqp.tsv  \
        -precursor_upper_mz_limit $params.oswAssayGenerator_precursor_upper_mz_limit \
        -product_lower_mz_limit $params.oswAssayGenerator_product_lower_mz_limit \
        -min_transitions $params.oswAssayGenerator_min_transitions \
        -max_transitions $params.oswAssayGenerator_max_transitions
        """
    else if( params.oswAssayGenerator_mode == 'IPF' )
	"""
        OpenSwathAssayGenerator -in $library
        -out pqp.tsv
        -enable_ipf 
        -unimod_file $params.unimod
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

