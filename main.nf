if(params.help) {

    log.info "PQP generation workflow"
    log.info "----------------------"
    log.info ""
    log.info "Options:"
    log.info " ** General parameters **"
    log.info "  --help:		Show this message and exit"
    log.info "  --dda_folder:	Path to DDA files ($params.dda_folder)."
    log.info "                  NOTE: pseudoDDA files will be saved here."
    log.info "  --dia_folder:	Path to DIA files ($params.dia_folder)"
    log.info "  --protein_db:	Protein database to search against ($params.protein_db)"
    log.info "  --unimod:	Unimod file location ($params.unimod)"
    log.info ""
    log.info " ** DIA-Umpire parameters **"
    log.info "  --diau_threads:		Number of threads to use for DIA-Umpire SE ($params.diau_threads)"
    log.info "  --diau_se_params:	DIA-Umpire parameter file for SE ($params.diau_se_params)"
    log.info ""
    log.info " ** Comet parameters **"
    log.info "  --comet_threads:	Threads to use for comet search ($params.comet_threads)"
    log.info "  --comet_params:		The comet params file ($params.comet_params)"    
    log.info ""
    log.info " ** TPP parameters **"
    log.info "  --peptideProphet_params:	PeptideProphet parameters ($params.peptideProphet_params)"
    log.info "  --decoy:			Decoy prefix in protein database ($params.decoy)"
    log.info "  --iprophet_threads:		Threads to use for iprophet ($params.iprophet_threads)"
    log.info ""
    log.info " ** easypqp library generation parameters **"
    log.info "  --easypqp_rt_lowess_fraction:	Fraction od data point to use for RT regression ($params.easypqp_rt_lowess_fraction)"
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


dda_filenames = Channel
    .fromPath("$params.dda_folder/*.mzXML")
    .flatMap{ it -> "${it.baseName}".replaceAll(/_Q[0-9]/,"") }

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
    val diaumpire_done from dda_filenames.collect()

    when:
    ! diaumpire_done.contains(dia_file.baseName)
        
    output:
    file '*.mgf' into diaUmpireOut

    script:
    """
    sed -i 's,Thread = 0,Thread = $params.diau_threads,' "$diau_se_params"
    java -jar -Xmx32G /DIA_Umpire_SE.jar "$dia_file" "$diau_se_params"
    """

}


process mgf2mzxml {
    // Convert mgfs generated by DIAUmpire to mzXML
    tag "$mgf_file"

    publishDir params.dda_folder, mode: 'link'

    cache 'lenient'

    input:
    // For each DIA file DIA-Umpire will generate 3 pseudo mgf files
    file mgf_file from diaUmpireOut.flatten()

    output:
    file '*.mzXML' into mgf2mzxmlOut
    
    """
    msconvert "$mgf_file" --mzXML
    """
}


mgf2mzxmlOut.into{mgf2mzxmlOut1; mgf2mzxmlOut2; mgf2mzxmlOut3}


process cometSearch {
    tag "$mzXML"
    
    cpus params.comet_threads
    // Human, 3 variable mods, semi, 2 missed cleavages and some margin for safety
    memory 30.GB

    publishDir 'Results/Comet', mode: 'link'

    cache 'lenient'
    
    input:
    file mzXML from Channel.fromPath("${params.dda_folder}/*.mzXML").concat(mgf2mzxmlOut1)
    file comet_params from file(params.comet_params)
    file protein_db from file(params.protein_db)

    output:
    file '*.pep.xml' into cometSearchPepxmlOut
    file mzXML into cometSearchMzxmlOut

    """
    # Set proteins DB
    sed -i s,db_path,$protein_db, $comet_params
    sed -i 's,num_threads = 0,num_threads = ${params.comet_threads},' $comet_params

    comet -P$comet_params $mzXML
    sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g ${mzXML.simpleName}.pep.xml
    sed -ri 's|<search_database local_path="|<search_database local_path="${workflow.launchDir}/Results/Comet/|'g ${mzXML.simpleName}.pep.xml
    """
}


// Run PeptideProphet on each search output individually
process peptideProphet {
    tag "$pepxml"

    publishDir 'Results/PeptideProphet', mode: 'link'
    
    input:
    file protein_db from file(params.protein_db)
    file pepxml from cometSearchPepxmlOut
    // TODO: figure out how to link only the required mzXML file, rather than all.
    file mzxml from cometSearchMzxmlOut.collect()

    output:
    file '*.pep.xml' into peptideProphetOut
    file '*.pep.xml.pIstats'
    file '*.pep-MODELS.html'
    
    script:
    """
    InteractParser "${pepxml.baseName}".pep.xml "${pepxml}"
    PeptideProphetParser "${pepxml.baseName}".pep.xml DECOY=$params.decoy $params.peptideProphet_params
    tpp_models.pl "${pepxml.baseName}".pep.xml
    """
}


// Run iProphet on all PeptideProphet files at once.
process iProphet {
    tag "$pepxmls"
    
    cpus params.iprophet_threads
    
    publishDir 'Results/iProphet', mode: 'link'

    input:
    file pepxmls from peptideProphetOut.collect()

    output:
    file 'iprophet.pep.xml' into iProphetOut
    file 'iprophet.pep-MODELS.html' into iProphetModelOut
    
    script:
    """
    InterProphetParser THREADS=$params.iprophet_threads DECOY=$params.decoy ${pepxmls} iprophet.pep.xml
    tpp_models.pl iprophet.pep.xml
    """
}


iProphetOut.into{ iProphetOut1; iProphetOut2}
iProphetModelOut.into{ iProphetModelOut1; iProphetModelOut2 }


// Run ProphetProphet on the output of iProphet
process proteinProphet {
    tag "$pepxml"
    
    publishDir 'Results/ProteinProphet', mode: 'link'
    
    input:
    file pepxml from iProphetOut1
    file protein_db from file(params.protein_db)

    output:
    file '*.prot.xml' into proteinProphetOut
    file '*.prot-MODELS.html' into proteinProphetModelOut
    
    script:
    """
    ProteinProphet $pepxml iprophet.prot.xml IPROPHET
    tpp_models.pl iprophet.prot.xml
    """
}


// Get a list of proteins filtered at 1%FDR.
//
// This will be used to filter the PSMs that go into library
// generation using ProteinProphet FDR, rather than pyprophet FDR.
process getProteinList {
    tag "$protxml"

    publishDir 'Results/ProteinProphet', mode: 'link'

    input:
    file protxml from proteinProphetOut
    file protxml_models from proteinProphetModelOut

    output:
    file 'protein_list.tsv' into getProteinListOut

    script:
    """
    PROB=\$(get_prophet_prob.py -i $protxml_models)
    xsltproc --param p_threshold \$PROB /usr/local/bin/get_protein_list.xsl $protxml > protein_list.tsv
    """
}


// This needs to run once for each mzXML file we have (even if we were to
// pool all search results into a single pepXML file).
//
// One psmpkl and peakpkl file will be generated for each mzXML passed
// to the process. The PSMs will be extracted from the single iprophet
// output file (the associated probability will be from iProphet). The
// peaks will come from the corresponding mzXML file.
process easypqpConvert {
    tag "$mzxml"
    
    memory = 10.GB
    
    input:
    file pepxml from iProphetOut2
    file mzxml from Channel.fromPath("${params.dda_folder}/*.mzXML").concat(mgf2mzxmlOut3)
    file unimod from file(params.unimod)

    output:
    tuple file('*.psmpkl'), file('*.peakpkl') into easypqpConvert
    
    script:
    """
    easypqp convert --unimod $unimod --pepxml $pepxml --psms ${mzxml.baseName}.psmpkl --spectra $mzxml --peaks ${mzxml.baseName}.peakpkl
    """
}


process getPepPThreshold {
    tag "$pepxml_models"
    
    input:
    file pepxml_models from iProphetModelOut1

    output:
    env PROB into iProphetPThresholdOut

    script:
    """
    PROB=\$(get_prophet_prob.py -i $pepxml_models)
    """
}


process filterPqp {
    tag "$psmfile"

    memory = 50.GB
    
    input:
    tuple file(psmfile), file(peakfile) from easypqpConvert
    file pepxml_models from iProphetModelOut2
    env PROB from iProphetPThresholdOut
    file protein_list from getProteinListOut

    output:
    file '*_filtered.psmpkl' into filterPqpPsmOut
    file '*_filtered.peakpkl' into filterPqpPeakOut
    
    script:
    """
    filterpqp.py -s $psmfile -k $peakfile -l $protein_list -p \$PROB
    """
}

filterPqpPsmOut.into{filterPqpPsmOut1; filterPqpPsmOut2}
filterPqpPeakOut.into{filterPqpPeakOut1; filterPqpPeakOut2}

// Create library for RT alignment
//
// We are only going to keep IDs from Q1 to minimize the risk of
// reducing the library to IDs that are not seen in the DIA runs.
filterPqpPsmRTOut = filterPqpPsmOut1.filter( ~/.*_Q1_filtered.psmpkl/ )
filterPqpPeakRTOut = filterPqpPeakOut1.filter( ~/.*_Q1_filtered.peakpkl/ )
process easypqpRT {
    tag "$psms"
    
    publishDir 'Results/easypqpLib', mode: 'link'
    
    input:
    file psms from filterPqpPsmRTOut.collect()
    file peakpkl from filterPqpPeakRTOut.collect()

    output:
    file "library_RT.tsv" into easypqpRTOut
    file "pyprophet_peptide_report_RT.pdf" optional true // Not generated with --nofdr
    file "pyprophet_protein_report_RT.pdf" optional true // Not generated with --nofdr
    file "easypqp_rt_alignment*.pdf"
    file "*_run_peaks.tsv"
    
    script:
    """
    easypqp library --out library_RT.tsv \
    --rt_lowess_fraction=$params.easypqp_rt_lowess_fraction \
    --peptide_plot=pyprophet_peptide_report_RT.pdf \
    --protein_plot=pyprophet_protein_report_RT.pdf \
    --consensus \
    --nofdr \
    $psms \
    $peakpkl
    """
}

// Assay generation for RT alignment library.
//
// NOTE we don't need decoys in this one since it is only used for
// alignment and not for searches
process oswAssayGeneratorRT {
    tag "$library"

    publishDir 'Results/easypqpLib', mode: 'link'
    
    input:
    file library from easypqpRTOut
    
    output:
    file "library_targets_RT.pqp"
    
    script:
    if( params.oswAssayGenerator_mode == 'OSW' )
        """
        OpenSwathAssayGenerator -in $library \
        -out library_targets_RT.pqp  \
        -precursor_lower_mz_limit $params.oswAssayGenerator_precursor_lower_mz_limit \
        -precursor_upper_mz_limit $params.oswAssayGenerator_precursor_upper_mz_limit \
        -product_lower_mz_limit $params.oswAssayGenerator_product_lower_mz_limit \
        -product_upper_mz_limit $params.oswAssayGenerator_product_upper_mz_limit \
        -min_transitions $params.oswAssayGenerator_min_transitions \
        -max_transitions $params.oswAssayGenerator_max_transitions\
        -swath_windows_file $params.oswAssayGenerator_swath_windows_file \
        -unimod $params.unimod
        """
    else if( params.oswAssayGenerator_mode == 'IPF' )
	"""
        OpenSwathAssayGenerator -in $library
        -out library_targets_RT.pqp
        -enable_ipf 
        -unimod_file $params.unimod
        -disable_identification_ms2_precursors 
        -disable_identification_specific_losses
        """
    else
	error "Invalid assay generation mode: ${params.oswAssayGenerator_mode}"
}

// Create Spectral Library
process easypqp {
    tag "$psms"
    
    publishDir 'Results/easypqpLib', mode: 'link'
    
    input:
    file psms from filterPqpPsmOut2.collect()
    file peakpkl from filterPqpPeakOut2.collect()

    output:
    file "library.tsv" into easypqpOut
    file "pyprophet_peptide_report.pdf" optional true // Not generated with --nofdr
    file "pyprophet_protein_report.pdf" optional true // Not generated with --nofdr
    file "easypqp_rt_alignment*.pdf"
    file "*_run_peaks.tsv"
    
    script:
    """
    easypqp library --out library.tsv \
    --rt_lowess_fraction=$params.easypqp_rt_lowess_fraction \
    --peptide_plot=pyprophet_peptide_report.pdf \
    --protein_plot=pyprophet_protein_report.pdf \
    --nofdr \
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
    file "library_targets.pqp" into assayGeneratorOut
    
    script:
    if( params.oswAssayGenerator_mode == 'OSW' )
        """
        OpenSwathAssayGenerator -in $library \
        -out library_targets.pqp  \
        -precursor_lower_mz_limit $params.oswAssayGenerator_precursor_lower_mz_limit \
        -precursor_upper_mz_limit $params.oswAssayGenerator_precursor_upper_mz_limit \
        -product_lower_mz_limit $params.oswAssayGenerator_product_lower_mz_limit \
        -product_upper_mz_limit $params.oswAssayGenerator_product_upper_mz_limit \
        -min_transitions $params.oswAssayGenerator_min_transitions \
        -max_transitions $params.oswAssayGenerator_max_transitions \
        -swath_windows_file $params.oswAssayGenerator_swath_windows_file \
        -unimod $params.unimod
        """
    else if( params.oswAssayGenerator_mode == 'IPF' )
	"""
        OpenSwathAssayGenerator -in $library
        -out library_targets.pqp
        -enable_ipf 
        -unimod_file $params.unimod
        -disable_identification_ms2_precursors 
        -disable_identification_specific_losses
        """
    else
	error "Invalid assay generation mode: ${params.oswAssayGenerator_mode}"
}


process oswDecoyGenerator {
    tag "$pqp"
    
    publishDir 'Results/easypqpLib', mode: 'link'

    input:
    file pqp from assayGeneratorOut

    output:
    file 'library.pqp'
    
    script:
    """
    OpenSwathDecoyGenerator -in $pqp \
    -out library.pqp
    """
}

