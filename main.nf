// Create pseudoMS2 with DIAUmpire
// Searches with Comet and Tandem
// Pooled PeptideProphet for Comet and Tandem separately
// InterProphet on pooled PeptideProphets
// ProteinProphet on InterProphet
// Mayu on InterProphet
// Filter 1%FDR Protein -> SpectraST

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
    // Generate pseudoDDA files for each DIA file
    input:
    file dia_file from file("${params.dia_folder}/*.mzXML")
    file diau_se_params from file(params.diau_se_params)

    output:
    file '*.mgf' into diaUmpireOut
    
    """
    dia_umpire_se $dia_file $diau_se_params
    """
}


process mgf2mzxml {
    // Convert mgfs generated by DIAUmpire to mzXML
    publishDir 'Results/DIAUmpire'

    input:
    // For each DIA file DIA-Umpire will generate 3 pseudo mgf files
    file mgf_file from diaUmpireOut.flatten()

    output:
    file '*.mzXML' into pDdaFiles
    
    """
    msconvert $mgf_file --mzXML
    """
}


pDdaFiles.into{pDdaFiles1; pDdaFiles2}

process cometSearch {
    // Search all mzXML files in $params.dda_folder and diaUmpire
    // extracted ones with Comet
    publishDir 'Results/Comet'
    
    input:
    file mzXML_comet from Channel.fromPath("${params.dda_folder}/*.mzXML").concat(pDdaFiles1)
    file comet_params from file(params.comet_params)
    file protein_db from file(params.protein_db)

    output:
    file '*.pep.xml' into cometOut
    file mzXML_comet

    """
    comet $mzXML_comet
    """
}


process pooledCometTpp {
    // Interact together all comet searches
    publishDir 'Results/Comet'
    
    input:
    file pepxmls from cometOut.collect()
    file protein_db from file(params.protein_db)

    output:
    file 'comet_merged.pep.xml' into tppPepOut_comet
    file 'comet_merged.pep-MODELS.html'
    file 'comet_merged.pep.xml.index'
    file 'comet_merged.pep.xml.pIstats'

    // xinteract and refactor links in prot.xml 
    """
    xinteract $params.tpp -d$params.decoy -Ncomet_merged.pep.xml $pepxmls
    """
}


process tandemSearch {
    // Search all mzXML files in $params.dda_folder with Tandem
    publishDir 'Results/Tandem'
    
    input:
    file mzXML_tandem from Channel.fromPath("${params.dda_folder}/*.mzXML").concat(pDdaFiles2)
    file tandem_params from file(params.tandem_params)
    file tandem_default_params from file(params.tandem_default_params)
    file taxonomy_params from file(params.tandem_taxonomy)
    file protein_db from file(params.protein_db)

    output:
    file '*_raw.pep.xml' into tandemOut
    // For Tandem we need to stage the mzXML files to the next step
    // (pooledTandemTpp) as well to enable spectrum visualization in
    // PepXMLViewer
    file mzXML_tandem into tandemOutMzxml
    file 'tandem_params*'
    file 'taxonomy.xml'

    script:
    """
    FNAME=\$(basename $mzXML_tandem .mzXML)
    
    # Set input and output
    sed -i s,full_mzXML_filepath,$mzXML_tandem, $tandem_params
    sed -i s,full_tandem_output_path,\${FNAME}.xtan.xml, $tandem_params
    mv $tandem_params tandem_params_\${FNAME}.xml

    # Set proteins DB
    sed -i s,db_path,$protein_db, $taxonomy_params

    # Perform search
    tandem tandem_params_\${FNAME}.xml 
    
    # Convert output to pep.xml
    Tandem2XML \${FNAME}.xtan.xml > \${FNAME}_raw.pep.xml
    """
}


process pooledTandemTpp {
    // Interact together all tandem searches
    publishDir 'Results/Tandem'
    
    input:
    file pepxmls from tandemOut.collect()
    file mzXML_tandem from tandemOutMzxml.collect()
    file protein_db from file(params.protein_db)

    output:
    file 'tandem_merged.pep.xml' into tppPepOut_tandem
    file 'tandem_merged.pep-MODELS.html'
    file 'tandem_merged.pep.xml.index'
    file 'tandem_merged.pep.xml.pIstats'

    // xinteract and refactor links in prot.xml 
    """
    xinteract $params.tpp -d$params.decoy -Ntandem_merged.pep.xml $pepxmls
    """
}


process interProphet {
    // Combine search engine results via InterProphet
    publishDir 'Results/SpectraST'
    
    input:
    file pepxmls from tppPepOut_comet
	.merge(tppPepOut_tandem)
	.collect()

    output:
    file 'iprophet.pep.xml' into interPepOut
    file 'iprophet.pep-MODELS.html'
    
    script:
    """
    InterProphetParser THREADS=8 $pepxmls iprophet.pep.xml
    tpp_models.pl iprophet.pep.xml
    """
}


interPepOut.into{interPepOut1; interPepOut2; interPepOut3; interPepOut4}


process proteinProphet {
    // Run ProteinProphet on combined search results
    publishDir 'Results/SpectraST'
    
    input:
    file pepxml from interPepOut1
    file protein_db from file(params.protein_db)

    output:
    file 'iprophet.prot-MODELS.html'
    file 'iprophet.prot.xml' into tppProtOut
    
    script:
    """
    RefreshParser $pepxml $protein_db
    ProteinProphet $pepxml iprophet.prot.xml IPROPHET
    tpp_models.pl iprophet.prot.xml
    sed -ri 's|/work/.{2}/.{30}|/Results/SpectraST|'g iprophet.prot.xml
    """
}


process mayu {
    // For each TPP analysis run Mayu
    publishDir 'Results/SpectraST'

    input:
    file pepxml from interPepOut2
    file protein_db from file(params.protein_db)
    
    output:
    file 'mayu_iprophet.pep.xml_main_1.07.txt'
    file 'mayu_iprophet.pep.xml_main_1.07.csv' into mayuOut
    
    """
    Mayu.pl -A $pepxml \
    -C $protein_db \
    -E $params.decoy \
    -M mayu_$pepxml \
    -P pepFDR=0.01:1
    """
}


process tppStat {
    // For each TPP analysis run calctppstat
    publishDir 'Results/SpectraST'
    
    input:
    file pepxml from interPepOut3
    file protxml from tppProtOut

    output:
    file '*.summary.txt'
    
    """
    /usr/local/tpp/cgi-bin/calctppstat.pl -i $pepxml -d $params.decoy --full > ${pepxml}.summary.txt
    """
}


process parseMayu {
    input:
    file mayu_csv from mayuOut

    output:
    stdout into parseMayuOut

    script:
    """
    parse_mayu.py $mayu_csv
    """
}


process spectraST {
    // Assemble consensus spectral library and export it to .mrm with
    // PTMs converted to UniMod
    publishDir 'Results/SpectraST'
    
    input:
    file pepxml from interPepOut4
    val probability from parseMayuOut
    file irt from file(params.rt_file)
    file fix_mods from file(params.st_fix_mods)

    output:
    file "SpecLib.splib"
    file "SpecLib.sptxt"
    file "SpecLib.pepidx"
    file "SpecLib_cons.splib"
    file "SpecLib_cons.sptxt"
    file "SpecLib_cons.pepidx"
    file "spectrast.log"
    file "SpecLib_cons_conv.mrm" into spectraST
    
    script:
    """
    spectrast -cNSpecLib $params.st_fragmentation -cf'Protein!~$params.decoy' -cP$probability -c_IRT$irt -c_IRR $pepxml
    spectrast -cNSpecLib_cons $params.st_fragmentation -cAC SpecLib.splib
    spectrast -cNSpecLib_cons_conv $params.st_fragmentation -cM SpecLib_cons.splib
    sed -i -f $fix_mods SpecLib_cons_conv.mrm
    """
}


process assayGenerator {
    // Optimize assays
    publishDir 'Results/SpectraST'

    input:
    file mrm_lib from spectraST
    file swath_window_file from file(params.swath_window_file)
    
    output:
    file "SpecLib_opt.pqp" into assayGeneratorOut
    
    script:
    """
    OpenSwathAssayGenerator -in $mrm_lib -out SpecLib_opt.pqp -swath_windows_file $swath_window_file
    """
}


process decoyGenerator {
    // Append decoys
    publishDir "Results/SpectraST"
    
    input:
    file spec_lib from assayGeneratorOut

    output:
    file "SpecLib_opt_dec.pqp"
    
    script:
    """
    OpenSwathDecoyGenerator -in $spec_lib -out SpecLib_opt_dec.pqp
    """
}

