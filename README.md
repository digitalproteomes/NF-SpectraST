# Nextflow workflow for Peptide Query Parameter generation

This workflow generates a spectral library in PQP format for OpenSWATH analysis.
Both DDA and DIA files (via DIAUmpire) are used for the library generation.
The results from multiple search engines are combined via InterProphet.
For library generation the list of proteins derived from ProteinProphet is filtered at 1%FDR.


## Usage

The workflow takes the following parameters:
* --help:              show this message and exit
* --dia_folder:        file to be DIAUmpired (default: Data/DIA)
* --diau_se_params:    parameter file for DIAUmpire (default: Params/diaumpire_se.params)
* --dda_folder:        files to be searched (default: Data/DDA)
* --comet_params:      comet parameter file (default: Params/comet.params)
* --protein_db:        comet parameter file (default: Results/Databases/proteome.fasta)
* --tpp:               TPP options (default: -OAPdlIw -PPM)
* --decoy:             decoy prefix (default: DECOY_)
* --rt_file:           RT normalization file (default: Params/irtkit.txt
* --st_fragmentation:  fragmentation type to use in SpectraST (default: -cIHCD)
* --st_fix_mods:       sed script file for conversion of TPP modifications to UniMod modifications (default: Params/fix_mods.txt)
* --swath_window_file: file with the definition of the SWATH windows (default: Params/swaths64.txt)

Example usage:

```bash
nextflow run digitalproteomes/NF-SpectraST
```
At the end of the workflow the prepared PQP library will be called SpecLib_opt_dec.pqp and can eb foudn in the *Results/SpectraST* folder.