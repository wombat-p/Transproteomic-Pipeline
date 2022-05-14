#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/tpp_workflow
========================================================================================
(NOT YET A nf-core!)

nextflow.enable.dsl=1

#### Homepage / Documentation
----------------------------------------------------------------------------------------
*/

import groovy.json.JsonOutput

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
    nextflow run main.nf  -profile docker

    For a test run with pre-given data, use:
    nextflow run main.nf -profile docker, test


    Mandatory arguments:
      --raws                            Path to input data (must be surrounded with quotes)
      --fasta                Fasta file for database search
      -profile                          Configuration profile to use. Can use multiple (comma separated)
                                        Available: standard, conda, docker, singularity, awsbatch, test

    Mass Spectrometry Search:
      --precursor_mass_tolerance        Mass tolerance of precursor mass (ppm)
      --fragment_mass_tolerance         Mass tolerance of fragment mass bin (da)
      --digest_mass_range               Mass range of peptides considered for matching
      --activation_method               Fragmentation method ('ALL', 'CID', 'ECD', 'ETD', 'PQD', 'HCD', 'IRMPD')
      --enzyme                          Enzymatic cleavage ('unspecific cleavage', 'Trypsin', see comet enzymes)
      --miscleavages                    Number of allowed miscleavages
      --number_mods                     Maximum number of modifications of PSMs
      --fixed_mods                      Not working yet. Fixed modification is always 'Carbamidomethyl (C)'
      --variable_mods                   Only variable modifications 'Oxidation of M' 'Acetylation of protein N-term' and 'Phosphorylation of STY' allowed. Separated by commas without additional spaces
      --num_hits                        Number of reported hits
      --min_charge                      Minimal precursor charge 
      --max_charge                      Maximal precursor charge 
      --skip_decoy_generation           Use a fasta database that already includes decoy sequences
      --comet_param_file                Parameter file for comet database search. This will overwrite all other parameters for the comet search
      
      
    Identification validation:
      --peptide_min_length              Minimum accepted peptide length
      --fdr_peptide_threshold           Threshold for FDR filtering
      
      
    Quantification options      
      --quantification_fdr              FDR threshold to accept peptides for quantification
      --quantification_min_prob         Specify a minimum probability cut off for quantification
      --experiment_design               Text-file containing 2 columns: first with raw file names and second with names for experimental conditions


          
    Other options:
      --outdir                          The output directory where the results will be saved
      --email                           Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                             Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                        The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                       The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Validate inputs
params.raws = params.raws ?: { log.error "No read data privided. Make sure you have used the '--raws' option."; exit 1 }()
params.fasta = params.fasta ?: { log.error "No fasta file privided. Make sure you have used the '--fasta' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()


/*
 * Define the default parameters
 */

//MS params
params.peptide_min_length = 8
//params.peptide_max_length = 12
params.fragment_mass_tolerance = 1.0005
params.precursor_mass_tolerance = 20
params.fragment_bin_offset = 0.4
params.fdr_peptide_threshold = 0.01
params.number_mods = 3

params.comet_param_file = "none"

params.num_hits = 100
params.pick_ms_levels = 2
params.run_centroidisation = false

params.min_charge = 2
params.max_charge = 3
params.activation_method = 'ALL'

params.enzyme = 'Trypsin'
params.miscleavages = 2
params.fixed_mods = 'Carbamidomethylation of C'
params.variable_mods = 'Oxidation of M'


params.skip_decoy_generation = false
if (params.skip_decoy_generation) {
  log.warn "Be aware: skipping decoy generation will prevent generating variants and subset FDR refinement"
  log.warn "Decoys have to be named with DECOY_ as prefix in your fasta database"
}

params.quantification_fdr = 0.01
params.quantification_min_prob = 0
params.min_num_peptides = 2

params.experiment_design = "none"

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


/*
 * Create a channel for input raw files
 */
input_raw = Channel.fromPath( params.raws )
input_raw.into { raws_convert; raws_merge_quant }

/*
 * Create a channel for fasta file
 */
input_fasta = Channel.fromPath( params.fasta )
input_fasta.into { fasta_search_comet; fasta_peptideprophet; fasta_stpeter; fasta_qc }

/*
 * Create a channel for comet parameter file
 */
input_comet_param = Channel.fromPath( params.comet_param_file )

if ((params.comet_param_file != "none") && !(file(params.comet_param_file).exists())) {
  log.error "Specified comet parameter file does not exit"; exit 1
}


/* 
 * Create a channel for experimental design file
 */
input_exp_design = Channel.fromPath(params.experiment_design)
if (params.experiment_design == "none") {
  log.warn "No experimental design! All raw files will be considered being from the one and the same experimental condition."
} else if(!(file(params.experiment_design).exists())) {
  log.error "File with experimental design does not exit"; exit 1  
}


/*
 * STEP 1 - convert raw files to mzml
 */
process convert_raw_mzml {
  label 'process_low'
  label 'process_single_thread'

  publishDir "${params.outdir}/mzMLs", mode:'copy'
  
  input:
  file rawfile from raws_convert
  
  output:
  file "${rawfile.baseName}.mzML" into (mzmls_search_comet, mzmls_peptideprophet)
  
  script:
  """
  ThermoRawFileParser.sh -i ${rawfile} -o ./  -f 2 -z
  """
}


/*
 * STEP 3 - set comet parameter file
 */
process set_comet_configuration {
  label 'process_very_low'
  label 'process_single_thread'

  publishDir "${params.outdir}/params", mode:'copy'
  
  input:
  file comet_param from input_comet_param.ifEmpty(file("none"))
  
  output:
  file "comet.params" into (comet_param_file)
  
  script: 
  if (comet_param.getName() == "none") {  
    // no param file provided
    ttrue = true
    enzymemap = ["Trypsin": 1, "Trypsin/P": 2, "Lys_C": 3, "Lys_N": 4, "Arg_C": 5, "Asp_N": 6, "CNBr": 7, "Glu_C": 8, "PepsinA": 9, "Chymotrypsin": 10, "Unspecified": 0]
    enzyme = enzymemap[params.enzyme]
    if (enzyme == null) {
      enzyme = enzymemap["Unspecified"]
    }
    skip_decoy = params.skip_decoy_generation
    modmap = ["Oxidation of M": "15.9949 M 0 3 -1 0 0 0.0", "Phosphorylation of STY": "79.966331 STY 0 3 -1 0 0 97.976896", "Acetylation of K": "42.010565 K 0 3 -1 0 0", "Acetylation of protein N-term": " 42.010565 n 0 3 0 0 0", "none": "0.0 X 0 3 -1 0"]
    mods = params.variable_mods.split(",")  
    modout = ""
    for (int i=1; i<9; i++) {
      if(mods.size() > i-1 ) {
        tmod =  modmap[mods[i-1]]
        if (tmod != null) {
          modout += "variable_mod0" + i + " = " + tmod + "\\n"
        } else {
          modout += "variable_mod0" + i + " = " + modmap["none"] + "\\n"
        }
      } else {
        modout += "variable_mod0" + i + " = " + modmap["none"] + "\\n"
      }
    }
    
    """
    comet -p
    mv comet.params.new comet.params
    if [ "${skip_decoy}" = "${true}" ]
    then
      sed -i 's/^decoy_search.*/decoy_search = 0/' comet.params
    else
      sed -i 's/^decoy_search.*/decoy_search = 1/' comet.params
    fi
    # given in ppm
    sed -i 's/^num_threads.*/'"num_threads = ${task.cpus}"/ comet.params
    sed -i 's/p^eptide_mass_tolerance.*/'"peptide_mass_tolerance = ${params.precursor_mass_tolerance}"/ comet.params
    sed -i 's/^search_enzyme_number.*/'"search_enzyme_number = ${enzyme}"/ comet.params
    # given in da and assuming high resolution
    sed -i 's/^fragment_bin_tol.*/'"fragment_bin_tol = ${params.fragment_mass_tolerance}"/ comet.params
    sed -i 's/^fragment_bin_offset.*/'"fragment_bin_offset = 0.0"/ comet.params
    # mass range
    sed -i 's/^activation_method.*/'"activation_method = ${params.activation_method}"/ comet.params
    sed -i 's/^allowed_missed_cleavage.*/'"allowed_missed_cleavage = ${params.miscleavages}"/  comet.params
    sed -i 's/^max_variable_mods_in_peptide.*/'"max_variable_mods_in_peptide = ${params.number_mods}"/ comet.params
    sed -i 's/^num_results.*/'"num_results = ${params.num_hits}"/ comet.params
    sed -i 's/^precursor_charge.*/'"precursor_charge = ${params.min_charge} ${params.max_charge}"/  comet.params  
    insert_line=\$(grep -n variable_mod0 comet.params | grep -Eo '^[^:]+' | head -n1)
    sed -i '/^variable_mod0/d' comet.params
    sed -i "\$insert_line i ${modout}" comet.params
    """ 
    // with a given comet file, all other paramter will be discarded
  } else {
    """
    cp ${comet_param} comet.params
    """
  }
}


/*
 * STEP 3 - run comet
*/ 
process db_search_comet {
  label 'process_medium'

  publishDir "${params.outdir}/comet", mode:'copy'
  
  input:
  each mzml_file from mzmls_search_comet
  file fasta from fasta_search_comet
  file comet_param from comet_param_file
  
  output:
  file "${mzml_file.baseName}.pep.xml" into pepxml
  
  script:
  """
  comet -P"${comet_param}" -N"${mzml_file.baseName}" -D"${fasta}" "${mzml_file}"
  """
}


/*
 * STEP 4 - run PeptideProphet to validate PSMs
*/ 
process run_peptideprophet {
  label 'process_medium'
  label 'process_single_thread'

  publishDir "${params.outdir}/peptideprophet", mode:'copy'
  
  input:
  each pepxml_file from pepxml
  file fasta from fasta_peptideprophet
  file mzml from mzmls_peptideprophet.collect()
  
  output:
  file "${pepxml_file.baseName}.interact.pep.xml" into interact_pepxml
  
  script:
  enzymemap = ["Trypsin": "", "Trypsin/P": "", "Lys_C": "-eN", "Lys_N": "-eL", "Arg_C": "-eN", "Asp_N": "-eA", "CNBr": "-eM", "Glu_C": "-eG", "PepsinA": "-eN", "Chymotrypsin": "-eC", "Unspecified": "-eN"]
  enzyme = enzymemap[params.enzyme]
  
  """
  ls ${mzml}
  xinteract -N"${pepxml_file.baseName}.interact.pep.xml" -p"${params.fdr_peptide_threshold}" ${enzyme} -l"${params.peptide_min_length}" -THREADS=${task.cpus} -PPM -O -D"${fasta}" "${pepxml_file}"
  """
}


/*
 * STEP 4 - run ProteinProphet to validate on protein level
*/ 
process run_proteinprophet {
  label 'process_medium'
  label 'process_single_thread'
  
  publishDir "${params.outdir}/proteinprophet", mode:'copy', pattern:'*.prot.xml'
  
  input:
  file pepxml_file from interact_pepxml
  
  output:
  tuple file("${pepxml_file.baseName}.prot.xml"), file(pepxml_file) into protxml
  
  script:
  """
  ProteinProphet "${pepxml_file}" "${pepxml_file.baseName}.prot.xml"
  """
}


/*
 * STEP 5 - run StPeter for protein quantification (label-free)
*/ 
process run_stpeter {
  label 'process_medium'
  label 'process_single_thread'
  
  publishDir "${params.outdir}/stpeter", mode:'copy', pattern:'*.prot.xml'
  
  input:
  tuple file(protxml), file(pepxml) from protxml 
  each file(fasta) from fasta_stpeter
  
  output:
  tuple file("${protxml.baseName}_stpeter.prot.xml"), file(pepxml) into stpeter_output
  
  script:
  """
  cp "${protxml}" stpeter_in.prot.xml
  StPeter -f ${params.quantification_fdr} -t ${params.fragment_mass_tolerance} stpeter_in.prot.xml
  # -i option not available anymore:    StPeter -i -f ${params.quantification_fdr}  "${protxml}"
  mv stpeter_in.prot.xml "${protxml.baseName}_stpeter.prot.xml"
  """
}


/*
 * STEP 6 - Convert StPeter prot.xml to csv
*/ 
process convert_stpeter {
  label 'process_low'
  label 'process_single_thread'
  
  publishDir "${params.outdir}/", mode:'copy'
  
  input:
  tuple file(stpeter), file(pepxml) from stpeter_output
  
  output:
  tuple file("${stpeter.baseName}_prot.csv"), file("${stpeter.baseName}_pep.csv") into protquant
  
  script:
  """
  cp "${stpeter}" StPeterOut.prot.xml 
  cp "${pepxml}" Sample.pep.xml
  Rscript $baseDir/scripts/proxml2csv.R --xml=StPeterOut.prot.xml --fdr=${params.quantification_fdr} --npep ${params.min_num_peptides}
  #    python $baseDir/scripts/protXML2csv.py 
  #    mv StPeterProts.csv "${stpeter.baseName}_prot.csv"
  #    mv StPeterPeps.csv "${stpeter.baseName}_pep.csv"
  mv StPeterOut.prot.xml.csv "${stpeter.baseName}_pep.csv"
  mv StPeterOut.prot.xml_summary.csv "${stpeter.baseName}_prot.csv"
  """
}


/*
 * STEP 7 - merge files according to experimental design
*/ 
process run_merge_quant {
  label 'process_medium'
  label 'process_single_thread'
  
  publishDir "${params.outdir}/", mode:'copy'
  
  input:
  file csv_files from protquant.collect()
  val raw_files from raws_merge_quant.collect()
  file exp_design_file from input_exp_design.ifEmpty(file("none"))
  
  output:
  file "all_prot_quant_merged.csv" into allprotquant
  file "all_pep_quant_merged.csv" into allpepquant
  file "stand_prot_quant_merged.csv" into stdprotquant
  file "stand_pep_quant_merged.csv" into stdpepquant
  file "exp_design.txt" into expdesign
  
  script:
  if (exp_design_file.getName() == "none") {
    // no experimental design file provided
    expdesign_text = "raw_file\texp_condition"
    for( int i=0; i<raw_files.size(); i++ ) {
      expdesign_text += "\n${raw_files[i].name}\tMain"
    }
    """
    touch exp_design.txt  
    echo "${expdesign_text}" >> exp_design.txt
    R CMD BATCH $baseDir/scripts/MergeOutput.R
    """
  } else {
    """
    cp "${exp_design_file}" exp_design.txt
    R CMD BATCH $baseDir/scripts/MergeOutput.R
    """
    }
}

/*
 * STEP 8 - Some QC
*/
process run_final_qc {
  label 'process_medium'
  label 'process_single_thread'
  
  publishDir "${params.outdir}/", mode:'copy'
  
    input:
        val foo from JsonOutput.prettyPrint(JsonOutput.toJson(params))
	file exp_design_file from expdesign
	file std_prot_file from stdprotquant
	file std_pep_file from stdpepquant
	file fasta_file from fasta_qc
 
  output:
   file "params.json" into parameters
   file "benchmarks.json" into benchmarks
  
  script:
  """
  echo '$foo' > params.json
  cp "${fasta_file}" database.fasta
  R CMD BATCH $baseDir/scripts/CalcBenchmarks.R

  """
}



workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}
