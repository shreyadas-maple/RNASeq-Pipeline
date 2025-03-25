
// Defining the input parameters for the pipeline

// First we have the gut fastq files 1 and 2
params.reads = "$projectDir/data/ggal/*_{1,2}.fq"

// Then we have the transcriptome file which we will index
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"

// This is the MultiQC report that will be generated
params.multiqc = "$projectDir/multiqc"

// This is an output directory
params.outdir = "results"

log.info """\
    RNASEQ-NF PIPELINE
    ==================
    transcriptome: ${params.transcriptome_file}
    reads: ${params.reads}
    outdir: ${params.outdir} 
    """
    .stripIndent(true)


/*
 * This is the INDEX process that will generate a
 * binary index given a transcriptome file
 */
process INDEX {

    // Directive used to specify number of CPUs to execute this process
    cpus 2

    input:
    path transcriptome

    /* 
     * This is the using Salmon tool to quantify the transcripts
     * through RNA-seq data. This generates a salmon_index, which 
     * is an indexed transcriptome that will be later passed in the
     * workflow index_ch channel
     */

    output:
    path 'salmon_index'

    /* NOTE: task.cpus is a way that we can define the number of cpus
     * used to run this workflow for different platforms, without chaning
     * the entire workflow
     */
    script:
    """
    salmon index -- threads $task.cpus -t $transcriptome -i salmon_index
    """
}

/*
 * This process is important for expression quantification. It takes
 * the index trascriptome from the index_ch and the read pairs from 
 * read_pair_ch. The read_pair_ch is composed of a tuple of 2 items: 
 * sample_id, reads;this is the match the structure of the item emitted 
 * by the fromFilePairs.
 */
process QUANTIFICATION {
    // Adding a tag on the process
    tag "Salmon on $sample_id"
    publishDir params.outdir, mode: 'copy'

    input:
    path salmon_index
    // Remember that the output for the read_pair_ch is a tuple, so we have
    // to match that here
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """

}

process FASTQC {
    // Adding a tag on the process
    tag "FASTQC on $sample_id"

    // Remember this is the output from the read_pair_ch
    input:
    tuple val(sample_id), path(reads)

    //
    output:
    path "fastqc_${sample_id}_logs"

    /*
     * The script here says that we are making a directory
     */
    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

/*
 * We are going to create process that creates a MultiQC report.
 * Because we running many different processes on the same data, 
 * it would be easier to add all the outputs together into a
 * multiQC report.
 */

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {

    /*
     * .set is operator to name a channel. It is the same as saying 
     * read_pairs_ch = ...
     *
     * checkIfExists will check if the file that is being inputed in 
     * the channel exist in the directory
     */

    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set{read_pairs_ch}

    /* 
     * Here we are providing the input transcriptome file 
     * to the INDEX process, which goes via the index_ch channel
     */
    index_ch = INDEX(params.transcriptome_file)

    /*
     * Here
     */
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)

    fastqc_ch = FASTQC(read_pairs_ch)

    // Here we are mixing the quant channel and fastqc channel
    // and creating a single element with all the outputs
    MULTIQC(quant_ch.mix(fastqc_ch).collect())


    // This print a tuple made up of 2 items: the first is the read pair
    // prefix and the second is a list representing the actual pair of files
    //read_pairs_ch.view()

    
}

