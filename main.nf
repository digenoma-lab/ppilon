#!/usr/bin/env nextflow


//help function for the tool
def show_help (){

  log.info tool_header()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --reads '*_R{1,2}.fastq.gz'  -profile singularity

    Mandatory arguments:
      --reads [file]                Path to input data

    References
      --fasta [file]                  Path to fasta reference


      -profile [str]              Configuration profile to use.
                                  Available: singularity

      --bwa                       default: bwa-mem2
                                  optional: bwa


      Test run:

      nextflow run  ppilon/main.nf --reads 'test_dataset/reads/*.R{1,2}.fastq.gz' --fasta test_dataset/genome.fa

      """.stripIndent()
}

// we star coding the pipeline
params.outdir="results"
params.nsplit=200
params.help = null
params.bwa="bwa-mem2"
params.cpu=20

// Show help message
if (params.help) exit 0, show_help()

//check star index or fasta and GTF
if (!params.fasta) exit 1, "Either specify an assembly fasta file!"

ch_fasta = Channel.value(file(params.fasta)).ifEmpty{exit 1, "Fasta file not found: ${params.fasta}"}



//print header and tool name
log.info tool_header()


//expect a file with header "label fwd_path rev_path"
//see file ./test_dataset/sample_fwrev.txt

    //expect a regular expresion like '*_{1,2}.fastq.gz'
    Channel.fromFilePairs(params.reads, size: 2 )
        .ifEmpty{exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
        .into{read_files_pilon; test_channel}


/*
 * Build STAR index
 */

process build_bwa_index {
    tag "bwa-index"
    memory '30 GB'
    publishDir params.outdir, mode: 'copy'

    input:
        file(fasta) from ch_fasta


    output:
        file("${fasta}*") into bwa_index

    //when: !(params.bwa_index)

    script:
     //adjust a variable for working with a smaller reference
    //def opt_test = params.test ? "--genomeSAindexNbases 8" : ''
    """

    ${params.bwa} index ${fasta}

    """
}

//ch_star_index = params.bwa_index ? Channel.value(file(params.star_index)).ifEmpty{exit 1, "STAR index not found: ${params.star_index}" } : star_index
ch_bwa_index = bwa_index
ch_bwa_index = ch_bwa_index.dump(tag:'ch_bwa_index')

//map the rna-seq reads to the genome

process bwa_mapping{
  tag "${sample}"
  cpus params.cpu
  memory '40 GB'

  //we can remove this to don't keep the bam files
  publishDir "${params.outdir}/bwa_mapping", mode: 'copy'

  input:
      set val(sample), file(reads) from read_files_pilon
      file(bwa_index) from ch_bwa_index
      file(fasta) from ch_fasta

  output:
      //star bam files
      set val(sample), file("${sample}_bwa.bam"), file("${sample}_bwa.bam.bai") into bwa_bam
      //set val(sample),  into bwa_bam_index
      //star mapping stats and gene counts *.{tsv,txt}
      file("${fasta}.fai") into ch_fai

  script:
  """
  ${params.bwa} mem -t ${params.cpu} ${fasta} ${reads} | samtools sort -@8 -o ${sample}_bwa.bam -
  samtools index ${sample}_bwa.bam
  samtools faidx ${fasta}
  """
  /*
  #The STAR gene counts coincide with those produced by htseq-count with default parameters.
  #The file colums are:
  #column 1: gene ID
  #column 2: counts for unstranded RNA-seq
  #column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
  #column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse
  */
}

process split_contigs {
    //tag "${fastaFai}"
    //tag "split_contigs"
    tag "${fastaFai}"
    publishDir params.outdir+"/parts", mode: 'copy'

    input:
        //file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        //file '*.targetlist' into split_parts
        file '*.targetlist' into split_parts mode flatten

    script:
    """
    perl ${baseDir}/split-faix.pl ${fastaFai} ${params.nsplit}

    """
}


//we combine the splitted part with bwa_mem
split_beds=bwa_bam.combine(split_parts)
//split_beds.view()


process pilon {
    tag "${sample}-${targetfiles}"
    memory '5 GB'

    publishDir params.outdir+"/pilon", mode: 'copy'

    input:
        set val(sample), file(frag),file(index), file(targetfiles) from split_beds
        file(fasta) from ch_fasta

    output:
       //file(targetfiles.basename+".polish.fa") into pilon_out
      file("${targetfiles.baseName}.pilon.fasta") into pilon_out


    //when: (!params.no_intervals) && step != 'annotate'

    script:
    bed_tag = targetfiles.baseName
    """
    #echo java -jar pilon-1.24.jar --genome ${fasta} --frags ${frag} --targets ${targetfiles} --output ${bed_tag}.pilon
    #touch ${bed_tag}.pilon.fa
    pilon --genome ${fasta} --frags ${frag} --targets ${targetfiles} --output ${bed_tag}.pilon
    """

}



//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        nextflow pipeline for polishing with pilon (${workflow.manifest.version})
        """
}
