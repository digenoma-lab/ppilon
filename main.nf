
nextflow.enable.dsl = 2

def show_help() {

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

def tool_header() {
        return """
        nextflow pipeline for polishing with pilon (${params.version})
        """
}

process build_bwa_index {
    tag "bwa-index"
    publishDir "$params.outdir/bwa_index", mode: "copy"

   input:
      file(fasta)

   output:
     path "${fasta}*", emit: bwa_index

  script:
    //def prefix=${fasta}.baseName()

    if(params.debug){
      """
      echo ${params.bwa} index ${fasta}
      touch ${fasta}.0123 ${fasta}.amb ${fasta}.ann  ${fasta}.pac
      """
    }else{
      """
      ${params.bwa} index ${fasta}
      """
    }
}

//ch_bwa_index = bwa_index
//ch_bwa_index = ch_bwa_index.dump(tag:'ch_bwa_index')

process bwa_mapping {
  tag "${sample}"

  //we can remove this to don't keep the bam files
  publishDir "${params.outdir}/bwa_mapping", mode: 'copy'

  input:
      tuple val(sample), file(reads)
      file(bwa_index)
      file(fasta)

  output:
      //star bam files
      tuple val(sample), file("${sample}_bwa.bam"), file("${sample}_bwa.bam.bai"), emit: bam
      path("${fasta}.fai"), emit: fai

  script:
  if(params.debug){
  """
  echo ${params.bwa} mem -t ${params.cpu} ${fasta} ${reads}  samtools sort -@8 -o ${sample}_bwa.bam -
  echo samtools index ${sample}_bwa.bam
  echo samtools faidx ${fasta}
  touch ${sample}_bwa.bam ${sample}_bwa.bam.bai ${fasta}.fai
  """
  }else{
  """
  ${params.bwa} mem -t ${params.cpu} ${fasta} ${reads} | samtools sort -@8 -o ${sample}_bwa.bam -
  samtools index ${sample}_bwa.bam
  samtools faidx ${fasta}
  """
  }
}

process split_contigs {
    tag "${fastaFai}"
    publishDir params.outdir+"/parts", mode: 'copy'

    input:
        file(fastaFai)

    output:
        path '*.targetlist' , emit: split_parts

    script:
    if (params.debug){
    """
    echo perl ${baseDir}/split-faix.pl ${fastaFai} ${params.nsplit}
    touch pp1.targetlist pp2.targetlist pp3.targetlist
    """
    }else{
    """
    perl ${baseDir}/split-faix.pl ${fastaFai} ${params.nsplit}
    """
    }
}


//we combine the splitted part with bwa_mem
//split_beds=bwa_bam.combine(split_parts)
//split_beds.view()


process pilon {
    tag "${targetfiles}"
    //memory '10 GB'

    publishDir params.outdir+"/pilon", mode: 'copy'

    input:
        tuple val(sample), file(frag), file(index), file(targetfiles), file(fasta)
        //file(fasta)
    output:
        tuple val("${sample}-${targetfiles}"), path("${targetfiles.baseName}.pilon.fasta"), emit: fasta

    script:
    def bed_tag = targetfiles.baseName
    if(params.debug){
    """
    echo pilon --genome ${fasta} --frags ${frag} --targets ${targetfiles} --output ${bed_tag}.pilon
    touch ${bed_tag}.pilon.fasta
    """
    }else{
    """
    pilon --genome ${fasta} --frags ${frag} --targets ${targetfiles} --output ${bed_tag}.pilon
    """
    }
}


//run pilon
workflow {
// Show help message
if (params.help) exit 0, show_help()
//check star index or fasta and GTF
if (!params.fasta) exit 1, "Either specify an assembly fasta file!"
ch_fasta = Channel.value(file(params.fasta)).ifEmpty{exit 1, "Fasta file not found: ${params.fasta}"}
//print header and tool name
log.info tool_header()

//expect a regular expresion like '*_{1,2}.fastq.gz'
sreads=Channel.fromFilePairs(params.reads, size: 2 )
    .ifEmpty{exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
sreads.view()

//we index the assembly file
build_bwa_index(ch_fasta)
//step 2: map the short-reads
bwa_mapping(sreads,build_bwa_index.out.bwa_index,ch_fasta)
//step 3 : split contigs
split_contigs(bwa_mapping.out.fai)
//step 4: run pilon
pilon_in=bwa_mapping.out.bam.combine(split_contigs.out.split_parts.flatten())
//pilon_in=pilon_in.combine(bwa_mapping.out.fai)
pilon_in=pilon_in.combine(ch_fasta)
//pilon_in.view()
pilon(pilon_in)
}
