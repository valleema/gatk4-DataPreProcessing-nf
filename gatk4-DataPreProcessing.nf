#! /usr/bin/env nextflow

// Copyright (C) 2018 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "    gatk4-DataPreProcessing-nf: From FASTQ to hg38 Callable BAM and gVCF      "
log.info "-------------------------------------------------------------------------"
log.info ""

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/gatk4-DataPreProcessing-nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input                  FASTQ FILE            FASTQ files (between quotes)"
    log.info "--output_dir             OUTPUT FOLDER         Output for gVCF file"
    log.info "--ref_fasta              FASTA FILE            Reference FASTA file"
    log.info "--dbsnp                  VCF FILE              dbSNP VCF file"
    log.info "--onekg                  VCF FILE              1000 Genomes High Confidence SNV VCF file"
    log.info "--mills                  VCF FILE              Mills and 1000 Genomes Gold Standard SID VCF file"
    log.info "--interval_list          INTERVAL_LIST FILE    GATK wgs calling regions interval_list file"
    log.info "--bait_bed               BED FILE              Regions of capture probes"
    log.info "--target_bed             BED FILE              Regions targeted"
    exit 1
}


//
// Parameters Init
//
params.input         = null
params.output_dir    = "."
params.ref_fasta     = null
params.dbsnp         = null
params.onekg         = null
params.mills         = null
params.interval_list = null
params.bait_bed      = null
params.target_bed    = null


//
// Parse Input Parameters
//
Channel.fromFilePairs( params.input )
        .ifEmpty{exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!" }
        .set{fastq_ch}
output     = file(params.output_dir)
ref        = file(params.ref_fasta)
interList  = file(params.interval_list)
bait_bed   = file(params.bait_bed)
target_bed = file(params.target_bed)
ref_dbsnp  = file(params.dbsnp)
ref_1kg    = file(params.onekg)
ref_mills  = file(params.mills)
ref_alt    = ref.parent / ref.baseName + ".fa.alt"
ref_amb    = ref.parent / ref.baseName + ".fa.amb"
ref_ann    = ref.parent / ref.baseName + ".fa.ann"
ref_bwt    = ref.parent / ref.baseName + ".fa.bwt"
ref_pac    = ref.parent / ref.baseName + ".fa.pac"
ref_sa     = ref.parent / ref.baseName + ".fa.sa"

//
// Prep BEDs, intervals, faidx, dict
//
process Prep_Files{
  cpus 1
  time '4h'

  input:
    file genome from ref

  output:
    file "${genome}.fai" into faidx_ch1,faidx_ch2,faidx_ch3
    file "${genome.baseName}.dict" into dict_ch1,dict_ch2,dict_ch3
    set "bait.interval_list", "target.interval_list" into metrics_intList_ch
  
  script:
  """
    samtools faidx ${genome}

    java -jar /opt/picard.jar \
      CreateSequenceDictionary \
      R=${genome} \
      O=${genome.baseName}.dict

    java -jar /opt/picard.jar \
      BedToIntervalList \
        I=${bait_bed} \
        O=bait.interval_list \
        SD=${genome.baseName}.dict

    java -jar /opt/picard.jar \
      BedToIntervalList \
        I=${target_bed} \
        O=target.interval_list \
        SD=${genome.baseName}.dict
  """
}

//
// Process Mapping With BWA MEM (AND MORE)
//
process BWA_mapping{
  tag "$sampleID"

  cpus 16
  // memory '72 GB'
  time '12h'
  
  input: 
    file genome from ref
    file "${genome.baseName}.fa.alt" from ref_alt
    file "${genome.baseName}.fa.amb" from ref_amb
    file "${genome.baseName}.fa.ann" from ref_ann
    file "${genome.baseName}.fa.bwt" from ref_bwt
    file "${genome.baseName}.fa.pac" from ref_pac
    file "${genome.baseName}.fa.sa" from ref_sa
	  set sampleID, file(fastqs) from fastq_ch
  
  output: 
    set sampleID, file("${sampleID}.aln.bam"), file("${sampleID}.aln.bam.bai") into aln_bam_ch

  script:    
  """
    run-bwamem -s -t ${task.cpus} -R '@RG\\tID:${sampleID}\\tSM:${sampleID}\\tPL:Illumina' -o ${sampleID} $genome $fastqs | sh
    
    sambamba index ${sampleID}.aln.bam
  """

}


//
// Process Mark Duplicates
//
process SAMBAMBA_markdup{
  tag "$sampleID"

  cpus 8
  // memory '64 GB'
  time '12h'
  
  input: 
	  set sampleID, file(bam), file(bai) from aln_bam_ch
  
  output: 
      set sampleID, file("${sampleID}.aln.dup.bam"), file("${sampleID}.aln.dup.bam.bai") into dup_bam_ch

  script:    
  """
	sambamba markdup -t ${task.cpus} --tmpdir=$TMPnf/$sampleID $bam ${sampleID}.aln.dup.bam
	
  """

}


//
// Process BaseRecalibrator + ApplyBQSR
//
process GATK_BaseRecalibrator_ApplyBQSR{
  tag "$sampleID"

  cpus 1
  // memory '32 GB'
  time '4h'
    
  publishDir "$output/BAM/$sampleID", mode: 'copy'
  
  input: 
      file genome from ref 
      set sampleID, file(bam), file(index) from dup_bam_ch

  output: 
      set sampleID, file("${sampleID}.recal.bam"), file("${sampleID}.recal.bai") into recal_bam_ch1,recal_bam_ch2

  script:    
  """
    java -jar /opt/picard.jar \
    CreateSequenceDictionary \
    R=${genome} \
    O=${genome.baseName}.dict

    samtools faidx ${genome}
    
    mkdir -p $TMPnf/$sampleID

    java -Xmx4g -Xms4g -Djava.io.tmpdir=$TMPnf/$sampleID -jar /opt/gatk4.jar \
    		BaseRecalibrator \
    		-R ${genome} \
    		-I ${bam} \
        --use-original-qualities \
        -O ${sampleID}.recal_data.table \
        -L ${interList} \
        -known-sites ${ref_dbsnp} \
        -known-sites ${ref_1kg} \
        -known-sites ${ref_mills}
			
    java -Xmx4g -Xms4g -Djava.io.tmpdir=$TMPnf/$sampleID -jar /opt/gatk4.jar \
    		ApplyBQSR \
    		-R ${genome} \
    		-I ${bam} \
        --use-original-qualities \
        -O ${sampleID}.recal.bam \
        --bqsr-recal-file ${sampleID}.recal_data.table
			
  """

}


//
// Process Post-Alignment QC
//
process BamQC{
  tag "$sampleID"

  cpus 8
  // memory '32 GB'
  time '6h'
      
  publishDir "$output/QC/$sampleID", mode: 'copy'
  
  input: 
	  set sampleID, file(bam), file(bai) from recal_bam_ch1
	  set file(bait), file(target) from metrics_intList_ch
    file genome from ref
    file faidx from faidx_ch1
    file dict from dict_ch1

  output: 
    set sampleID, file("${sampleID}.hs_metrics.txt"), file("${sampleID}.stats.txt") into recal_bam_stats_ch

  script:    
  """
    samtools flagstat $bam > ${sampleID}.stats.txt
    
    java -jar /opt/picard.jar \
      CollectHsMetrics \
        I=$bam \
        O=${sampleID}.hs_metrics.txt \
        R=${genome} \
        BAIT_INTERVALS=${bait} \
        TARGET_INTERVALS=${target}

  """
}

//
// Process Split Intervals, to scatter the load
//
process SplitIntervals {
	cpus 1
	// memory '4 GB' 
	time '1h'

  input:
    file genome from ref
	  interList
    file faidx from faidx_ch2
    file dict from dict_ch2

	output:
    file "scatter/*-scattered.interval_list" into interval_ch

	script:
	"""
    java -Xmx4g -Xms4g -jar /opt/gatk4.jar \
      SplitIntervals \
      -R ${genome} \
      -L ${interList} \
      --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
      --scatter-count 10 \
      -O scatter
	
	"""
}


//
// Process launching HC
//
process HaplotypeCaller {

	cpus 4 // --native-pair-hmm-threads GATK HC argument is set to 4 by default
	// memory '64 GB'
	time '12h'

	tag { bamID+"-"+int_tag }
	
  input:
    file genome from ref
    file faidx from faidx_ch3
    file dict from dict_ch3
	  set bamID, file(bam), file(bai), file(Interval) from recal_bam_ch2.spread(interval_ch)

	output:
    set bamID, file("${bamID}.${int_tag}.g.vcf") , file("${bamID}.${int_tag}.g.vcf.idx") into gvcf_ch
	
  script:
	  int_tag = Interval.baseName

	"""
    java -Xmx8g -Xms8g -jar /opt/gatk4.jar \
      HaplotypeCaller \
      -R ${genome} \
      -I ${bam} \
      -O ${bamID}.${int_tag}.g.vcf \
      -L ${Interval} \
      -contamination 0 \
      -ERC GVCF		
    """
}	


//
// Process Merge and Sort gVCF
//
process MergeGVCFs {

	cpus 4
	// memory '24 GB'
	time '6h'

	tag { bamID }
	
  publishDir "$output/gVCF/$bamID", mode: 'copy'

  input:
    set bamID, file (gvcfs), file (gvcfidxs) from gvcf_ch.groupTuple()

	output:
    set bamID, file("${bamID}.g.vcf") , file("${bamID}.g.vcf.idx") into merged_gvcf_ch
	
    script:
	  """
      java -Xmx24g -Xms24g -jar /opt/gatk4.jar \
        SortVcf \
        ${gvcfs.collect { "--INPUT $it " }.join()} \
        --OUTPUT ${bamID}.g.vcf
    """
}





