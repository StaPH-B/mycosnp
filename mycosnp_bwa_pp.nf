#!/usr/bin/env nextflow

//Deshellion: Workflow that prepares a reference FASTA file for BWA alignment and 
//GATK variant calling by masking repeats in the reference and generating the BWA index.
//Author: Lynn Dotrang
//eMail: thoai.dotrang@dgs.virginia.gov


// Before completetion, need to edit for all full paths and all hard coded values (CPU usage options, etc)
Channel
    .fromFilePairs(  "/home/ldotrang/C_auris_testing/SampleTestData/CombinedPairsTest/SRR*{1,2}.fastq.gz", size: 2)
    .ifEmpty { exit 1, "Cannot find any reads matching" }
    .set {raw_reads}

Channel
    .fromPath(  "/home/ldotrang/C_auris_testing/mycosnp-demo/test-data/reference/*.fasta" )
    .ifEmpty { exit 1, "Cannot find any reference fasta" }
    .view()
    .into {reference_seq; reference_seq2; reference_seq3}

// Am I going to need this if I expect the correct file time every time?
//process concat_fastq_lanes{
//     input: 
//     files

//     output: 
//     file("masked_ref.bed") into masked_ref_bed
//     shell:
//     // Very long script goes here
    
//     """
//     #!/bin/bash -ue
//     ##SCRIPT THAT DOES THE CONCAT STUFF HERE
//     """
// }

process seqkit_fastq_pair{
    input: 
    tuple val(name), file(reads) from raw_reads

    output: 
    tuple name, file("${name}{_1,_2}.paired.fastq.gz") into filtered_reads
    shell:
    
    """
    seqkit pair -1 ${reads[0]} -2 ${reads[1]} -u -j 4
    """
}

process seqtk_downsample{
    
    tag "$name"
    input:
    tuple val(name), file(reads), file(reference) from filtered_reads.combine(reference_seq)
    

    output:
    tuple name, file("${name}{_1,_2}.fastq") into downsampled_reads, bwa_downsampled_reads, fastqc_downsampled_reads
    shell:
    """
    #!/bin/bash -ue

    #Set some variables first 
    COVERAGE=30
    SEED=\${RANDOM}
  

    REFERENCE_LEN=\$(awk '!/^>/ {len+=length(\$0)} END {print len}' < ${reference}) 
    echo \${REFERENCE_LEN} > seqtktestoutput.txt
    READS_LEN=\$(zcat ${reads[0]} ${reads[1]} | awk '/^@/ {getline; len+=length(\$0)} END {print len}') 
    echo \${READS_LEN} >> seqtktestoutput.txt
    SAMPLE_RATE=\$(echo "\${COVERAGE} \${READS_LEN} \${REFERENCE_LEN}" | awk '{x=\$1/(\$2/\$3); x=(1<x?1:x)} END {print x}')
    echo \${SAMPLE_RATE} >> seqtktestoutput.txt
    NUM_READS=\$(zcat ${reads[0]} | awk '/^@/ {lines+=1} END {print lines}')
    echo \${NUM_READS} >> seqtktestoutput.txt
    SAMPLED_NUM_READS=\$(echo "\${NUM_READS} \${SAMPLE_RATE}" | awk '{x=\$1*\$2} END {printf "%.0f", x}')
    echo \${SAMPLED_NUM_READS} >> seqtktestoutput.txt
   
    #need to run seqtk twice, once per read
    seqtk sample -s \${SEED} ${reads[0]} \${SAMPLED_NUM_READS} > ${name}_1.fastq
    seqtk sample -s \${SEED} ${reads[1]} \${SAMPLED_NUM_READS} > ${name}_2.fastq
    """
}
process faqcs{
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_preprocess_output/qc_reports", mode: 'copy'
    input: 
    set val(name), file(reads) from downsampled_reads

    output:
    file('*.pdf') into qc_printouts, mqc_faqcs

    shell:
    """
    FaQCs --prefix ${name} -1 ${reads[0]} -2 ${reads[1]} -d .
    """
}

process qc_report{
    // this part tuple does not work correctly, but it will still proceed ok 
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_preprocess_output/qc_reports", mode: 'copy'
    tag "$name"
    input:
    tuple name, file(report), file(reference) from qc_printouts.combine(reference_seq2)
    output:
    
    file("${name.baseName}.formattedreport.csv") into qc_report
  
    shell:
    //this villain again 
    
    """

    REFERENCE_LEN=\$(awk '!/^>/ {len+=length(\$0)} END {print len}' < ${reference})
    printf \"Sample Name\\t# Reads Before Trimming\\tGC Before Trimming\\tAverage Phred Before Trimming\\tCoverage Before Trimming\\t# Reads After Trimming\\t# Paired Reads After Trimming\\t# Unpaired Reads After Trimming\\tGC After Trimming\\tAverage Phred After Trimming\\tCoverage After Trimming\\n\" > ${name.baseName}.formattedreport.csv;
    pdftotext -raw \"${name}\" \"${name}.txt\"|| continue;
    printf \"${name.baseName}\" | awk -F'_qc_report' '{printf \$1\"\\t\"}' >> ${name.baseName}.formattedreport.csv
          cat \"${name}.txt\" | grep 'Reads #:' | sed -z 's/\\n/|/g' | sed 's/Reads #: //g' | sed 's/Paired //g' | sed 's/Unpaired //g' | awk -F'|' '{printf \$1\"\\t\"}'  >> ${name.baseName}.formattedreport.csv
          cat \"${name}.txt\" | grep 'GC' | grep -v 'Reads' | sed -z 's/\\n/|/g' | sed 's/GC //g' | sed 's/ ± /|/g' | awk -F'|' '{printf \$1\"\\t\"}'  >> ${name.baseName}.formattedreport.csv
          cat \"${name}.txt\" | grep 'Average' | grep -v 'Reads' | sed -z 's/\\n/|/g' | sed 's/Average: //g' | awk -F'|' '{printf \$1\"\\t\"}' >> ${name.baseName}.formattedreport.csv
          cat \"${name}.txt\" | grep 'Total bases:' | sed -z 's/\\n/|/g' | sed 's/Total bases: //g' | awk -F'|' '{x=\$1/\$REFERENCE_LEN} END {printf x\"\\t\"}' >> ${name.baseName}.formattedreport.csv
          cat \"${name}.txt\" | grep 'Reads #:' | sed -z 's/\\n/|/g' | sed 's/Reads #: //g' | sed 's/Paired //g' | sed 's/Unpaired //g' | awk -F'|' '{printf \$2\"\\t\"\$3\"\\t\"\$4\"\\t\"}' >> ${name.baseName}.formattedreport.csv
          cat \"${name}.txt\" | grep 'GC' | grep -v 'Reads' | sed -z 's/\\n/|/g' | sed 's/GC //g' | sed 's/ ± /|/g' | awk -F'|' '{printf \$3\"\\t\"}' >> ${name.baseName}.formattedreport.csv
          cat \"${name}.txt\" | grep 'Average' | grep -v 'Reads' | sed -z 's/\\n/|/g' | sed 's/Average: //g' | awk -F'|' '{printf \$2\"\\t\"}'  >> ${name.baseName}.formattedreport.csv
          cat \"${name}.txt\" | grep 'Total bases:' | sed -z 's/\\n/|/g' | sed 's/Total bases: //g' | sed 's/ /|/g' | awk -F'|' '{x=\$2/\$REFERENCE_LEN} END {printf x\"\\n\"}'  >> ${name.baseName}.formattedreport.csv
          

    """
}

process bwa_align{
    
    tag "$name"
    input:
    tuple name, file(reads), file(reference) from bwa_downsampled_reads.combine(reference_seq3)

    output:
    tuple name, file("${name}.sam") into bwa_sams
  
    shell:
    """

    bwa index ${reference}
    bwa mem -t 4 ${reference} ${reads[0]} ${reads[1]} > ${name}.sam

    """
}

process bam_sort{
    tag "$name"
    //maybe move the samtools step above to this process
    //default is to sort by coordinates
    //-n option is added to sort by query number
    input:
    tuple name, file(sam) from bwa_sams
    output:
    tuple name, file("${name}_sorted.bam") into sorted_bams
    shell:
    """
    samtools view -b ${name}.sam -o ${name}.bam
    samtools sort ${name}.bam > ${name}_sorted.bam
    """
}

process picard{
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_preprocess_output", mode: 'copy'
    tag "$name"
   
    //removed duplicates is default to true
    //assume_sort_order is default to coordinate
    //validation_stringency is default to lenient in og wf, but default is strict
    //original wf has all of these as changeable params

    //for addorreplacerreads there are a lot of params
    // library
    // barcode
    // platform
    // and MORE
    // VALIDATION_STRINGENCY="LENIENT"
    // REMOVE_DUPLICATES="false"
    // ASSUME_SORT_ORDER="coordinate"
    // EXEC_METHOD="auto"
    // EXEC_INIT=":"
    input:
    tuple name, file(bam) from sorted_bams
    

    output:
    tuple name, file("${name}_picard.bam") into picard_bams, qualimap_bams
    file ("*-metrics.txt") into mqc_picardmd
    shell:
    """
    picard MarkDuplicates I=${name}_sorted.bam O=${name}_marked_dup.bam METRICS_FILE=${name}-metrics.txt \
    VALIDATION_STRINGENCY="LENIENT" \
    ASSUME_SORT_ORDER="coordinate"

    picard CleanSam I= ${name}_marked_dup.bam O= ${name}_cleansam.bam VALIDATION_STRINGENCY=LENIENT

    picard FixMateInformation I= ${name}_cleansam.bam O= ${name}_fixmate.bam

    picard AddOrReplaceReadGroups I= ${name}_fixmate.bam O= ${name}_picard.bam ID=1 LB=lib1 PL=ILLUMINA PU=unit1 SM=sample
    """
}
process bam_index{
    
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_preprocess_output", mode: 'copy'
    tag "$name"
    input:
    tuple name, file (bam) from picard_bams
    output:
    tuple name, file("${name}_picard.bam.bai") into indexed_bam
    shell:
    """
    samtools index ${name}_picard.bam
    """
}


// additional QC steps

process fastqc{
    //same as hickory

    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_preprocess_output/qc_reports", mode: 'copy'
    tag "$name"

    input:
    tuple name, file(reads) from fastqc_downsampled_reads

    output: 
    file("*_fastqc.{zip,html}") into fastqc_results

    script:
    """
    fastqc -q  ${reads}
    """

}

process qualimap{
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_preprocess_output/qc_reports", mode: 'copy'
    tag "$name"

    input: 
    tuple name, file (bam) from qualimap_bams

    output:
    tuple name, file ("*") into qualimap_results

    script:
    """
    qualimap bamqc -bam ${bam}
    """

}

process multiqc{
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_preprocess_output/qc_reports", mode: 'copy'
    input:
    file faqcs from mqc_faqcs
    file picardmd from mqc_picardmd
    file fastqc from fastqc_results.collect()
    file qualimap from qualimap_results

    output:
    file("*multiqc_report.html") into multiqc_report
    file("*_data")

    script:
    """
    multiqc . >/dev/null 2>&1
    """
}
