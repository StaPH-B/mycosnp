#!/usr/bin/env nextflow

//Deshellion: Workflow that prepares a reference FASTA file for BWA alignment and 
//GATK variant calling by masking repeats in the reference and generating the BWA index.
//Author: Lynn Dotrang
//eMail: thoai.dotrang@dgs.virginia.gov

Channel
    .fromPath(  "/home/ldotrang/C_auris_testing/mycosnp-demo/test-data/reference/*.fasta" )
    .ifEmpty { exit 1, "Cannot find any reference fasta" }
    .into { reference_seq; reference_seq2; reference_seq3 }


// Mask repeats
process maskRepeats_1{
    input: 
    file (reference) from reference_seq

    output: 
    file("masked_ref.bed") into masked_ref_bed
    shell:
    // There is some stuff with awk here idk what it's doing exactly
    // Here is the basic idea of what should happen in this step
    """
    #!/bin/bash -ue
    nucmer --maxmatch --nosimplify -t ${task.cpus} ${reference} ${reference}
    show-coords -r -T -H out.delta > masked_ref_BEFORE_ORDER.bed
    awk '{if (\$1 != \$3 && \$2 != \$4) print \$0}' masked_ref_BEFORE_ORDER.bed > masked_ref_BEFORE_ORDER2.bed
    awk '{print \$8\"\\t\"\$1\"\\t\"\$2}' masked_ref_BEFORE_ORDER2.bed > masked_ref.bed
    """
}

process maskRepeats_2 {
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_reference_output", mode: 'copy'
    input: 
    file (reference) from reference_seq2
    file (bed) from masked_ref_bed

    output: 
    file("masked_reference_seq.fasta") into masked_reference_seq_dict, masked_reference_seq_idx
    shell:
    // separated this part because different tool
    """
   
    bedtools maskfasta -fi ${reference} -bed ${bed} -fo masked_reference_seq.fasta

    """
}
// BWA index
process bwaIndex {
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_reference_output/bwa_index", mode: 'copy'
    input:
    file (reference) from reference_seq3

// How do I know if I am handling a group of files correctly?
// bwa index will produce a set of 5 files
// I removed the prefix option for now (-p)
    output:
    file("*") into bwa_index

    shell:
    """
    bwa index ${reference}
    """
}
// Index reference 
// Here we are using one input to do multiple things
process indexReference_1{
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_reference_output", mode: 'copy'
    input:
    file (masked_ref) from masked_reference_seq_dict

    output:
    file ("sequence_dictionary.dict") into reference_sequence_dictionary

    shell:
    """
    picard CreateSequenceDictionary R= ${masked_ref}  O= sequence_dictionary.dict
    """
}

//separated this step for picard and samtools
process indexReference_2{
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/bwa_reference_output", mode: 'copy'
    //generates a .fasta and a .fasta.fai file 
    input:
    file (masked_ref) from masked_reference_seq_idx

    output:
    file("*") into indexed_reference_fasta

    shell:
    """
    samtools faidx ${masked_ref}
    """
}
