#!/usr/bin/env nextflow

//Deshellion: Workflow that prepares a reference FASTA file for BWA alignment and 
//GATK variant calling by masking repeats in the reference and generating the BWA index.
//Author: Lynn Dotrang
//eMail: thoai.dotrang@dgs.virginia.gov

bam_ch = Channel
    .fromPath("/home/ldotrang/C_auris_testing/mycosnp-demo/mycosnp-gatk-variants/data/bam_index/all_bams/*.bam")
    .ifEmpty { exit 1, "Cannot find bam files" }
    .map {file -> tuple(file.baseName, file)}
// Channel
//     .fromPath("/home/ldotrang/C_auris_testing/mycosnp-demo/mycosnp-gatk-variants/data/bam_index/SRR13710812/SRR13710812.bam.bai")
//     .ifEmpty { exit 1, "Cannot find bam files" }
//     .set {bam_index}

Channel
    .fromPath(  "/home/ldotrang/C_auris_testing/mycosnp-demo/mycosnp-gatk-variants/data/indexed_reference/indexed_reference.*" )
    .ifEmpty { exit 1, "Cannot find any reference fasta" }
    .into {reference_fasta; reference_fasta2}

// Channel
//     .fromPath(  "/home/ldotrang/C_auris_testing/mycosnp-demo/mycosnp-gatk-variants/data/indexed_reference/indexed_reference.fasta.fai" )
//     .ifEmpty { exit 1, "Cannot find any reference fasta" }
//     .set {reference_index}

// Channel
//     .fromPath(  "/home/ldotrang/C_auris_testing/mycosnp-demo/mycosnp-gatk-variants/data/indexed_reference/indexed_reference.fasta.dict" )
//     .ifEmpty { exit 1, "Cannot find any reference fasta" }
//     .set {reference_dict}



process gatk_HC{
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/gatk_output/gvcf", mode: 'copy'
    tag "$name"

    //There are some default params here:
    // sample ploidy: 2
    // ref confidence mode (NONE, BP_RESOLUTION, GVCF): 'NONE'
    // threads: 4
    input:
   
    file reference from reference_fasta.collect()
    // file ref_index from reference_index.view()
    // file ref_dict from reference_dict.view()
        tuple name, file(bams) from bam_ch
    

    output:
    tuple name, file("${name}_gatkHC.g.vcf") into gvcfs
    tuple name, file("${name}_gatkHC.g.vcf.idx") into gvcf_idx
  
    shell:
    """
    gatk HaplotypeCaller -R ${reference[1]} -I ${bams} --sample-ploidy 2 --emit-ref-confidence GVCF --native-pair-hmm-threads 4 \
    -O ${name}_gatkHC.g.vcf
    """
}

process gatk_CB{
    publishDir "/home/ldotrang/C_auris_testing/nfBWAref/gatk_output/gvcf_combined", mode: 'copy'
    // removed the tuple and name tag and added a chunk of script to make a list of variants 

    input:
    file reference from reference_fasta2.collect()
    file (gvcf) from gvcfs.collect()

    output:
    file("*.g.vcf") into combined_gvcfs
    file("genotyped.vcf") into genotyped_vcf
    file("variant_filter.vcf") into filtered_vcf
    file("selected.vcf") into selected_vcf
    // file("selected-clean.vcf") into selected_clean_vcf
  

    script:
    //vcfs = gvcf.join(' --variant ') this did not end up working 
    // filter expression is default from MycoSNP
    // select type to include default is SNP
    // there is a script step if the type selected is SNP, we are going to assume it is
    // script not added yet!
    """
    for file in ${gvcf}; do
            if [ \${file##*.} = "vcf" ]; then
                echo "-V \${file}";
            fi
        done > gvcf.list

   
    gatk CombineGVCFs --output combined.g.vcf \
    --reference ${reference[1]}  --arguments_file gvcf.list \

    gatk GenotypeGVCFs --output genotyped.vcf --reference ${reference[1]} --variant combined.g.vcf

    gatk VariantFiltration --output variant_filter.vcf --reference ${reference[1]} --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "filter" --variant genotyped.vcf

    gatk SelectVariants --output selected.vcf --reference ${reference[1]} \
    --select-type-to-include "SNP" --variant variant_filter.vcf 

    
    """
}
// awk '\$5 !~\/[GATC],\*\/' selected.vcf | awk '\$5 !~\/\*\/' >> selected-clean.vcf
