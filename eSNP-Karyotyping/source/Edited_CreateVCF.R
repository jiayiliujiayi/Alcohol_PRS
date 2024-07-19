#' Jiayi rewrite Jun 2024
#' CreateVCF
#'
#'
#' Create the VCF (Variant Call Format) file
#' @param bamdir The Path of to the BAM Directory
#' @param Genome_Fa Path for whole genome FASTQ file
#' @param DIR_Genome_DICT
#' @param Picard_Path Path to the Picard directory, with all the JAR files
#' @param CMD_gatk
#' @param CMD_samtools
#' @param base_name
#' @param out_dir
#' @export
#' @return None


CreateVCF <- function(input_bam_dir = "./temp_dir",
                      Genome_Fa = Genome_Fa,
                      DIR_Genome_DICT = DIR_Genome_DICT, 
                      Picard_Path = Picard_Path,
                      CMD_gatk = CMD_gatk,
                      CMD_samtools = CMD_samtools, 
                      base_name,
                      out_dir = "./temp_dir") {
    
  ### CVCF1: AddOrReplaceReadGroups
    #print("AddOrReplaceReadGroups")
  CMD_AddOrReplaceReadGroups <- sprintf(
    "java -jar %s/picard.jar AddOrReplaceReadGroups I=%s/reheaded_%s.bam O=%s/%saccepted_hits_rg.bam ID=%s LB=lb PL=ILLUMINA PU=pu SM=%s",
    Picard_Path, input_bam_dir, base_name, out_dir, base_name, base_name, base_name
  )
  system(CMD_AddOrReplaceReadGroups)


  ### CVCF2: ReorderSam
    #print("ReorderSam")
  CMD_ReorderSam <- sprintf(
    "java -jar %s/picard.jar ReorderSam INPUT=%s/%saccepted_hits_rg.bam OUTPUT=%s/%saccepted_hits_rg_sorted.bam SEQUENCE_DICTIONARY=%s",
    Picard_Path, input_bam_dir, base_name, input_bam_dir, base_name, DIR_Genome_DICT
  )
  system(CMD_ReorderSam)

  ### CVCF3: BuildBamIndex ### THIS IS ESSENTIALLY THE SAME AS "Samtools index", but significantly slower because not supporting multithreading. 
    #print("BuildBamIndex")
#  CMD_BuildBamIndex <- sprintf(
#    "java -jar %s/picard.jar BuildBamIndex #I=%s/%saccepted_hits_rg_sorted.bam #O=%s/%saccepted_hits_rg_sorted.bai",
#    Picard_Path, input_bam_dir, base_name, input_bam_dir, base_name
#  )
#  system(CMD_BuildBamIndex)

  ### CVCF4: Samtools Index
    #print("Samtools Index")
  CMD_samtools_index <- sprintf(
    "%s index -@ 20 %s/%saccepted_hits_rg_sorted.bam -o %s/%saccepted_hits_rg_sorted.bai",
    CMD_samtools, input_bam_dir, base_name, input_bam_dir, base_name
  )
  system(CMD_samtools_index)
  #setwd("../")

  ### CVCF5: SplitNCigarReads
    #print("SplitNCigarReads")
  CMD_SplitNCigarReads <- sprintf(
    "%s SplitNCigarReads -R %s -I %s/%saccepted_hits_rg_sorted.bam -O %s/%s_split.bam",
    CMD_gatk, Genome_Fa, input_bam_dir, base_name, out_dir, base_name
  )
  system(CMD_SplitNCigarReads)

  ### CVCF6: HaplotypeCaller
    #print("HaplotypeCaller")
  CMD_HaplotypeCallerSpark <- sprintf(
    "%s HaplotypeCallerSpark -R %s -I %s/%s_split.bam --dont-use-soft-clipped-bases true -stand-call-conf 10.0 --native-pair-hmm-threads 24 -O %s/%s.vcf --spark-master local[20]",
    CMD_gatk, Genome_Fa, input_bam_dir, base_name, out_dir, base_name
  )
  system(CMD_HaplotypeCallerSpark)
    #print("done")

  # Example of additional command that might be used (commented out for now)
  # CMD_HaplotypeCaller <- sprintf(
  #   "%s HaplotypeCaller -R %s -I %s/%s_split.bam --dont-use-soft-clipped-bases true -stand-call-conf 10.0 --native-pair-hmm-threads 24 -O %s/%s.vcf",
  #   CMD_gatk, Genome_Fa, input_bam_dir, base_name, out_dir, base_name
  # )
  # system(CMD_HaplotypeCaller)

  # Cleanup (commented out for safety; uncomment if cleanup is intended)
  # unlink(x = input_bam_dir, recursive = TRUE)
}
