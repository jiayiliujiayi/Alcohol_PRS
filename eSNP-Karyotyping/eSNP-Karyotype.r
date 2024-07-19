# library

library(zoo)
library(gplots)
library(stringr)
library(data.table)
library(optparse)

# Parse options
## Set up command line argument parsing
option_list = list(
  make_option(c("--bamdir"), type="character", default=NULL, help="Directory containing BAM files", metavar="directory")
)

## Parse the command line arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Check if bamdir argument was provided
if (is.null(opt$bamdir)) {
  stop("Error: Please specify the BAM directory using --bamdir", call. = FALSE)
} else {
    bamdir = opt$bamdir
}

## Print the provided directory
print(paste("BAM directory:", opt$bamdir))

## Here you could add code to process BAM files found in the directory
# For instance, listing BAM files:
bam_files = list.files(path = opt$bamdir, pattern = "\\.bam$", full.names = TRUE)
print("BAM files found:")
print(bam_files)

# function sources

function_list_to_source <-
  list.files(path = "./source",
             pattern = "*.R",
             full.names = T)
for (i in 1:length(function_list_to_source)) {
  source(function_list_to_source[i])
}

# PATHS

bin = "/home/jl2791/.conda/envs/alcohol/bin" ## set your bin
CMD_samtools = sprintf("%s/samtools", bin)
CMD_java = sprintf("%s/java", bin)
Picard_Path <- "/home/jl2791/.conda/envs/alcohol/share/picard-2.27.5-0"
CMD_gatk <- sprintf("%s/gatk", bin)

Genome_Fa <- "/projects/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
DIR_Genome_DICT = "/projects/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.dict"
Organism <- "Human"

# DIRs

temp_dir <- "./temp_dir"
# test if temp_dir exists
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
}

plot_output_directory <- "./plot_output"
if (!dir.exists(plot_output_directory)) {
  dir.create(plot_output_directory)
}


# make a eSNPKaryotyping-compatible dbSNP file (one-off) #####
# ! Use $Databases/Genomes/hg38/common_snp150_by_chr ! #####
#Edit_dbSNP_Files(Directory = "snp",
#                 File_Name = "snp151Common",
#                 Organism = "Human")

bam_file_list <-
  list.files(path = opt$bamdir,
             full.names = F, 
             recursive = F,
             pattern = "\\.bam$")
base_names = gsub("\\.bam", "", bam_file_list)
print(bam_file_list)

# prepare sed for creating new header file

# rehead bam files
## Read the chromosome alias table
aliases <- fread("source/hg38.chromAlias.txt")

# Prepare a column for sed replacements for alias to sequenceName
aliases[, alias_to_seq := paste0("s/SN:", `alias names`, "/SN:", `# sequenceName`, "/g")]

# Prepare a column for sed replacements for sequenceName to V1
aliases[, seq_to_v1 := paste0("s/SN:", `# sequenceName`, "/SN:", V1, "/g")]

# Specific replacements for "SN:X" to "SN:chrX" and "SN:Y" to "SN:chrY"
specific_replacements <- c(
  "s/SN:X/SN:chrX/g",
  "s/SN:Y/SN:chrY/g", 
  "s/SN:KI270752.1/SN:chrUn_KI270752v1/g"
)

# Combine all sed commands into a single character vector
all_sed_commands <- c(aliases$alias_to_seq, aliases$seq_to_v1, specific_replacements)

# Write the commands to a shell script
writeLines(all_sed_commands, "./temp_dir/update_header.sed")

# loop through bam list

bam_file_list <-
  list.files(path = opt$bamdir,
             full.names = F, 
             recursive = F,
             pattern = "\\.bam$")
base_names = gsub("\\.bam", "", bam_file_list)

k <- 1
for (k in 1:length(bam_file_list))
    {
    print(sprintf("Started %s at %s ", base_names[k], print(date())))
    ###########################
    ### REHEAD the bam file ###
    ###########################
    
    ## Extract the current BAM header into a file
    
    system(
        sprintf("%s view -H %s/%s > %s/%s_header.sam", 
                CMD_samtools, bamdir, bam_file_list[k], temp_dir, base_names[k]
               )
    )
    
    ## using sed to edit the header file
    
    system(
        sprintf("%s %s %s/%s_header.sam", 
                "sed -i -f", "./temp_dir/update_header.sed", temp_dir, base_names[k]               
               )
    )
    
    ## rehead the bam file using edited header file
    
    system(
        sprintf("%s reheader %s/%s_header.sam %s/%s > %s/reheaded_%s.bam",
                CMD_samtools, temp_dir, base_names[k], bamdir, bam_file_list[k],temp_dir, base_names[k]
               )
    )
    
    ## index the bam file
    
    system(
        command = 
        sprintf("%s index -@ 20 %s/reheaded_%s.bam -o %s/reheaded_%s.bai", 
                CMD_samtools, temp_dir, base_names[k], temp_dir, base_names[k]))
    
    ###########################
    ###      Create VCF     ###
    ###########################
    
    CreateVCF(input_bam_dir = "./temp_dir",
              Genome_Fa,
              DIR_Genome_DICT, 
              Picard_Path,
              CMD_gatk,
              CMD_samtools, 
              base_name=base_names[k],
              out_dir = "./temp_dir")

    ###########################
    ###      Edit VCF       ###
    ###########################
    
    EditVCF(vcfdir = "./temp_dir", 
            base_name = base_names[k], 
            Organism = "Human", 
            out_dir = "./temp_dir") 
    
    ###########################
    ###  Read VCF table in ###
    ###########################
    
    VCF_table <-
    read.delim(file = sprintf('%s/%s.csv',temp_dir, base_names[k]), 
               header = T, sep = "\t",
               quote = "",
               dec = ".")
    VCF_table$chr <- as.numeric(VCF_table$chr)
    VCF_table <- VCF_table[order(VCF_table$chr,
                                 VCF_table$position), ]
    VCF_table <- VCF_table[VCF_table$chr > 0, ]
    
    #######################################
    ###  Return MajorMinorCalc results  ###
    #######################################
    
    MajorMinorCal_results <-
    MajorMinorCalc(Table = VCF_table,
                   minDP = 20, # ! change back to 20 later !
                   maxDP = 1000,
                   minAF = 0.2)
    
    ###########################
    ###  set the file name  ###
    ### and title for the PDF##
    ###########################
    
    pdf_series_name <- base_names[k]
    
    ###########################################
    #####   Plot Allelic ratio along the  #####
    #####  genome for duplication detection ###
    ###########################################
    ## noxy
    pdf(file = sprintf("%s/noxy_%s_genome.pdf", 
                       plot_output_directory, base_names[k]), 
        width=12, height=5, 
        paper = "USr",
        title = base_names[k])
    try({
        PlotGenome_noxy(orderedTable = MajorMinorCal_results,
                   Window = 151,
                   Ylim = 3,
                   PValue = TRUE,
                   Organism = Organism)
        dev.off()
    })
    ## withxy
    pdf(file = sprintf("%s/%s_genome.pdf", 
                       plot_output_directory, base_names[k]), 
        width=12, height=5, 
        paper = "USr",
        title = base_names[k])
    try({
        PlotGenome(orderedTable = MajorMinorCal_results,
                   Window = 151,
                   Ylim = 3,
                   PValue = TRUE,
                   Organism = Organism)
        dev.off()
    })
    
    # tbl_DeletionTable_output <-
    # DeletionTable(Directory = "./temp_dir",
    #               Table = MajorMinorCal_results,
    #               base_name = base_names[k], 
    #               dbSNP_Data_Directory = "./snp",
    #               dbSNP_File_Name = "Edited_snp151Common_chr",
    #               Genome_Fa_dict = Genome_Fa,
    #               Organism = "Human",
    #               source_bam = sprintf("./temp_dir/%s.bam", base_names[k]),
    #               temp_dir = temp_dir)
    # PATH_tbl_deletion = sprintf("%s/%s.txt", "./temp_dir", 
    #                     base_names[k])
    # fwrite(tbl_DeletionTable_output, 
    #        PATH_tbl_deletion, 
    #        quote = F, col.names = T, row.names = F, sep = "\t")
    
    # #######################
    # #####   Plot LOH  #####
    # #######################
    # ## noxy
    # pdf(file = paste(plot_output_directory,
    #                  '/noxy_',
    #                base_names[k],
    #                  "_zygosity_blocks.pdf", sep = ""),
    #     title = base_names[k])
    # Plot_Zygosity_Blocks_noxy(tbl_DeletionTable_output, 
    #                      Window = 1500000, Max = 6, Max2 = 60, Organism = "Human")
    # dev.off()
    
    # ## withxy
    # pdf(file = paste(plot_output_directory,
    #                  '/',
    #                base_names[k],
    #                  "_zygosity_blocks.pdf", sep = ""),
    #     title = base_names[k])
    # Plot_Zygosity_Blocks(tbl_DeletionTable_output, 
    #                      Window = 1500000, Max = 6, Max2 = 60, Organism = "Human")
    # dev.off()
    
    # tbl_DeletionTable_output$snp = 
    # write.table(
    #     sprintf("%s", ), 
    #     quote = F, 
    #     col.names = T, row.names = F, 
    #     sep = "\t"
    # )
    print(sprintf("Finished %s at %s ", base_names[k], print(date())))
    
    }
