#' DeletionTable
#'
#' Intersect the observed SNPs with the common SNPs table from dbSNPs, Creates file with the LOH data called Deletions.txt
#' @param Directory The Path of to the vcf Directory
#' @param Table The variable containing the output of the MajorMinorCalc function
#' @param dbSNP_Data_Directory The path for the directory where the edited dbSNP file are (created previously by the Edit_dbSNP_Files function)
#' @param dbSNP_File_Name The edited dbSNP file names, without the chromosome number
#' @param Genome_Fa_dict the path for the dictionary file created for the whole genome FASTA file.
#' @param Organism "Human" or "Mouse"
#' @param source_bam original BAM used to generate the VCF file
#' @export
#' @return None

DeletionTable <-
  function(Directory,
           Table,
           base_name, 
           dbSNP_Data_Directory,
           dbSNP_File_Name,
           Genome_Fa_dict,
           Organism,
           source_bam,
           temp_dir) {
    #print("Reading SNPs table")
    i <- 1

    if (Organism == "Human") {mx = 25}
    if (Organism == "Mouse") {mx = 22}

    #print(getwd())

    Genome_Fa_dict <-
      str_replace(string = Genome_Fa_dict,
                  pattern = "\\.fa$",
                  replacement = "\\.dict")
    dict <-
      read.csv(Genome_Fa_dict,
               as.is = T)
    dict_type = grep(pattern = "chr",
                     x = dict[1, 1])
    #print(dict_type)
    if (length(dict_type) == 0) {
      dict_type = 0
    }


    while (i < mx) {
      #print(paste("Chromosome:", i))
      chrTable = read.delim(paste(dbSNP_Data_Directory,
                                  '/',
                                  dbSNP_File_Name,
                                  i,
                                  ".txt", 
                                  sep = ""))
      if (i == 1 ) {
        snpTable = chrTable
      } else if (i > 1) {
        snpTable = rbind(snpTable,
                         chrTable)
      }
      i = i + 1
    }


    table2 <- Table
    colnames(table2)[2] <- "start"
    # print("Merging Tables")
    x <- merge(snpTable,
               table2,
               by = c("chr", "start"),
               all.x = T)
    x <- x[order(x$chr, x$start),]


    # print("Indexing BAM File")
    # system("samtools index accepted_hits.bam")

    i <- 1
    while (i < mx) {
      #print(paste("Chromosome ",
      #            i,
      #            "| ",
      #            Sys.time(),
      #            sep = ""))
      x1 <- x[x$chr == i,]


      if (dict_type == 0) {
        loc <-
          paste(x1$chr[1],
                ":",
                x1$start[1],
                "-",
                x1$start[dim(x1)[1]],
                sep = "")
        if (Organism == "Human") {
          if (i == 23) {
            loc <-
              paste("X:",
                    x1$start[1],
                    "-",
                    x1$start[dim(x1)[1]],
                    sep = "")
          }
          if (i == 24) {
            loc <-
              paste("Y:",
                    x1$start[1],
                    "-",
                    x1$start[dim(x1)[1]],
                    sep = "")
          }}

        if (Organism=="Mouse") {
          if (i == 20) {
            loc <-
              paste("X:",
                    x1$start[1],
                    "-",
                    x1$start[dim(x1)[1]],
                    sep = "")
          }
          if (i == 21) {
            loc <-
              paste("Y:",
                    x1$start[1],
                    "-",
                    x1$start[dim(x1)[1]],
                    sep = "")}}
      }

      if (dict_type == 1) {
        loc <-
          paste("chr",
                x1$chr[1],
                ":",
                x1$start[1],
                "-",
                x1$start[dim(x1)[1]],
                sep = "")
        if (Organism == "Human") {
          if (i == 23) {
            loc <-
              paste("chrX:",
                    x1$start[1],
                    "-",
                    x1$start[dim(x1)[1]],
                    sep = "")}
          if (i == 24) {
            loc <-
              paste("chrY:",
                    x1$start[1],
                    "-",
                    x1$start[dim(x1)[1]],
                    sep = "")
          }}
        if (Organism == "Mouse") {
          if (i == 20) {
            loc <-
              paste("chrX:",
                    x1$start[1],
                    "-",
                    x1$start[dim(x1)[1]],
                    sep = "")}
          if (i == 21) {
            loc <-
              paste("chrY:",
                    x1$start[1],
                    "-",
                    x1$start[dim(x1)[1]],
                    sep = "")
          }}
        #print(loc)
      }



      command <- 
        sprintf("%s depth -r %s %s > %s/%s_reads-per-position.txt",
                CMD_samtools, 
                loc,
                source_bam, 
                temp_dir,
                base_name
               )

      system(command)

cmd_awk = sprintf(
    "awk -F ' ' '$3 > 5 {print $0}' %s/%s_reads-per-position.txt > %s/%s_reads-per-position2.txt", 
    temp_dir,
    base_name, 
    temp_dir,
    base_name
)
system(cmd_awk)
        
cmd_awk2 <- sprintf(
    "awk -F ' ' '$3 > 5 {print $0}' %s/%s_reads-per-position.txt > %s/%s_reads-per-position2.txt", 
    temp_dir,
    base_name, 
    temp_dir,
    base_name
)

system(cmd_awk2)

        
      chr <-
        read.delim(
            sprintf("%s/%s_reads-per-position2.txt", 
                   temp_dir, base_name), 
            header = F)
      colnames(chr) <-
        c("chr", "start", "Depth")
      x2 <-
        merge(x1, chr,
              by = "start")
      x2 <-
        cbind(x2,
              Depth_group = x2$Depth)
      x2$Depth_group[x2$Depth < 50] <- "20-50"
      x2$Depth_group[x2$Depth > 49 & x2$Depth < 100] <- "50-100"
      x2$Depth_group[x2$Depth > 99 & x2$Depth < 200] <- "100-200"
      x2$Depth_group[x2$Depth > 199 & x2$Depth < 500] <- "200-500"
      x2$Depth_group[x2$Depth > 499] <- ">500"
      x2$Depth_group <- as.factor(x2$Depth_group)
      if (i == 1 ) {tbl <- x2}
      if (i > 1 ) {
        tbl <- rbind(tbl,
                    x2)
      }
      i = i + 1
    }


    #print("Writing Table")

    tbl[is.na(tbl)] = 0
    # tbl <-
    #   tbl[complete.cases(tbl[, 1:ncol(tbl)]), ]
    write.table(tbl, 
                sprintf("%s/%s_Deletions.txt", 
                        temp_dir, base_name), 
                sep = "\t",
                row.names = F,
                quote = F)
    return(tbl)
  }
