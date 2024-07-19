#' EditVCF
#'
#' Edit the VCF (Variant Call Format) file, Creates file with SNPs data at the BAM directory called variantTable.csv
#' @param vcfdir
#' @param base_name The Path of to the BAM Directory, also the VCF file
#' @param Organism "Human" or "Mouse"
#' @param out_dir "temp directory
#' @export
#' @return None

EditVCF <-
function (vcfdir, base_name, Organism, out_dir = "./temp_dir") 
{
    path = sprintf("%s/%s.vcf", vcfdir, base_name)
    read_vcf <- function(file_path) {
        num_lines_to_skip <- sum(readLines(file_path, warn = FALSE) %>% 
            startsWith("##"))
        vcf_data <- fread(file_path, skip = num_lines_to_skip, 
            header = TRUE)
        setnames(vcf_data, names(vcf_data), sub("^#", "", names(vcf_data)))
        return(vcf_data)
    }
    readData = read_vcf(path)
    if (Organism == "Human") {
        readData <- readData[CHROM %in% paste0("chr", c(1:22, 
            "X", "Y")), ]
    }
    else if (Organism == "Mouse") {
        readData <- readData[CHROM %in% paste0("chr", c(1:19, 
            "X", "Y")), ]
    }
    chrRegex = "^chr(\\w+)$"
    infoRegex = "^([01])\\/([01]):(\\d+)\\,(\\d+):(\\d+):\\d+:\\d+\\,\\d+\\,\\d+$"
    chrVector = readData[["CHROM"]]
    posVector = readData[["POS"]]
    infoVector = readData[[base_name]]
    chrNum <- gsub(chrRegex, "\\1", chrVector)
    infoRegex = "^([01])\\/([01]):(\\d+)\\,(\\d+):(\\d+):\\d+:\\d+\\,\\d+\\,\\d+$"
    chrNum <- gsub(chrRegex, "\\1", chrVector)
    if (Organism == "Human") {
        chrNum[chrNum == "X"] = "23"
        chrNum[chrNum == "Y"] = "24"
    }
    if (Organism == "Mouse") {
        chrNum[chrNum == "X"] = "20"
        chrNum[chrNum == "Y"] = "21"
    }
    chrNum <- as.numeric(chrNum)
    Karyotape <- 10 * abs(as.numeric(gsub(infoRegex, "\\1", infoVector)) - 
        as.numeric(gsub(infoRegex, "\\2", infoVector)))
    AD1 <- as.numeric(gsub(infoRegex, "\\3", infoVector))
    AD2 <- as.numeric(gsub(infoRegex, "\\4", infoVector))
    DP <- as.numeric(gsub(infoRegex, "\\5", infoVector))
    posVector <- as.numeric(posVector)
    table <- data.frame(chr = chrNum, position = posVector, AD1 = AD1, 
        AD2 = AD2, DP = DP, Karyotape = Karyotape)
    table <- table[complete.cases(table[, 1:6]), ]
    fileName <- paste0(out_dir, "/", base_name, ".csv")
    write.table(table, file = fileName, sep = "\t", row.names = F, 
        col.names = T, quote = F)
}