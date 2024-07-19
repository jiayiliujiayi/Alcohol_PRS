Calculate_chrom_ratios_and_pvalues <- function(Table, Window, Max, Max2, Organism) {
  
  max = Max
  max2 = Max2
  window = Window
  tbl = Table
  results = data.frame(chr = integer(), part = integer(), rat = numeric(), p_value = numeric())
  
  if (Organism == "Human") {
    centromere_pos <- c(123400000,  93900000,  90900000,  50000000,  48800000,  59800000,  60100000,  45200000,
                        43000000,  61000000,  10400000,  39800000,  53400000,  35500000,  17700000,  17200000,
                        19000000,  36800000,  25100000,  18500000,  26200000,  28100000,  12000000,  15000000)
    chr_size <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                  133797422, 135086622, 133275309, 114364328, 107043718, 101991189,  90338345,  83257441,  80373285,
                  58617616,  64444167,  46709983,  50818468, 156040895,  57227415)
    mx = 24
  }
  
  if (Organism == "Mouse") {
    chr_size = c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110,
                 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639,
                 61431566, 171031299, 91744698)
    centromere_pos = rep(200000000, 21)
    mx = 21
  }
  
  tbl$snp = tbl$AF1 > 0
  total_homo = sum(tbl$snp == FALSE)
  total_hetro = sum(tbl$snp == TRUE)
  rat = NULL

  for (i in 1:mx) {
    if (i < (mx - 1)) {
      chr = i
      tbl2 = tbl[tbl$chr.x == chr, ]
      
      d = tbl2[tbl2$start < centromere_pos[i], ]
      homo = sum(d$snp == FALSE)
      hetro = sum(d$snp == TRUE)
      rat = c(rat, homo / hetro)
      
      d = tbl2[tbl2$start > centromere_pos[i], ]
      homo = sum(d$snp == FALSE)
      hetro = sum(d$snp == TRUE)
      rat = c(rat, homo / hetro)
    }
  }

  rat[is.infinite(rat)] = total_homo / total_hetro * 5.5
  pval = sapply(rat, function(r) {
    if (!is.na(r)) {
      t.test(rat, mu = r)$p.value
    } else {
      NA
    }
  })

  pval = p.adjust(pval, method = "fdr")
  pval = log10(pval) * (-1)
  rat[is.infinite(rat)] = total_homo / total_hetro
  rat[is.na(rat)] = total_homo / total_hetro

  for (i in 1:(mx )) {
    for (j in 1:2) {
      results = rbind(results, data.frame(chr = i, part = j, rat = rat[(i - 1) * 2 + j], p_value = pval[(i - 1) * 2 + j]))
    }
  }

  return(results)
}


calculate_nsnps <- function(Table, Window, Max, Max2, Organism) {
  
  max = Max
  max2 = Max2
  window = Window
  tbl = Table
  results = data.frame(chromosome = integer(), relative_window_start = integer(), relative_window_end = integer(), num_heterozygous_snps = integer(), num_homozygous_snps = integer())

  if (Organism == "Human") {
    centromere_pos <- c(123400000,  93900000,  90900000,  50000000,  48800000,  59800000,  60100000,  45200000,
                        43000000,  61000000,  10400000,  39800000,  53400000,  35500000,  17700000,  17200000,
                        19000000,  36800000,  25100000,  18500000,  26200000,  28100000,  12000000,  15000000)

    chr_size <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                  133797422, 135086622, 133275309, 114364328, 107043718, 101991189,  90338345,  83257441,  80373285,
                  58617616,  64444167,  46709983,  50818468, 156040895,  57227415)
    mx = 24
  }
  
  if (Organism == "Mouse") {
    chr_size = c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110,
                 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639,
                 61431566, 171031299, 91744698)
    centromere_pos = rep(200000000, 21)
    mx = 21
  }
  
  tbl$snp = tbl$AF1 > 0

  for (i in 1:mx) {
    chr = i
    tbl2 = tbl[tbl$chr.x == chr, ]
    tr = tbl2[tbl2$snp == TRUE, ]
    fl = tbl2[tbl2$snp == FALSE, ]
    
    win = ceiling(chr_size[i] / window)
    pos = centromere_pos[i]
    
    for (j in 1:win) {
      top = j * window
      bottom = (j - 1) * window
      num_fl = sum(fl$start >= bottom & fl$start <= top)
      num_tr = sum(tr$start >= bottom & tr$start <= top)
      
      if (num_tr > max) {
        num_tr = max
      }
      if (num_fl > max2) {
        num_fl = max2
      }
      
      results = rbind(results, data.frame(chromosome = chr, relative_window_start = bottom, relative_window_end = top, num_heterozygous_snps = num_tr, num_homozygous_snps = num_fl))
      
      pos = pos - window
    }
  }

  return(results)
}

