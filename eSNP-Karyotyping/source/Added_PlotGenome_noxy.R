# PlotGenome_noxy
#
# Plot Allelic ratio along the genome for duplication detection.
# @param orderedTable The variable containing the output of the MajorMinorCalc function
# @param Window An odd number determining the window of the moving average plot, usually 151
# @param Ylim The plot y axes maximal limit, usually 3
# @param PValue Determine if to add the P-Value bar to the graph, accepts "TRUE" or "FALSE" values
# @param Organism "Human" or "Mouse"
# @export
# @return None

PlotGenome_noxy <-
    function (orderedTable, Window, Ylim, PValue, Organism) 
{
    par(mar = c(5, 4, 4, 2))
    window = Window
    if (Organism == "Human") {
        centromere_pos = c(123400000, 93900000, 90900000, 5e+07, 
            48800000, 59800000, 60100000, 45200000, 4.3e+07, 
            6.1e+07, 10400000, 39800000, 53400000, 35500000, 
            17700000, 17200000, 1.9e+07, 36800000, 25100000, 
            18500000, 26200000, 28100000)
        chr_size = c(248956422, 242193529, 198295559, 190214555, 
            181538259, 170805979, 159345973, 145138636, 138394717, 
            133797422, 135086622, 133275309, 114364328, 107043718, 
            101991189, 90338345, 83257441, 80373285, 58617616, 
            64444167, 46709983, 50818468)
    }
    if (Organism == "Mouse") {
        chr_size = c(195471971, 182113224, 160039680, 156508116, 
            151834684, 149736546, 145441459, 129401213, 124595110, 
            130694993, 122082543, 120129022, 120421639, 124902244, 
            104043685, 98207768, 94987271, 90702639, 61431566, 
            171031299, 91744698)
    }
    chr_total = 0
    for (i in 1:(length(chr_size) - 1)) {
        chr_total <- c(chr_total, sum(chr_size[1:i]))
    }
    genome_size = sum(chr_size)
    orderedTable$position <- orderedTable$position + chr_total[orderedTable$chr]
    col <- rollmedian(orderedTable$chr, window)%%2 + 1
    col[col == 1] <- "dodgerblue4"
    col[col == 2] <- "dodgerblue1"
    plot(rollmedian(orderedTable$position, window), rollmedian(orderedTable$MajorMinor, 
        window), col = "dimgrey", pch = 15, cex = 0.4, ylim = c(1, 
        Ylim), ylab = "Allelic Ratio", typ = "l", xlab = "", 
        xaxt = "n", xlim = c(1, genome_size))
    points(rollmedian(orderedTable$position, window), rollmedian(orderedTable$MajorMinor, 
        window), col = col, pch = 15, cex = 0.4, ylim = c(1, 
        Ylim), xlab = "Chromosomal Position")
    if (isTRUE(PValue)) {
        g = orderedTable$MajorMinor
        ttest = function(x) {
            ttt = t.test(x, g, alternative = "greater")$p.value
            return(ttt)
        }
        tt = rollapply(g, width = window, FUN = ttest)
        tt <- p.adjust(tt, "fdr")
        tt <- 0 - 1 * log10(tt)
        cl <- colorpanel(50000, low = "white", mid = "grey", 
            high = "red")
        pmax <- ceiling(max(tt))
        if (pmax < 2) {
            pmax = 2
        }
        l <- length(rollmedian(orderedTable$MajorMinor, window))
        xx <- rollmedian(orderedTable$position, window)
        for (i in 1:l) {
            lines(c(xx[i], xx[i]), c(2.5, 2.7), col = cl[round(tt[i] * 
                (50000/pmax))])
            if (tt[i] > 2) {
                lines(c(xx[i], xx[i]), c(2.85, 2.75), col = "gray34")
            }
        }
        left <- 5e+07
        for (i in 1:100) {
            rect(left, 3.1, left + 5e+06, 3.3, xpd = T, col = cl[i * 
                500], border = NA)
            left <- left + 5e+06
        }
        left <- 5e+07
        lines(c(left, 5.5e+08), c(3.1, 3.1), xpd = T)
        lines(c(left, 5.5e+08), c(3.3, 3.3), xpd = T)
        lines(c(left, left), c(3.1, 3.3), xpd = T)
        lines(c(5.5e+08, 5.5e+08), c(3.1, 3.3), xpd = T)
        mtext(side = 3, at = 0, "0", line = 0.8, cex = 0.8)
        mtext(side = 3, at = 6e+08, pmax, line = 0.8, cex = 0.8)
        mtext(side = 3, at = 6e+08/2, "-Log10(P-Value)", line = 2.5, 
            cex = 0.8)
    }
    if (Organism == "Human") {
        for (i in 1:24) {
            if (i > 1) {
                lines(c(chr_total[i], chr_total[i]), c(1, Ylim), 
                  col = "gray48")
            }
            lines(c(chr_total[i] + centromere_pos[i], chr_total[i] + 
                centromere_pos[i]), c(1, Ylim), col = "gray55", 
                lty = 4)
            text(chr_total[i] + chr_size[i]/2, Ylim - 0.1, i, 
                cex = 0.8)
        }
    }
    if (Organism == "Mouse") {
        for (i in 1:21) {
            if (i > 1) {
                lines(c(chr_total[i], chr_total[i]), c(1, Ylim), 
                  col = "gray48")
            }
            text(chr_total[i] + chr_size[i]/2, Ylim - 0.1, i, 
                cex = 0.8)
        }
    }
    sum = 0
    for (i in 1:(length(chr_size) - 1)) {
        lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1, 
            1), lwd = 10, col = "black")
        lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1, 
            1), lwd = 8, col = "gray50")
        lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1, 
            1), lwd = 7, col = "gray53")
        lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1, 
            1), lwd = 6, col = "gray59")
        lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1, 
            1), lwd = 4, col = "gray75")
        lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1, 
            1), lwd = 2, col = "gray85")
        lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1, 
            1), lwd = 1, col = "gray90")
        if (Organism == "Human") {
            lines(c(sum + centromere_pos[i], sum + centromere_pos[i]), 
                c(1.01, 0.99), col = "grey13", lwd = 2)
        }
        sum = sum + chr_size[i]
    }
    lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 10, 
        col = "black")
    lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 8, 
        col = "gray50")
    lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 7, 
        col = "gray53")
    lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 6, 
        col = "gray59")
    lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 4, 
        col = "gray75")
    lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 2, 
        col = "gray85")
    lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 1, 
        col = "gray90")
}