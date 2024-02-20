setwd("XXX")
set.seed(42)
library(limma)

# import raw
counts_mat = read.delim("PATH_TO_COUNTS")
meta = read.delim("PATH_TO_METADATA") # should contain: PRS, tx, line

## sanity check
# all.equal(rownames(meta), colnames(counts_mat))

# init DGE list
x <- DGEList(counts_mat)

## update sample info
all.equal(rownames(x$samples), meta$id)
x$samples$group = as.factor(paste0(meta$PRS, '.', meta$tx))
x$samples$PRS = as.factor(meta$PRS)
x$samples$tx = as.factor(meta$tx)
x$samples$line = as.factor(meta$line)
x$samples

## update gene info 
geneid <- rownames(x)
genes = data.frame(SYMBOL = geneid)
x$genes <- genes
x

## remove 0 genes
group = x$samples$group
keep.exprs <- filterByExpr(x, group=group) ### remove low expressed genes by group aka PRS.tx
x <- x[keep.exprs,, keep.lib.sizes=FALSE]

## norm counts
x <- calcNormFactors(x, method = "TMM")

# DE
line = x$samples$line
NESTEDLINE = meta %>% select(line, PRS) %>% distinct()
NESTEDLINE_split = split(NESTEDLINE, f = NESTEDLINE$PRS)
NESTEDLINE_split = lapply(NESTEDLINE_split, function(X) mutate(X, nested_line = 1:nrow(X)))
NESTEDLINE = Reduce(rbind, NESTEDLINE_split)
NESTEDLINE$PRS.nested_line = factor(paste0(NESTEDLINE$PRS, ".", NESTEDLINE$nested_line))
NESTEDLINE = merge(meta %>% select(PRS, line), NESTEDLINE, by = c("PRS", "line"), all.x = T)

PRS.nested_line = factor(NESTEDLINE$PRS.nested_line)

design <- model.matrix(~0 + group+PRS.nested_line)
colnames(design) <- gsub("group", "", colnames(design))

## make contrast
contr.matrix <- makeContrasts(HvsL_0 =  High.0 - Low.0, 
                              HvsL_75 = High.75 - Low.75,
                              tx75vs0_L = Low.75 - Low.0, 
                              tx75vs0_H = High.75 - High.0,  
                              levels = colnames(design))
                              
## remove homoscedasticity
v <- voom(x, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)

## limma pipeline for DE
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit)) # summary

## weight to fold-changes in the gene ranking
tfit <- treat(vfit, lfc=0) # treat tests whether log-fold-changes are greater than athreshold rather than merely different to zero
dt <- decideTests(tfit)
summary(dt)

## extract DE genes
HvsL_0 = topTreat(tfit, coef=1, n=Inf) %>% filter(adj.P.Val < 0.05)
HvsL_75 = topTreat(tfit, coef=2, n=Inf) %>% filter(adj.P.Val < 0.05)
tx75vs0_L = topTreat(tfit, coef=3, n=Inf) %>% filter(adj.P.Val < 0.05)
tx75vs0_H = topTreat(tfit, coef=4, n=Inf) %>% filter(adj.P.Val < 0.05)

HPRS_UP_0 = rownames(HvsL_0)[HvsL_0$logFC > 0]
HPRS_UP_75 = rownames(HvsL_75)[HvsL_75$logFC > 0]

HPRS_DOWN_0 = rownames(HvsL_0)[HvsL_0$logFC < 0]
HPRS_DOWN_75 = rownames(HvsL_75)[HvsL_75$logFC < 0]

ETOH75_UP_L = rownames(tx75vs0_L)[tx75vs0_L$logFC > 0]
ETOH75_UP_H = rownames(tx75vs0_H)[tx75vs0_H$logFC > 0]

ETOH75_DOWN_L = rownames(tx75vs0_L)[tx75vs0_L$logFC < 0]
ETOH75_DOWN_H = rownames(tx75vs0_H)[tx75vs0_H$logFC < 0]

HPRS_UP = intersect(HPRS_UP_0, HPRS_UP_75)
HPRS_DOWN = intersect(HPRS_DOWN_0, HPRS_DOWN_75)

ETOH75_UP = intersect(ETOH75_UP_L, ETOH75_UP_H)
ETOH75_DOWN = intersect(ETOH75_DOWN_L, ETOH75_DOWN_H)

# write out
save.image("./99.processed/limmaOut.RData")