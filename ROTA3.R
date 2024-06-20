# Rotation 3
# Author: Ahona Roy, National Institute of Biomedical Genomics
# Supervisor: Dr. Saroj K Mohapatra
# September 2020

# Loading the necessary libraries
library(Biobase)
library(GEOquery)
library(genefilter)
library(limma)
library(gplots)
library(org.Hs.eg.db)
library(VennDiagram)
library(RColorBrewer)

# Loading dataset1 : 1ng/kg endotoxin dose
studyid = "GSE36177"
fn = paste0(studyid, ".rda")
if(file.exists(file = fn)) {
  cat("Reading data from file ...")
  load(fn)
  cat(" done!\n")
} else {
  eset1 = getGEO(studyid, GSEMatrix = TRUE, getGPL = FALSE)
  eset1 = eset1[[1]] 
  save(eset1, file = fn)
}

# Setting the group information
eset1$Group = c("Endotoxin","Zcontrol")[as.integer(eset1$`time:ch1`=="0 hours")+1]
rm(fn, studyid)

# Loading dataset2 : 2ng/kg endotoxin dose
studyid = "GSE108685"
fn = paste0(studyid, ".rda")
if(file.exists(file = fn)) {
  cat("Reading data from file ...")
  load(fn)
  cat(" done!\n")
} else {
  eset2 = getGEO(studyid, GSEMatrix = TRUE, getGPL = FALSE)
  eset2 = eset2[[1]] 
  save(eset2, file = fn)
}

# Subsetting the 0 hr and 4 hr timepoints from placebo group in dataset2 
dose = which(eset2$`stem cell dose:ch1` == "placebo")
eset2 = eset2[, dose]

t1 = which(eset2$`timepoint pre(-)/post(+) lps:ch1`=="0 hours")
t2 = which(eset2$`timepoint pre(-)/post(+) lps:ch1`=="+4 hours")

time = c(t1, t2)
eset2 = eset2[, time]

eset2$Group = c("Endotoxin", "Zcontrol")[as.integer(eset2$`timepoint pre(-)/post(+) lps:ch1`=="0 hours")+1]
rm(fn, studyid, dose, t1, t2, time)

# Performing a t-test for each gene between endotoxin and 
# control for dataset1
rtt1 = rowttests(eset1, fac = "Group")

# Selecting genes with p-value less than 0.05 and log-fold 
# change of 1
lfc = rtt1[, "dm"] 
pval = rtt1[, "p.value"]
stat = rtt1[, "statistic"]

sel = pval < 0.05 & abs(lfc) >= 1

# Number of significantly altered genes
sum(sel)

# 1052

# Correcting the p-values for multiple testing with fdr method 
# Selecting genes with adjusted p-value less than 0.05 
# and log-fold change of 1
pfdr = p.adjust(pval, method = "fdr")
selfdr = pfdr < 0.05 & abs(lfc) >= 1

# Number of significantly altered genes with adjusted p value
sum(selfdr)

# 1051

rm(lfc, pval, stat, sel, pfdr, selfdr)

# Performing a t-test for each gene between endotoxin and 
# control for dataset2
rtt2 = rowttests(eset2, fac = "Group")

# Selecting genes with p-value less than 0.05 and log-fold 
# change of 1 
lfc = rtt2[, "dm"] 
pval = rtt2[, "p.value"]
stat = rtt2[, "statistic"]

sel = pval < 0.05 & abs(lfc) >= 1

# Number of significantly altered genes
sum(sel)

# 4786

# Correcting the p-values for multiple testing with fdr method 
# Selecting genes with adjusted p-value less than 0.05 
# and log-fold change of 1
pfdr = p.adjust(pval, method = "fdr")
selfdr = pfdr < 0.05 & abs(lfc) >= 1

# Number of significantly altered genes with adjusted p value
sum(selfdr)

# 4768

rm(rtt1, rtt2, lfc, pval, stat, sel, pfdr, selfdr)

# Reading the feature data for both the datasets
fdata1 = read.table(file = "GSE36177_ID.txt", header=TRUE, stringsAsFactors = FALSE, 
                    sep = "\t", fill = TRUE)

fdata2 = read.table(file = "GSE108685_ID.txt", header=TRUE, stringsAsFactors = FALSE, 
                    sep = "\t", fill = TRUE)

# Finding the common probe IDs between the assay data 
# and the feature data for both the datasets
pid1 = intersect(rownames(eset1), fdata1$ID)
pid2 = intersect(rownames(eset2), fdata2$ID)

# Subsetting the feature data for both the datasets 
# based on the common probe IDs having an Entrez ID
data1 = subset(fdata1, ID %in% pid1)
egid1 = data1$Entrez_Gene_ID
egid1 = egid1[!is.na(egid1)]
data1 = subset(data1, Entrez_Gene_ID %in% egid1)
rownames(data1) = data1$ID

data2 = subset(fdata2, ID %in% pid2)
egid2 = as.integer(data2$Entrez.Gene)
egid2 = egid2[!is.na(egid2)]
data2 = subset(data2, Entrez.Gene %in% egid2)
rownames(data2) = data2$ID

# Subsetting the assay data for both the datasets based on the 
# common probe IDs
exp1 = exprs(eset1)
exp1 = subset(exp1, rownames(exp1) %in% pid1)
e1 = rownames(exp1)

exp2 = exprs(eset2)
exp2 = subset(exp2, rownames(exp2) %in% pid2)
e2 = rownames(exp2)

# Mapping the probe ID to the Entrez ID for both the datasets
for(i in 1 : length(egid1)) {
  m = match(rownames(data1)[i], e1, nomatch = NA)
  e1 = replace(e1, m, egid1[i])
} 
e1 = as.integer(e1)

for(i in 1 : length(egid2)) {
  m = match(rownames(data2)[i], e2, nomatch = NA)
  e2 = replace(e2, m, egid2[i])
}
e2 = as.integer(e2)

rm(pid1, pid2, egid1, egid2)

# Adding the Entrez ID dataframe to the matrix and removing the probes 
# with no Entrez ID
exp1 = cbind(exp1, e1)
exp2 = cbind(exp2, e2) 

exp1 = na.omit(exp1)
e1 = e1[!is.na(e1)]
exp2 = na.omit(exp2)
e2 = e2[!is.na(e2)]

# Removing the duplicate probe IDs having the same Entrez ID and keeping the 
# one with the highest variance between the samples
for(i in 1 : length(e1)) {
  pid = which((exp1)[,"e1"] == e1[i])
  if (length(pid) == 1) {
    exp1 = exp1
  } else {
    x = exp1[pid, 1:(ncol(exp1) - 1)]
    rv = rowVars(x)
    sel = names(rv)[-which.max(rv)]
    exp1 = exp1[-which(rownames(exp1) %in% sel),]
  }
}

for(i in 1 : length(e2)) {
  pid = which((exp2)[,"e2"] == e2[i])
  if (length(pid) == 1) {
    exp2 = exp2
  } else {
    x = exp2[pid, 1:(ncol(exp2) - 1)]
    rv = rowVars(x)
    sel = names(rv)[-which.max(rv)]
    exp2 = exp2[-which(rownames(exp2) %in% sel),]
  }
}

rm(x, m, i, rv, sel, pid)

# Replacing the probe IDs with Entrez ID
e1 = exp1[,"e1"]
rownames(exp1) = e1
exp1 = exp1[, -31]

e2 = exp2[,"e2"]
rownames(exp2) = e2
exp2 = exp2[, -18]

rm(e1, e2)

# Venn diagram showing the genes in the two datasets
common=intersect(rownames(exp1), rownames(exp2))

grid.newpage() 
draw.pairwise.venn(area1 = nrow(exp1), area2 = nrow(exp2), 
                   cross.area = length(common), category = c("1 ng/kg", "2 ng/kg"), 
                   cat.pos = c(335, 25), cat.dist = rep(0.025, 2), cat.cex = 2.5, 
                   scaled = FALSE, cex = 4)
grid.text("Number of genes in the 2 datasets of different endotoxin dosage", 
          y = 0.9, gp = gpar(col = "black", fontsize = 35))

# Subsetting the two datasets based on the common genes between them 
exp1.common = exp1[common,]
rownames(exp1.common) = common 

n = seq(from = 1, to = ncol(exp1.common), by = 2)
for(i in n) {
  colnames(exp1.common)[i] = "Zcontrol"
  colnames(exp1.common)[i+1] = "Endotoxin"
}
rm(i, n)

exp2.common = exp2[common,]
rownames(exp2.common) = common
colnames(exp2.common)[1:10] = "Zcontrol"
colnames(exp2.common)[11:17] = "Endotoxin"

# Performing a t-test for each of the gene between endotoxin and 
# control for the subsetted dataset1 
rtt1 = rowttests(exp1.common, factor(colnames(exp1.common)))

# Selecting genes with p-value less than 0.05 and log-fold 
# change of 1
lfc1 = rtt1[, "dm"] 
pval = rtt1[, "p.value"]
stat = rtt1[, "statistic"]

sel = pval < 0.05 & abs(lfc1) >= 1

# Number of significantly altered genes
sum(sel)

# 308

# Correcting the p-values for multiple testing with fdr method 
# Selecting genes with adjusted p-value less than 0.05 
# and log-fold change of 1
pfdr = p.adjust(pval, method = "fdr")
selfdr = pfdr < 0.05 & abs(lfc1) >= 1

# Number of significantly altered genes with adjusted p value
sum(selfdr)

# 307

# Number of significantly upregulated and downregulated genes with 
# adjusted p value in dataset1
signif = pfdr < 0.05
up = lfc1 >= 1
down = lfc1 <= -1

signif.up1 = signif & up
sum(signif.up1)

# 218

signif.dn1 = signif & down
sum(signif.dn1)

# 89

rm(pval, stat, selfdr, pfdr, sel, signif, up, down)

# Performing a t-test for each of the gene between endotoxin and 
# control for the subsetted dataset2 
rtt2 = rowttests(exp2.common, factor(colnames(exp2.common)))

# Selecting genes with p-value less than 0.05 and log-fold 
# change of 1
lfc2 = rtt2[, "dm"] 
pval = rtt2[, "p.value"]
stat = rtt2[, "statistic"]

sel = pval < 0.05 & abs(lfc2) >= 1

# Number of significantly altered genes
sum(sel)

# 957

# Correcting the p-values for multiple testing with fdr method 
# Selecting genes with adjusted p-value less than 0.05 
# and log-fold change of 1
pfdr = p.adjust(pval, method = "fdr")
selfdr = pfdr < 0.05 & abs(lfc2) >= 1

# Number of significantly altered genes with adjusted p value
sum(selfdr)

# 952

# Number of significantly upregulated and downregulated genes with 
# adjusted p value in dataset2
signif = pfdr < 0.05
up = lfc2 >= 1
down = lfc2 <= -1

signif.up2 = signif & up
sum(signif.up2)

# 343

signif.dn2 = signif & down
sum(signif.dn2)

# 609

rm(pval, stat, selfdr, pfdr, sel, signif, up, down)

# Subsetting the significantly altered genes in dataset1 
exp1.up = which(signif.up1)
exp1.dn = which(signif.dn1)

# Subsetting the significantly altered genes in dataset2
exp2.up = which(signif.up2)
exp2.dn = which(signif.dn2)

# Creating a matrix for the two datasets showing change in expression
result = matrix(data = 0, nrow = nrow(exp1.common), ncol = 2)
rownames(result) = rownames(exp1.common)
colnames(result) = c("dataset1", "dataset2")

result[exp1.up, 1] = 1
result[exp1.dn, 1] = -1

result[exp2.up, 2] = 1
result[exp2.dn, 2] = -1

# Venn diagram showing differentially expressed genes in both datasets
vennDiagram(result, include = c("up","down"), counts.col = c("red", "blue"),
            names = c("1 ng/kg", "2 ng/kg"), cex = 4)
title(main = "Differentially expressed genes in the two datasets", cex.main = 3)

rm(result, exp1.up, exp1.dn, exp2.up, exp2.dn)

# Genes showing monotonic increase in expression
deg.up = which(signif.up1 & signif.up2)
m.up = which((lfc2[deg.up] - lfc1[deg.up]) > 0)
m.up = deg.up[m.up]
length(m.up)

# 127

# Genes showing monotonic decrease in expression
deg.down = which(signif.dn1 & signif.dn2)
m.down = which((lfc2[deg.down] - lfc1[deg.down]) < 0)
m.down = deg.up[m.down]
length(m.down)

# 69

rm(common, signif.up1, signif.up2, signif.dn1, signif.dn2, deg.up, deg.down)


# Genes showing a linear monotonic increase in expression
x = 2*lfc1[m.up]
y = lfc2[m.up]

l.up = which(y >= (x-0.2*x))
l.up = m.up[l.up]
length(l.up)

# 36 

# Genes showing a linear monotonic decrease in expression
x = 2*lfc2[m.down]
y = lfc1[m.down]

l.down = which(y >= (x-0.2*x))
l.down = m.down[l.down]
length(l.down)

# 1

rm(x, y)

# Entrez IDs of genes showing linear increase in expression
rownames(exp1.common)[l.up]

# [1] "118788" "115123" "23253"  "27071"  "51363"  "3301"   "57162"  "3717"   "116071" "306"   
# [11] "51738"  "2040"   "23001"  "2633"   "91694"  "83463"  "80315"  "838"    "83666"  "7357"  
# [21] "1432"   "3433"   "9111"   "1378"   "3422"   "55509"  "9895"   "64135"  "23586"  "9943"  
# [31] "3240"   "25932"  "6648"   "4683"   "2181"   "2907" 

# Entrez IDs of genes showing linear decrease in expression
rownames(exp1.common)[l.down]

# [1] "8807"

# Subsetting the log fold change (lfc) values of the genes showing
# linear change in expression from both the datasets
lfc.d1 = c(lfc1[l.up], lfc1[l.down])
lfc.d2 = c(lfc2[l.up], lfc2[l.down])

# Creating a matrix containing the lfc values of those genes from both the datasets
result = matrix(data = 0, nrow = 37, ncol = 2)
rownames(result) = c(rownames(exp1.common)[l.up], rownames(exp1.common)[l.down])
colnames(result) = c("1 ng/kg", "2 ng/kg")

result[, 1] = lfc.d1
result[, 2] = lfc.d2

# Getting the gene symbol and name for each Entrez ID 
gsym = links(org.Hs.egSYMBOL[rownames(result)])[,2]
gname = links(org.Hs.egGENENAME[rownames(result)])[,2]
rname = paste(gsym, " [", gname, "]", " | ", rownames(result), sep = "")

# Ordering the gene names 
dx = apply(X = result, MARGIN = 1, FUN = diff)
o = order(result[,1] + sign(dx))
plotdata = result[o,]
rname = rname[o]

# Heatmap showing linear change in expression of the genes
mycol = colorRampPalette(c("white","darkblue"))(24)

heatmap.2(plotdata, col = mycol,
          cellnote = formatC(plotdata, digits = 2), notecol = "black",
          Rowv = FALSE, dendrogram = "none",labRow = rname,
          margins = c(4, 35), cexRow = 1, cexCol = 1, srtCol = 0, adjCol = c(NA, 1),
          key = TRUE, keysize = 1, key.title = "Colour Key", key.ylab =  NA,
          key.xlab = "lfc value", density.info = "none", trace = "none")

rm(mycol, plotdata, dx, o, rname, gsym, gname, l.up, l.down, 
   m.up, m.down, lfc.d1, lfc.d2, lfc1, lfc2)


