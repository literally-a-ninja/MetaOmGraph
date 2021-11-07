library(edgeR)


## READING THE COUNT MATRIX FROM THE FILE

counts <- read.table("/Users/harshavk/Documents/Harsha/Fall2021/RFiles/limma_counts.csv", header = TRUE, sep = ",")
counts2 <- counts[,-1]
rownames(counts2) <- counts[,1]
counts <- data.matrix(counts2)
head(counts)


## CONVERTING THE COUNT MATRIX TO A DGELIST OBJECT AND CALCULATING THE NORMALIZATION FACTORS

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
d0
dim(d0)


## FILTERING LOW-EXPRESSED GENES BASED ON A CUTOFF

cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
d



## READING GROUPS AND PLOTTING MDS CHART

snames <- colnames(counts) # Sample names
groups <- read.table("/Users/harshavk/Documents/Harsha/Fall2021/RFiles/limma_groups.csv", header = TRUE, sep = ",")
group <- interaction(groups['Groups'])
group

png("/Users/harshavk/Documents/Harsha/Fall2021/RFiles/mds-plot.png")
plotMDS(d, col = as.numeric(group))
dev.off()



## PERFORMING THE LIMMA VOOM OPERATION

mm <- model.matrix(~0 + group)
png("/Users/harshavk/Documents/Harsha/Fall2021/RFiles/voom-plot.png")
y <- voom(d, mm, plot = T)
dev.off()




## Finding differentially expressed genes

fit <- lmFit(y, mm)
head(coef(fit))
contr <- makeContrasts(groupgroup1 - groupgroup2, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

length(which(top.table$adj.P.Val < 0.05))

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "/Users/harshavk/Documents/Harsha/Fall2021/RFiles/differentialexpression.tsv", row.names = F, sep = "\t", quote = F)

