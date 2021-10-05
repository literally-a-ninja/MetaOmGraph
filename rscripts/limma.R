library(edgeR)

counts <- read.delim("all_counts.xtx", row.names = 1)
d0 <- DGEList(counts)

cutoff <- 1
drop <- which(apply(cpm(d0, 1, max) < cutoff)
d <- d0[-drop,]
dim(d)

snames <- colnames(counts)

cultivar <- substr(snames, 1, nchar(snames) -2)
time <- substr(snames, nchar(snames) -1, nchar(snames) -1)

group <- interaction(cultivar, time)

mm <- model.matrix(~0 + group)

y <- voom(d, mm, plot = T)

tmp <- voom(d0, mm, plot = T)

fit <- lmFit(y, m)

contr <- makeContrasts(groupI5.9 - groupI5.6, levels = colnames(coef(fit)))

tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf)

