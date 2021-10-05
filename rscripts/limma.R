library(edgeR)

counts <- read.delim(mog_data, row.names = 1)

snames <- colnames(counts) # Sample names

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d)

cultivar <- substr(snames, 1, nchar(snames) - 2) 
time <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
group <- interaction(cultivar, time)
group <- interaction(cultivar, time)

#Need to send data to MOG to plot (MDS)

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)

#User input for x, y
contr <- makeContrasts(groupI5.x - groupI5.y, levels = colnames(coef(fit)))

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

#Batch analysis

batch <- factor(rep(rep(1:2, each = 2), 6))
mm <- model.matrix(~0 + group + batch)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
