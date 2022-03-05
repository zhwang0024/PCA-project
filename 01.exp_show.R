library("ggplot2")
library("reshape2")
library("maptools")

outdir = "exp_show";
if(!dir.exists(outdir)){
  dir.create(outdir)
}

mycol <-c("#CD0000","#3A89CC","#769C30","#D99536","#7B0078","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")
mypch <-c(21:25,3,4,7,9,8,10,15:18,0:14)
pch=21
col="#1E90FF"

for(xx in c("CO","LO")){
  GRP <- read.table(file = paste0("data/",xx,".grp.xls"), header = T, sep = "\t")
  row.names(GRP) <- GRP$Sample
  
  DAT <- read.table(file = paste0("data/",xx,".gene_fpkm.xls"), header = T, sep = "\t", row.names = 1, comment.char = "")
  EXP <- t(DAT)
  EXP <- log2(EXP+0.01)
  #EXP <- EXP[,colMeans(EXP)>=0]
  EEE <- data.frame(Sample=rownames(EXP),EXP)
  MMM <- melt(data = EEE, id.vars = "Sample")
  
  p <- ggplot(MMM,aes(x=value, color=Sample))
  p <- p + geom_density()
  
  pdf(file=paste0(outdir,"/",xx,".","samples.explevel_density.pdf"),width=10,height=6)
  plot(p)
  dev.off()
  
  A <- data.frame(Sample = GRP$Sample, Group = GRP$G1)
  B <- data.frame(Sample = rownames(t(DAT)), t(DAT))
  C <- merge(x = A, y = B, by = "Sample")
  M <- melt(data = C, id.vars = c("Sample","Group"))
  M$value = log2(M$value+0.01)
  MM <- aggregate(value ~ Group + variable, data = M, FUN = "mean")
  
  p <- ggplot(MM,aes(x=value, color=Group))
  p <- p + geom_density()
  
  pdf(file=paste0(outdir,"/",xx,".","groups.explevel_density.pdf"),width=7,height=6)
  plot(p)
  dev.off()
}
  



