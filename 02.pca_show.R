library("ggplot2")
library("reshape2")
library("maptools")

outdir = "pca_show";
if(!dir.exists(outdir)){
  dir.create(outdir)
}

mycol <-c("#CD0000","#3A89CC","#769C30","#D99536","#7B0078","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")
mypch <-c(21:25,3,4,7,9,8,10,15:18,0:14)
pch=21
col="#1E90FF"

for(zz in c(0.5,1,2,3,4)){
  if(!dir.exists(paste0(outdir,"/",zz))){
    dir.create(paste0(outdir,"/",zz))
  }
  if(!dir.exists(paste0(outdir,"/",zz,"/pca_time"))){
    dir.create(paste0(outdir,"/",zz,"/pca_time"))
  }
  if(!dir.exists(paste0(outdir,"/",zz,"/pca_2d"))){
    dir.create(paste0(outdir,"/",zz,"/pca_2d"))
  }
  if(!dir.exists(paste0(outdir,"/",zz,"/pca_exp"))){
    dir.create(paste0(outdir,"/",zz,"/pca_exp"))
  }
  if(!dir.exists(paste0(outdir,"/",zz,"/pca_contribution"))){
    dir.create(paste0(outdir,"/",zz,"/pca_contribution"))
  }
  for(xx in c("CO","LO")){
    DAT <- read.table(file = paste0("data/",xx,".gene_fpkm.xls"), header = T, sep = "\t", row.names = 1, comment.char = "")
    EXP <- DAT[apply(DAT,1,mean)>zz,]
    
    EXP <- t(EXP)
    EXP <- log2(EXP+0.01)
    GRP <- read.table(file = paste0("data/",xx,".grp.xls"), header = T, sep = "\t")
    row.names(GRP) <- GRP$Sample
    
    
    pca   <- prcomp(x = EXP, scale = T, center = T)
    percent_variance <- summary(pca)$importance[2,]*100
    
    d <- data.frame(Sample = row.names(pca$x),pca$x)
    d <- merge(x = GRP, y = d)
    
    m <- melt(data = d, id.vars = c("Sample","Group"))
    m$Group <- factor(m$Group,levels = c("E10","E14","E18","1D","1W","3W","5W"))
    m$Time  <- 1
    m$Time[m$Group == "E14"] = 2
    m$Time[m$Group == "E18"] = 3
    m$Time[m$Group == "1D"]  = 4
    m$Time[m$Group == "1W"]  = 5
    m$Time[m$Group == "3W"]  = 6
    m$Time[m$Group == "5W"]  = 7
    
    for(a in 1:length(percent_variance)){
      pdf(file=paste0(outdir,"/",zz,"/pca_time/",xx,".pca_",a,".pdf"),width=6,height=4.5)
      p <- ggplot(m[m$variable==paste0("PC",a),], aes(x = Time, y = value, group = "variable"))
      p <- p + geom_smooth(col="grey90", size = 2, se = F)
      p <- p + geom_point(shape = 15, color = "#ff33cc", size = 3)
      p <- p + scale_x_continuous(breaks=c(1,2,3,4,5,6,7), labels=c("E10","E14","E18","1D","1W","3W","5W"))
      p <- p + ylab(paste0("Projection on PC",a," axis"))
      p <- p + theme_bw()
      p <- p + theme(
        panel.grid = element_blank(),
        axis.text  = element_text(size = 12),
        axis.title = element_text(size = 15)
        )
      print(p)
      dev.off()
    }
    pdf(file=paste0(outdir,"/",zz,"/pca_time/",xx,".pca_all.pdf"),width=16,height=13)
    p <- ggplot(m, aes(x = Time, y = value, group = variable))
    p <- p + geom_smooth(col="grey90", size = 2, se = F)
    p <- p + geom_point(shape = 15, color = "#ff33cc", size = 3)
    p <- p + scale_x_continuous(breaks=c(1,2,3,4,5,6,7), labels=c("E10","E14","E18","1D","1W","3W","5W"))
    p <- p + ylab(paste0("Projection on PC axis"))
    p <- p + facet_wrap(vars(variable), scales = "free_y")
    p <- p + theme_bw()
    p <- p + theme(
      strip.text = element_text(size = 15),
      panel.grid = element_blank(),
      axis.text  = element_text(size = 12),
      axis.title = element_text(size = 15)
    )
    print(p)
    dev.off()
    
    # for(w in 1:4){
    #   plot(x = colMeans(EXP), y = pca$rotation[,w])
    #   cat(cor(colMeans(EXP),pca$rotation[,w]),"\n")
    # }
    for(ll in 1:5){
      
      zt <- data.frame(gene = row.names(pca$rotation),contribution = pca$rotation[,ll])
      nn <- ceiling(nrow(zt) * 0.05)
      gene_up   <- zt[order(zt$contribution,decreasing = T),][1:nn,]
      gene_down <- zt[order(zt$contribution,decreasing = F),][1:nn,]
      cut_up    <- gene_up$contribution[nn]
      max_up    <- gene_up$contribution[1]
      cut_down  <- gene_down$contribution[nn]
      min_down  <- gene_down$contribution[1]
      
      cat(paste(collapse = "\n",gene_up$gene),file=paste0(outdir,"/",zz,"/pca_contribution/",xx,".pca_",ll,".up.list"))
      cat(paste(collapse = "\n",gene_down$gene),file=paste0(outdir,"/",zz,"/pca_contribution/",xx,".pca_",ll,".down.list"))
      
      pdf(file=paste0(outdir,"/",zz,"/pca_contribution/",xx,".pca_",ll,".pdf"),width=6,height=4.5)
      p <- ggplot(zt, aes(x = contribution))
      #p <- p + geom_histogram(binwidth = 0.0001)
      p <- p + geom_density(alpha = 0.5, color = 'orange',size=1)
      p <- p + annotate("rect", xmin=cut_up, xmax=max_up, ymin = 0, ymax = Inf, alpha = .2)
      p <- p + annotate("rect", xmin=min_down, xmax=cut_down, ymin = 0, ymax = Inf, alpha = .2)
      p <- p + geom_vline(aes(xintercept = 0), linetype="dashed",colour="blue",size=1)
      #p <- p + geom_rect(aes(xmin=cut_up, xmax=max_up, ymin=0, ymax=Inf), fill='grey',alpha=0.7)
      #p <- p + geom_rect(aes(xmin=min_down, xmax=cut_down, ymin=0, ymax=Inf),fill='grey',alpha=0.7)
      #p <- p + geom_vline(aes(xintercept = cut_up), linetype="dashed",colour="blue",size=0.7)
      #p <- p + geom_vline(aes(xintercept = gene_down))
      p <- p + theme_bw()
      p <- p + theme(
        panel.grid = element_blank(),
        axis.text  = element_text(size = 12),
        axis.title = element_text(size = 15)
      )
      print(p)
      dev.off()
    }
    
    r1 <- data.frame(gene = row.names(pca$rotation),pca$rotation)
    r2 <- data.frame(gene = row.names(t(EXP)),mean = colMeans(EXP))
    rr <- merge(x = r2, y = r1, by = "gene")
    rr <- melt(data = rr, id.vars = c("gene","mean"))
    
    pdf(file=paste0(outdir,"/",zz,"/pca_exp/",xx,".scatter.pdf"),width=16,height=13)
    p <- ggplot(rr, aes(x = mean, y = value, group = variable))
    p <- p + geom_smooth(col="grey90", size = 2, se = F)
    p <- p + geom_point(size = 0.5)
    p <- p + xlab(paste0("Average gene expression"))
    p <- p + ylab(paste0("Component of the principle axis"))
    p <- p + facet_wrap(vars(variable), scales = "free")
    p <- p + theme_bw()
    p <- p + theme(
      strip.text = element_text(size = 15),
      panel.grid = element_blank(),
      axis.text  = element_text(),
      axis.title = element_text(size = 15)
    )
    print(p)
    dev.off()
    
    nn <- ceiling(ncol(EXP) * 0.05)
    u  <- t(EXP)[c(1,2),]
    u  <- data.frame(PC = "PC", gene = "gene", u)
    for(a in 1:length(percent_variance)){
      o   <- t(EXP)[rownames(pca$rotation)[order(pca$rotation[,a],decreasing = T)[1:nn]],]
      ord <- hclust( dist(o, method = "euclidean"), method = "ward.D" )$order
      o   <- o[rownames(o)[ord],]
      u   <- rbind(u,data.frame(PC = paste0("PC",a),gene = row.names(o), o))
    }
    u  <- u[c(-1,-2),]
    uu <- melt(data = u, id.vars = c("PC","gene"))
    
    pdf(file=paste0(outdir,"/",zz,"/pca_exp/",xx,".heatmap.pdf"),width=16,height=16)
    p <- ggplot(data = uu, aes(x = variable, y = gene))
    p <- p + geom_tile(aes(fill = value))
    p <- p + scale_fill_gradient2(low="blue", high="red")
    p <- p + facet_wrap(vars(PC), scales = "free_y")
    p <- p + theme_bw()
    p <- p + theme(
      axis.text.x=element_text(angle=45, vjust=0, hjust=0,size=10),
      axis.text.y=element_blank(), 
      strip.text=element_text(size=12))
    print(p)
    dev.off()
    
    
    # for(a in 1:length(percent_variance)){
    #   for(b in 1:nrow(EXP)){
    #     cat(rownames(EXP)[b],",PC",a," : ",sum((EXP[b,]-pca$center)/pca$scale * pca$rotation[,a]),"\n")
    #   }
    # }
    
    for(i in 1:3){
      for(j in (i+1):4){
        for(label_group in colnames(GRP)[-1]){
          
          legend       <- as.matrix(unique(GRP[,label_group]))
          class_count  <- as.matrix(table(GRP[,label_group]))
          class_color  <- mycol[1:(length(class_count))]
          class_pch    <- mypch[1:(length(class_count))]
          class        <- data.frame(count=class_count,color=as.character(class_color),pch=class_pch)
          col          <- as.character(class[GRP[rownames(pca$x),][,label_group],]$color)
          pch          <- class[GRP[rownames(pca$x),][,label_group],]$pch  
          lcol         <- as.vector(class[legend,]$color)
          lpch         <- as.vector(class[legend,]$pch)
          
          pdf(file=paste0(outdir,"/",zz,"/pca_2d/","/",xx,".pca_",i,"-",j,".pdf"),width=8,height=8)
          mex<-0.2*abs(max(pca$x[,i])-min(pca$x[,i]))
          mey<-0.2*abs(max(pca$x[,j])-min(pca$x[,j]))
          plot(x = pca$x[,i], y = pca$x[,j], 
               xlim=c(min(pca$x[,i]) - mex, max(pca$x[,i]) + mex),
               ylim=c(min(pca$x[,j]) - mey, max(pca$x[,j]) + mey),
               xlab=paste0("PC",i," : ",round(percent_variance[i],2),"%"),
               ylab=paste0("PC",j," : ",round(percent_variance[j],2),"%"),
               main="PCA",cex=1.2,las=1,pch=pch,col=col)
          pointLabel(x=pca$x[,i],y=pca$x[,j],labels=paste("\n   ",rownames(pca$x),"    \n",sep=""),cex=0.7,col=col)
          legend("topright",legend=legend,col=lcol,pch=lpch,pt.bg=paste(lcol,"FF",sep=""))
          dev.off()
        }
      }
    }
  }
}
