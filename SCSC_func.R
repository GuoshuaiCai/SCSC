weight.calc<-function(count, K=2, epsilon=1000, nc=1, celltype=NULL,type="dropout", plot=FALSE){

  if(is.null(celltype)){
    celltype<-rep("cell1",ncol(count))
    covar<-1
  }else{
    covar<-"cluster"
  }

  coldata<-data.frame(cluster=factor(celltype))
  rownames(coldata) <- colnames(count)

  if(type=="dropout"){
    
    library("zinbwave")
    library("DESeq2")

    #zinb.wave
    exp.obj <- DESeqDataSetFromMatrix(countData = count, colData=coldata, design=as.formula(paste0("~",covar)))
    zinb.fit <- zinbFit(exp.obj, K=K, epsilon=epsilon,BPPARAM=BiocParallel::MulticoreParam(nc))
    zinb.wave <- zinbwave(exp.obj, fitted_model = zinb.fit, K = 2,
                             epsilon=1000,
                             observationalWeights = TRUE)
    w<- assay(zinb.wave, "weights")
  }

  if(type=="precision"){
    
    library("limma")

    d0                <- DGEList(count)
    d0                <- calcNormFactors(d0)
    design            <- model.matrix(as.formula(paste0("~",covar)), data=coldata)
    voom.fit         <- voom(d0, design, plot = plot,lib.size = colSums(count))
    w<-voom.fit$weights
  }

  if(type=="sclink"){

    library("scImpute")
    source("scLink.R")
    library("parallel")

    totalCounts_by_cell = colSums(count)
    totalCounts_by_cell[totalCounts_by_cell == 0] = 1
    count_lnorm = sweep(count, MARGIN = 2, 10^6/totalCounts_by_cell, FUN = "*")
    count_lnorm = log10(count_lnorm + 1.01)
    parslist<-get_mix_parameters(count = count_lnorm, point = log10(1.01), ncores = nc)
    w<-matrix(NA, nrow=nrow(count_lnorm),ncol=ncol(count_lnorm))
    for(i in 1:nrow(w)){
      if(is.na(parslist[[i]][1])){
        parslist[[i]][1]<-0
      }
      w[i,] = calculate_weight(count_lnorm[i,], parslist[[i]])[,2]
    }
  }

  return(w)
}

normalize.quant<-function(count, method="totalcount", scale.factor=NULL){
  if(is.null(scale.factor)){
    scale.factor<-median(colSums(count))
  }
  if(method=="totalcount"){
    data<-log2(t(t(count)/colSums(count)*scale.factor)+1)
  }
  return(data)
}

data.impute<-function(count, method="saver", nc=1, k=15){
  if(method=="saver"){
    library("SAVER")
    saver.res <- saver(exp.count, ncores = nc)
    exp.count.saver<-saver.res$estimate
    exp.data <- normalize.quant(exp.count.saver)
    return(list("data"=exp.data,"res"=saver.res))
  }
  if(method=="knn-smoothing"){
    source("knn_smooth.R")
    exp.count.knn <- knn_smoothing(exp.count,k=k)
    exp.data <- normalize.quant(exp.count.knn)
    return(list("data"=exp.data))
  }
  if(method=="magic"){
    library("Rmagic")
    exp.data <- t(magic(t(exp.data))$result)
    return(list("data"=exp.data))
  }
}

wcorr.calc<-function(dat.exp, dat.gen, w.exp, w.gen, alpha=3, method="pearson", celltype=NULL, type="GT", mode="weighted", exp.saver.res=NULL, gen.saver.res=NULL){

  library("wCorr")

  #match matrix
  if(type=="GT"){
    gene.gen<-sub(".+_","",rownames(dat.gen))
  }
  if(type=="T"){
    gene.gen<-rownames(dat.gen)
  }

  in.gen<-which(gene.gen %in% rownames(dat.exp))
  dat.gen<-dat.gen[in.gen,,drop=FALSE]

  match.exp<-match(gene.gen,rownames(dat.exp))
  dat.exp<-dat.exp[match.exp,,drop=FALSE]


  if(mode=="equal"){
    w.exp<-matrix(1,dim(dat.exp)[1],dim(dat.exp)[2])
    w.gen<-matrix(1,dim(dat.gen)[1],dim(dat.gen)[2])
  }
  if(mode=="nonzero"){
    w.exp<-(dat.exp>0)
    mode(w.exp)<-"numeric"
    w.gen<-(dat.gen>0)
    mode(w.gen)<-"numeric"
  }
  if(length(grep("imputed",mode))>0){
    w.exp<-matrix(1,dim(dat.exp)[1],dim(dat.exp)[2])
    w.gen<-matrix(1,dim(dat.gen)[1],dim(dat.gen)[2])
  }
  
  if(mode=="imputed.saver"){  
    gen.saver.res.estimate<-gen.saver.res$estimate[in.gen,,drop=FALSE]
    gen.saver.res.se<-gen.saver.res$se[in.gen,,drop=FALSE]
    exp.saver.res.estimate<-exp.saver.res$estimate[match.exp,,drop=FALSE]
    exp.saver.res.se<-exp.saver.res$se[match.exp,,drop=FALSE]
  }
  
  w.gen<-w.gen[in.gen,,drop=FALSE]
  w.exp<-w.exp[match.exp,,drop=FALSE]

  w2<-(w.gen*w.exp)^alpha

  if(is.null(celltype)){
    celltype<-factor(rep("cell1",ncol(dat.exp)))
  }else{
    celltype<-factor(celltype)
  }
  celltype.level<-levels(celltype)

  out<-list()

  for(l in 1:length(celltype.level)){
    celltype.l<-celltype.level[l]
    pos<-which(celltype==celltype.l)

    out.m<-matrix(NA, nrow(dat.exp),2)

    for(i in 1:nrow(dat.exp)){
      exp.i<-dat.exp[i,]
      gen.i<-dat.gen[i,]
      w2.i<-w2[i,]

      if(mode=="imputed.saver"){
        adj.gen.i <- sqrt(var(gen.saver.res.estimate[i, ], na.rm = TRUE)/
                    (var(gen.saver.res.estimate[i, ], na.rm = TRUE) +
                    mean(gen.saver.res.se[i, ]^2, na.rm = TRUE)))
        adj.exp.i <- sqrt(var(exp.saver.res.estimate[i, ], na.rm = TRUE)/
                    (var(exp.saver.res.estimate[i, ], na.rm = TRUE) +
                    mean(exp.saver.res.se[i, ]^2, na.rm = TRUE)))
        adj.i <- adj.gen.i * adj.exp.i
      }else{
        adj.i <- 1
      }

      out.m[i,1]<-weightedCorr(exp.i[pos],gen.i[pos],weight=w2.i[pos], method=method)*adj.i
      out.m[i,2]<-length(which(w2.i[pos]>0.5))
      colnames(out.m)<-paste(c("Rho","n"),celltype.l,sep="_")
    }
    out[[celltype.l]]<-out.m

  }

  out<-do.call("cbind",out)
  out<-cbind(rownames(dat.gen),rownames(dat.exp),out)
  colnames(out)[1:2]<-c("gen.peak","gene")

  return(out)
}

wcorr.calc.allpairs<-function(dat.exp, dat.gen, w.exp, w.gen, alpha=3, method="pearson", celltype=NULL, type="GT", mode="weighted", exp.saver.res=NULL, gen.saver.res=NULL){

  library("wCorr")

  if(type=="GT"){
    gene.gen<-sub(".+_","",rownames(dat.gen))
  }
  if(type=="T"){
    gene.gen<-rownames(dat.gen)
  }

  if(is.null(celltype)){
    celltype<-factor(rep("cell1",ncol(dat.exp)))
  }else{
    celltype<-factor(celltype)
  }
  celltype.level<-levels(celltype)

  if(mode=="equal"){
    w.exp<-matrix(1,dim(dat.exp)[1],dim(dat.exp)[2])
    w.gen<-matrix(1,dim(dat.gen)[1],dim(dat.gen)[2])
  }

  if(mode=="nonzero"){
    w.exp<-(dat.exp>0)
    mode(w.exp)<-"numeric"
    w.gen<-(dat.gen>0)
    mode(w.gen)<-"numeric"
  }

  if(length(grep("imputed",mode))>0){
    w.exp<-matrix(1,dim(dat.exp)[1],dim(dat.exp)[2])
    w.gen<-matrix(1,dim(dat.gen)[1],dim(dat.gen)[2])
  }
  
  if(mode=="imputed.saver"){  
    gen.saver.res.estimate<-gen.saver.res$estimate
    gen.saver.res.se<-gen.saver.res$se
    exp.saver.res.estimate<-exp.saver.res$estimate
    exp.saver.res.se<-exp.saver.res$se
  }

  out<-list()

  for(l in 1:length(celltype.level)){
    celltype.l<-celltype.level[l]
    pos<-which(celltype==celltype.l)

    out.m<-matrix(NA, nrow(dat.exp),nrow(dat.gen))
    rownames(out.m)<-rownames(dat.exp)
    colnames(out.m)<-rownames(dat.gen)

    for(i in 1:nrow(dat.exp)){
      for(j in 1:nrow(dat.gen)){
        exp.i<-dat.exp[i,]
        gen.j<-dat.gen[j,]
        w2.ij<-(w.exp[i,]*w.gen[j,])^alpha

        if(mode=="imputed.saver"){
          adj.gen.i <- sqrt(var(gen.saver.res.estimate[j, ], na.rm = TRUE)/
                    (var(gen.saver.res.estimate[j, ], na.rm = TRUE) +
                    mean(gen.saver.res.se[j, ]^2, na.rm = TRUE)))
          adj.exp.i <- sqrt(var(exp.saver.res.estimate[i, ], na.rm = TRUE)/
                    (var(exp.saver.res.estimate[i, ], na.rm = TRUE) +
                    mean(exp.saver.res.se[i, ]^2, na.rm = TRUE)))
          adj.i <- adj.gen.i * adj.exp.i
        }else{
          adj.i <- 1
        }

        out.m[i,j]<-weightedCorr(exp.i[pos],gen.j[pos],weight=w2.ij[pos], method=method)*adj.i
      }
    }
    out[[celltype.l]]<-out.m
  }

  return(out)
}


wcorr.plot<-function(dat.exp, dat.gen, w.exp, w.gen, alpha=3, celltype=NULL, cell=NULL, exp.feature=NULL, gen.feature=NULL, color="red", pch=20, cex=1, breaks=50){

  #match matrix

  if(is.null(exp.feature)|is.null(gen.feature)){
    stop("the features need to be specified")
  }

  if(is.null(cell)){
    pos<-1:ncol(dat.exp)
  }else{
    if(is.null(celltype)){
      stop("the celltype information need to be provided")
    }
    pos<-which(celltype==cell)
  }

  row.gen<-which(rownames(dat.gen)==gen.feature)
  gen.i<-dat.gen[row.gen,pos]
  w.gen.i<-w.gen[row.gen,pos]

  row.exp<-which(rownames(dat.exp)==exp.feature)
  exp.i<-dat.exp[row.exp,pos]
  w.exp.i<-w.exp[row.exp,pos]

  w2.i<-(w.gen.i*w.exp.i)^alpha

  cat.i<-as.numeric(cut(w2.i,breaks = breaks))

  color.i<-colorRampPalette(c("grey",color))(breaks)[cat.i]

  plot(gen.i,exp.i, pch=pch, cex=cex, col=color.i,xlab="Genetic feature abundance",ylab="Expression", cex.lab=1.5, cex.axis=1.3)
  abline(lm(exp.i~gen.i,weight=w2.i),col = color,lwd=2)
}