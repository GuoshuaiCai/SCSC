# SCSC for <ins>S</ins>ingle <ins>C</ins>ell <ins>S</ins>equencing Data <ins>C</ins>orrelation Analysis

## Author: Guoshuai Cai


### Platform/data Compatibility
The proposed package is versatile and is logically adaptable to single-cell multi-omic measures. It can be used for assessing the gene-gene co-exression and genetic feature-gene expression correlation. Four strategies (classical, non-zero, dropout-weighted, imputation) were enabled.


### Core tool dependencies:

zinbwave, https://github.com/drisso/zinbwave

DESeq2, https://bioconductor.org/packages/release/bioc/html/DESeq2.html

wcorr, https://cran.r-project.org/web/packages/wCorr/index.html

Other tools for weight calculation (scLink) or data imputation (SAVER, MAGIC, knn-smoothing) are also implemented.

Others: Tool dependencies of all above tools.

## Example
```{r}
load("SCSC_func.R")
load("data_demo.RData")
```

### weight calculation (for weighted strategy)
```{r}

w.exp.do<-weight.calc(exp.count,type="dropout",nc=4)
w.exp.sl<-weight.calc(exp.count,type="sclink",nc=4)
```

### normalization
```{r}
exp.data<-normalize.quant(exp.count)
```

### Data imputation (for imputation strategy), data will be imputated and normalized
```{r}
exp.saver<-data.impute(exp.count, method="saver", nc=4)
exp.saver.res<-exp.saver$res
exp.sv<-exp.saver$data

exp.mg<-data.impute(exp.count, method="magic")$data

exp.ks<-data.impute(exp.count, method="knn-smoothing", k=15)$data

exp.zw<-data.impute(exp.count, method="zinbwave", k=2, nc=4)$data
```

### corrlation calculation
```{r}
#all data (classical)
wcorr.eq<-wcorr.calc.allpairs(exp.data,exp.data,method="pearson",mode="equal")

#non-zero data
wcorr.nz<-wcorr.calc.allpairs(exp.data,exp.data,method="pearson",mode="nonzero")

#dropout-weighted, weights are estiamted by zinbwave
wcorr.dw<-wcorr.calc.allpairs(exp.data,exp.data,w.exp.do,w.exp.do,method="pearson",mode="weighted",alpha=3)

#dropout-weghted, weights are estiamted by scLink
wcorr.sl<-wcorr.calc.allpairs(exp.data,exp.data,w.exp.sl,w.exp.sl,method="pearson",mode="weighted",alpha=1)

#imputed data, by SAVER
wcorr.sv<-wcorr.calc.allpairs(exp.sv,exp.sv,method="pearson",mode="imputed.saver",
                              exp.saver.res=exp.saver.res, gen.saver.res=exp.saver.res)

#imputed data, by MAGIC
wcorr.mg<-wcorr.calc.allpairs(exp.mg,exp.mg,method="pearson",mode="imputed")

#imputed data, by knn-smoothing
wcorr.ks<-wcorr.calc.allpairs(exp.ks,exp.ks,method="pearson",mode="imputed")

#imputed data, by zinbwave
wcorr.zw<-wcorr.calc.allpairs(exp.zw,exp.zw,method="pearson",mode="imputed")
```

### plot
```{r}
wcorr.plot(exp.data, exp.data,w.exp.do,w.exp.do,gen.feature="Gene3599",exp.feature="Gene3408",
           alpha=2,color="red",breaks=100)
```
<img src="https://github.com/GuoshuaiCai/SCSC/blob/main/plot_demo.png" width="500" height="500">
